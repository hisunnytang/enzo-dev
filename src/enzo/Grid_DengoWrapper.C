/***********************************************************************
  /
  /  GRID CLASS (WRAP THE DENGO CHEMISTRY SOLVER)
  /
  /  written by: Sunny Tang
  /  date:       April, 2020
  /  modified1:
  /
  /  PURPOSE: Solve chemistry and cooling with dengo.
  /
  /  RETURNS:
  /    SUCCESS or FAIL
  /
 ************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
        float *TemperatureUnits, float *TimeUnits,
        float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

#ifdef USE_DENGO
int grid::DengoWrapper()
{


    if (ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

    LCAPERF_START("grid_DengoWrapper");

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    int H2_1Num;
    int H2_2Num;
    int H_1Num;
    int H_2Num;
    int H_m0Num;
    int He_1Num;
    int He_2Num;
    int He_3Num;
    int deNum;

    double dt_cool = dtFixed;

    // Radiative Transfer is not supported by dengo as of now...
    /*
#ifdef TRANSFER
dt_cool = (grackle_data->radiative_transfer_intermediate_step == TRUE) ? dtPhoton : dtFixed;
#endif
     */

    /* Compute the size of the fields. */

    int i;
    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];

    Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
    g_grid_dimension = new Eint32[GridRank];
    g_grid_start = new Eint32[GridRank];
    g_grid_end = new Eint32[GridRank];
    for (i = 0; i < GridRank; i++) {
        g_grid_dimension[i] = (Eint32) GridDimension[i];
        g_grid_start[i] = (Eint32) GridStartIndex[i];
        g_grid_end[i] = (Eint32) GridEndIndex[i];
    }

    /* Find fields: density, total energy, velocity1-3. */
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                Vel3Num, TENum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

    /* Find Multi-species fields. */
    // Define the index here by dengo
    H2_1Num = 0;
    H2_2Num = 0;
    H_1Num = 0;
    H_2Num = 0;
    H_m0Num = 0;
    He_1Num = 0;
    He_2Num = 0;
    He_3Num = 0;
    deNum = 0;

    // TODO: fixed the identifyspeciesfields
    /* Find Multi-species fields. */
    if (IdentifyDengoSpeciesFields(H2_1Num,H2_2Num,H_1Num,H_2Num,H_m0Num,He_1Num,He_2Num,He_3Num,deNum) == FAIL)
    {
        ENZO_FAIL("Error in grid->IdentifyDengoSpeciesFields.\n");
    }
    /* Get easy to handle pointers for each variable. */

    float *density     = BaryonField[DensNum];
    float *totalenergy = BaryonField[TENum];
    float *gasenergy   = BaryonField[GENum];
    float *velocity1   = BaryonField[Vel1Num];
    float *velocity2   = BaryonField[Vel2Num];
    float *velocity3   = BaryonField[Vel3Num];

    float *volumetric_heating_rate = NULL;
    float *specific_heating_rate   = NULL;

    /* Compute the cooling time. */

    FLOAT a = 1.0, dadt;
    float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
          VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time);
    if (ComovingCoordinates) {
        CosmologyComputeExpansionFactor(Time+0.5*dt_cool, &a, &dadt);

        aUnits = 1.0/(1.0 + InitialRedshift);
    } else if (RadiationFieldRedshift > -1){
        a        = 1.0 / (1.0 + RadiationFieldRedshift);
        aUnits   = 1.0;
    }
    float afloat = float(a);

    /* Update units. */

    code_units dengo_units;
    dengo_units.comoving_coordinates = (Eint32) ComovingCoordinates;
    dengo_units.density_units        = (double) DensityUnits;
    dengo_units.length_units         = (double) LengthUnits;
    dengo_units.time_units           = (double) TimeUnits;
    dengo_units.velocity_units       = (double) VelocityUnits;
    dengo_units.a_units              = (double) aUnits;
    dengo_units.a_value              = (double) a;

    // no metal cooling yet in dengo...
    /*
    // Metal cooling codes. 

    int MetalNum = 0, SNColourNum = 0;
    int MetalFieldPresent = FALSE;

    // First see if there's a metal field (so we can conserve species in
    // the solver)
    MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
    SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
    MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

    // Double check if there's a metal field when we have metal cooling
    if (MetalCooling && MetalFieldPresent == FALSE) {
    if (debug)
    fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
    }

    // If both metal fields (Pop I/II and III) exist, create a field
    that contains their sum 

    float *MetalPointer = NULL;
    float *TotalMetals = NULL;

    if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
    TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
    } // ENDIF both metal types
    else {
    if (MetalNum != -1)
    MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
    MetalPointer = BaryonField[SNColourNum];
    } // ENDELSE both metal types
     */

    int temp_thermal = FALSE;
    float *thermal_energy;
    if ( UseMHD ){
        iBx = FindField(Bfield1, FieldType, NumberOfBaryonFields);
        iBy = FindField(Bfield2, FieldType, NumberOfBaryonFields);
        iBz = FindField(Bfield3, FieldType, NumberOfBaryonFields);  
    }

    if (HydroMethod==Zeus_Hydro) {
        thermal_energy = BaryonField[TENum];
    }
    else if (DualEnergyFormalism) {
        thermal_energy = BaryonField[GENum];
    }
    else {
        temp_thermal = TRUE;
        thermal_energy = new float[size];
        for (i = 0; i < size; i++) {
            thermal_energy[i] = BaryonField[TENum][i] - 
                0.5 * POW(BaryonField[Vel1Num][i], 2.0);
            if(GridRank > 1)
                thermal_energy[i] -= 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
            if(GridRank > 2)
                thermal_energy[i] -= 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

            if( UseMHD ) {
                thermal_energy[i] -= 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                        POW(BaryonField[iBy][i], 2.0) + 
                        POW(BaryonField[iBz][i], 2.0)) / 
                    BaryonField[DensNum][i];
            }
        } // for (int i = 0; i < size; i++)
    }

    //
    // Put code here to assign fields to volumetric or specific
    // heating rate pointers
    //


    /* set up the my_fields */
    dengo_field_data my_fields;
    // instead of sending the grid start and end idx
    // we offset the pointer already by the 
    // pointer to these scalars are sent in
    int cell_counts = 1;
    for (i = 0; i < GridRank; i++) {
        cell_counts *= g_grid_end[i] - g_grid_start[i];
    }
    my_fields.ncells = (unsigned long int) cell_counts;
    // TODO: make readable from the input.enzo file
    my_fields.dengo_data_file = dengo_data_file;

    Eint32 xdim = g_grid_dimension[0];
    Eint32 ydim = g_grid_dimension[1];
    Eint32 zdim = g_grid_dimension[2];

    Eint32 xs = g_grid_start[0];
    Eint32 ys = g_grid_start[1];
    Eint32 zs = g_grid_start[2];

    my_fields.grid_start     = g_grid_start;
    my_fields.grid_end       = g_grid_end;
    my_fields.grid_dimension = g_grid_dimension;

    my_fields.reltol = dengo_reltol;
    my_fields.floor_value = 1.0e-30;

    // I guess the offset are like arr[k][j][i]
    Eint32 offset = 0;// xdim*ydim* zs + ydim* ys + xs;

    // DENGO should spit out the required arrays in here!
    /* now add in the baryon fields */

    my_fields.density         = ((double *)density + offset);
    my_fields.ge_density      = thermal_energy; // 0 offset by construct

    my_fields.density = ((double*)density + offset);
    my_fields.H2_1_density = ((double *) BaryonField[H2_1Num] + offset);
    my_fields.H2_2_density = ((double *) BaryonField[H2_2Num] + offset);
    my_fields.H_1_density = ((double *) BaryonField[H_1Num] + offset);
    my_fields.H_2_density = ((double *) BaryonField[H_2Num] + offset);
    my_fields.H_m0_density = ((double *) BaryonField[H_m0Num] + offset);
    my_fields.He_1_density = ((double *) BaryonField[He_1Num] + offset);
    my_fields.He_2_density = ((double *) BaryonField[He_2Num] + offset);
    my_fields.He_3_density = ((double *) BaryonField[He_3Num] + offset);
    my_fields.de_density = ((double *) BaryonField[deNum] + offset);

    // inlude metal, heating, photoionization
    // metal, heating rate, photoionization rate are not supported yet...
    //my_fields.metal_density   = MetalPointer;
    //my_fields.volumetric_heating_rate  = volumetric_heating_rate;
    //my_fields.specific_heating_rate    = specific_heating_rate;

    /* WOULD BE A GREAT NEXT STEP! no radiative transfer yet for dengo.... */
    /*
#ifdef TRANSFER
    // Find RT fields
    int kphHINum, kphHeINum, kphHeIINum, kdissH2INum,
    gammaNum, kphHMNum, kdissH2IINum;

    IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
    kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

    // unit conversion from Enzo RT units to CGS
    float rtunits = erg_eV / TimeUnits;

    if( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate   = BaryonField[kphHINum];

    if (RadiativeTransferHydrogenOnly == FALSE){
    my_fields.RT_HeI_ionization_rate  = BaryonField[kphHeINum];
    my_fields.RT_HeII_ionization_rate = BaryonField[kphHeIINum];
    }

    if (MultiSpecies > 1)
    my_fields.RT_H2_dissociation_rate = BaryonField[kdissH2INum];

    // need to convert to CGS units
    for( i = 0; i < size; i++) BaryonField[gammaNum][i] *= rtunits;

    my_fields.RT_heating_rate = BaryonField[gammaNum];


    }
#endif // TRANSFER
     */

    /* Call the chemistry solver. */

    if (primordial_solve_chemistry_enzo(&dengo_units, &my_fields, (double) dt_cool) == FAIL){
        fprintf(stderr, "Error in dengo solve_chemistry.\n");
        return FAIL;
    }

    if (HydroMethod != Zeus_Hydro) {
        for (i = 0; i < size; i++) {
            BaryonField[TENum][i] = thermal_energy[i] +
                0.5 * POW(BaryonField[Vel1Num][i], 2.0);
            if(GridRank > 1)
                BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
            if(GridRank > 2)
                BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

            if( UseMHD ) {
                BaryonField[TENum][i] += 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                        POW(BaryonField[iBy][i], 2.0) + 
                        POW(BaryonField[iBz][i], 2.0)) / 
                    BaryonField[DensNum][i];
            }

        } // for (int i = 0; i < size; i++)
    } // if (HydroMethod != Zeus_Hydro)

    if (temp_thermal == TRUE) {
        delete [] thermal_energy;
    }


    delete [] g_grid_dimension;
    delete [] g_grid_start;
    delete [] g_grid_end;

    LCAPERF_STOP("grid_DengoWrapper");

    return SUCCESS;
}
#endif
