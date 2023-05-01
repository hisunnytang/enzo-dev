/***********************************************************************
  /
  /  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
  /
  /  written by: Greg Bryan
  /  date:       April, 1995
  /  modified1:  Elizabeth Harper-Clark, August 2009
  /              added in CoolingModel parameter
  /
  /  PURPOSE:
  /
  /  RETURNS:
  /
 ************************************************************************/

// Compute the cooling time

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "phys_constants.h"

/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */

extern int RadiationFieldRecomputeMetalRates;

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
        float *TemperatureUnits, float *TimeUnits,
        float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);


/* BeginDengo */
int grid::ComputeCoolingTime(float *cooling_time, int CoolingTimeOnly)
{
    if (RadiativeCooling == 0) return SUCCESS;
    if (ProcessorNumber != MyProcessorNumber)
        return SUCCESS;
    int H2_1Num = 0;
    int H2_2Num = 0;
    int H_1Num = 0;
    int H_2Num = 0;
    int H_m0Num = 0;
    int He_1Num = 0;
    int He_2Num = 0;
    int He_3Num = 0;
    int deNum = 0;/* Compute the size of the fields. */
    int i;
    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    /* Find fields: density, total energy, velocity1-3. */
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                Vel3Num, TENum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

    /* Find Multi-species fields. */
    if (IdentifyDengoSpeciesFields(H2_1Num,H2_2Num,H_1Num,H_2Num,H_m0Num,He_1Num,He_2Num,He_3Num,deNum) == FAIL)
    {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }

    /* Find photo-ionization fields */

    int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
    int gammaNum;
    IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
            kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

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
        CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);

        aUnits = 1.0/(1.0 + InitialRedshift);
    } else if (RadiationFieldRedshift > -1){
        a       = 1.0 / (1.0 + RadiationFieldRedshift);
        aUnits  = 1.0;
    }
    float afloat = float(a);

    /* Metal cooling codes. */

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

    /* If both metal fields (Pop I/II and III) exist, create a field
       that contains their sum */

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
    
    // obtain the grid dimension
    Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
    g_grid_dimension = new Eint32[GridRank];
    g_grid_start = new Eint32[GridRank];
    g_grid_end = new Eint32[GridRank];
    for (i = 0; i < GridRank; i++) {
      g_grid_dimension[i] = (Eint32) GridDimension[i];
      g_grid_start[i] = 0;
      g_grid_end[i] = (Eint32) GridDimension[i]-1;
    }

    /* Update units. */
    code_units dengo_units;
    dengo_units.comoving_coordinates = (Eint32) ComovingCoordinates;
    dengo_units.density_units        = (double) DensityUnits;
    dengo_units.length_units         = (double) LengthUnits;
    dengo_units.time_units           = (double) TimeUnits;
    dengo_units.velocity_units       = (double) VelocityUnits;
    dengo_units.a_units              = (double) aUnits;
    dengo_units.a_value              = (double) a;

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

    /* set up the my_fields */
    dengo_field_data my_fields;

    // instead of sending the grid start and end idx
    // we offset the pointer already by the 
    // pointer to these scalars are sent in

    Eint32 cell_counts = 1;
    for (i = 0; i < GridRank; i++) {
      cell_counts *= g_grid_dimension[i];
    }
    my_fields.ncells = (unsigned long int) cell_counts;
    my_fields.dengo_data_file = dengo_data_file;
    
    Eint32 xdim = g_grid_dimension[0];
    Eint32 ydim = g_grid_dimension[1];
    Eint32 zdim = g_grid_dimension[2];

    Eint32 xs = g_grid_start[0];
    Eint32 ys = g_grid_start[1];
    Eint32 zs = g_grid_start[2];
    
    // I guess the offset are like arr[k][j][i]
    Eint32 offset = xdim*ydim* zs + ydim* ys + xs;

    // DENGO should spit out the required arrays in here!
    /* now add in the baryon fields */
    my_fields.density         = ((double *)density + offset);
    my_fields.ge_density      = thermal_energy; // 0 offset by construct
    my_fields.H2_1_density = ((double *) BaryonField[H2_1Num] + offset);
    my_fields.H2_2_density = ((double *) BaryonField[H2_2Num] + offset);
    my_fields.H_1_density = ((double *) BaryonField[H_1Num] + offset);
    my_fields.H_2_density = ((double *) BaryonField[H_2Num] + offset);
    my_fields.H_m0_density = ((double *) BaryonField[H_m0Num] + offset);
    my_fields.He_1_density = ((double *) BaryonField[He_1Num] + offset);
    my_fields.He_2_density = ((double *) BaryonField[He_2Num] + offset);
    my_fields.He_3_density = ((double *) BaryonField[He_3Num] + offset);
    my_fields.de_density = ((double *) BaryonField[deNum] + offset);

    my_fields.grid_start     = g_grid_start;
    my_fields.grid_end       = g_grid_end;
    my_fields.grid_dimension = g_grid_dimension;
    
    my_fields.CoolingTime = cooling_time;
    if (dengo_estimate_cooling_time_enzo(&dengo_units, &my_fields) == FAIL) {
      ENZO_FAIL("Error in Dengo calculate_cooling_time.\n");
    }

    for (i = 0; i < size; i++) {
      cooling_time[i] = fabs(cooling_time[i]);
    }

    if (temp_thermal == TRUE) {
      delete [] thermal_energy;
    }
    
    delete [] TotalMetals;
    delete [] g_grid_dimension;
    delete [] g_grid_start;
    delete [] g_grid_end;

    return SUCCESS;
}
/* EndDengo */
