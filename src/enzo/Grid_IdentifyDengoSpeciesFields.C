/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE SPECIES FIELDS FOR SIMON GLOVER'S COOLING)
/
/  written by: Britton Smith
/  date:       September, 2007
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);


int grid::IdentifyDengoSpeciesFields(int &H2_1Num,int &H2_2Num,int &H_1Num,int &H_2Num,int &H_m0Num,int &He_1Num,int &He_2Num,int &He_3Num,int &deNum)
{
    H2_1Num = 0;
    H2_2Num = 0;
    H_1Num = 0;
    H_2Num = 0;
    H_m0Num = 0;
    He_1Num = 0;
    He_2Num = 0;
    He_3Num = 0;
    deNum = 0;
    H2_1Num = FindField(H2_1Density, FieldType, NumberOfBaryonFields);
    H2_2Num = FindField(H2_2Density, FieldType, NumberOfBaryonFields);
    H_1Num = FindField(H_1Density, FieldType, NumberOfBaryonFields);
    H_2Num = FindField(H_2Density, FieldType, NumberOfBaryonFields);
    H_m0Num = FindField(H_m0Density, FieldType, NumberOfBaryonFields);
    He_1Num = FindField(He_1Density, FieldType, NumberOfBaryonFields);
    He_2Num = FindField(He_2Density, FieldType, NumberOfBaryonFields);
    He_3Num = FindField(He_3Density, FieldType, NumberOfBaryonFields);
    deNum = FindField(deDensity, FieldType, NumberOfBaryonFields);

    if ((H2_1Num < 0)||(H2_2Num < 0)||(H_1Num < 0)||(H_2Num < 0)||(H_m0Num < 0)||(He_1Num < 0)||(He_2Num < 0)||(He_3Num < 0)||(deNum < 0)
    ){
    
    ENZO_VFAIL("De=%"ISYM", HI=%"ISYM", HII=%"ISYM", HeI=%"ISYM", HeII=%"ISYM", HeIII=%"ISYM"\n", deNum, H_1Num, H_2Num, He_1Num, He_2Num);

    ENZO_VFAIL("Error identifying species for DengoChemistryModel = %"ISYM".\n",
        DengoChemistryModel)

    }

    return SUCCESS;
}
 