/***************************************************************************/
/* (C) 2016 Elettra - Sincrotrone Trieste S.C.p.A.. All rights reserved.   */
/*                                                                         */
/*                                                                         */
/* This file is part of Pore3D, a software library for quantitative        */
/* analysis of 3D (volume) images.                                         */
/*                                                                         */
/* Pore3D is free software: you can redistribute it and/or modify it       */
/* under the terms of the GNU General Public License as published by the   */
/* Free Software Foundation, either version 3 of the License, or (at your  */
/* option) any later version.                                              */
/*                                                                         */
/* Pore3D is distributed in the hope that it will be useful, but WITHOUT   */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   */
/* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    */
/* for more details.                                                       */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with Pore3D. If not, see <http://www.gnu.org/licenses/>.          */
/*                                                                         */
/***************************************************************************/

//
// Author: Francesco Brun
// Last modified: Sept, 28th 2016
//

#include <omp.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <stdio.h>

#include "p3dBlob.h"
#include "p3dTime.h"

#include "Common/p3dUIntList.h"
#include "Common/p3dDoubleList.h"

#define BACKGROUND 0

int p3dREVEstimation(
        unsigned char* in_rev, // IN: binary volume
        double** porosity, // OUT: array of porosity for the related cube side
        unsigned int** cubeEdges, // OUT: array of length of the side of the REV cube
        unsigned int* numel, // OUT: length of the arrays
        const int dimx,
        const int dimy,
        const int dimz,
        const int stepsize, // IN: stepsize for porosity
        const int centerx, // IN: X coordinate for the center of REV (-1 for 
        //     automatic determination of the center)
        const int centery, // IN: Y coordinate for the center of REV
        const int centerz, // IN: Z coordinate for the center of REV
        int (*wr_log)(const char*, ...)
        ) {
    
    int i, j, k;
    int porosityct;
    int step;
    int cen_x, cen_y, cen_z;

    double_list_t porosity_list;
    uint_list_t edges_list;
    int arrayct;



    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dREVEstimation");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Initialize lists:
    double_list_init(&porosity_list);
    uint_list_init(&edges_list);

    // Initialize step:
    step = stepsize;
    arrayct = 0;

    // Determine the center of the REV:
    if (centerx == -1) {
        // Automatic determination of the center:
        cen_x = dimx / 2;
        cen_y = dimy / 2;
        cen_z = dimz / 2;
    } else {
        // Assign values from input parameters:
        cen_x = (int) centerx;
        cen_y = (int) centery;
        cen_z = (int) centerz;
    }
    
    
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Estimating the REV on a porosity-basis...");
        wr_log("\tCenter: [ %d, %d, %d].", cen_x, cen_y, cen_z);
        wr_log("\tStepsize: %d [voxels].", stepsize);
    }

    // Loop until boundaries reached:
    while (((cen_x - step / 2) > 0) &&
            ((cen_y - step / 2) > 0) &&
            ((cen_z - step / 2) > 0) &&
            ((cen_z + step / 2) < dimz) &&
            ((cen_y + step / 2) < dimy) &&
            ((cen_x + step / 2) < dimx)
            ) {
        porosityct = 0;

        for (k = (cen_z - step / 2); k < (cen_z + step / 2); k++)
            for (j = (cen_y - step / 2); j < (cen_y + step / 2); j++)
                for (i = (cen_x - step / 2); i < (cen_x + step / 2); i++) {
                    if (in_rev[ I(i, j, k, dimx, dimy) ] == 0) {
                        porosityct++;
                    }
                }


        // Put porosity in the list:
        double_list_add(&porosity_list, (double) (porosityct / (step * step * step * 1.0)));

        // Put stepsize in the list:
        uint_list_add(&edges_list, (unsigned int) (step + 0));

        // Increment array counter:
        arrayct++;

        // Increase stepsize:
        step = step + stepsize;
    }


    // Return lists converted to array:
    *porosity = double_list_toarray(&porosity_list, arrayct);
    *cubeEdges = uint_list_toarray(&edges_list, arrayct);
    *numel = arrayct;


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - REV estimation performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

