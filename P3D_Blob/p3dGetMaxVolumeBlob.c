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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dBlob.h"
#include "p3dTime.h"

#include "Common/p3dConnectedComponentsLabeling.h"
#include "Common/p3dUtils.h"

int p3dGetMaxVolumeBlob3D(
        unsigned char* in_rev,
        unsigned char* out_rev,
        const int dimx,
        const int dimy,
        const int dimz,
        int conn,
        int (*wr_log)(const char*, ...)
        ) {
    unsigned short* lbl_rev;

    unsigned int* volumes;
    unsigned short num_el, lbl;

    int i, j, k;

    int lbl_max;
    int vol_max;    


    /*char auth_code;

    //
    // Authenticate:
    //
    //auth_code = authenticate("p3dGetMaxVolumeBlob3D");
    //if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Extracting blob having maximum volume...");
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
    }

    // Initialize output cloning input:
    memcpy(out_rev, in_rev, dimx * dimy * dimz * sizeof (unsigned char));

    // Allocate memory for labels:
    P3D_TRY( lbl_rev = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));
   
    // Perform connected component labeling:
    P3D_TRY( p3dConnectedComponentsLabeling_ushort(in_rev, lbl_rev, &num_el, &volumes, NULL, dimx,
             dimy, dimz, conn, P3D_FALSE, P3D_FALSE));   


    lbl_max = 0;
    vol_max = 0;

    // For each connected component labeled:
    for (lbl = 0; lbl < num_el; lbl++) {
        if (volumes[lbl] > vol_max) {
            lbl_max = lbl;
            vol_max = volumes[lbl];
        }
    }

    // Remove connected component:
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                // Labels start from 3:
                if (lbl_rev[ I(i, j, k, dimx, dimy) ] != (lbl_max + 3)) {
                    out_rev[ I(i, j, k, dimx, dimy) ] = 0;
                }
            }

   // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tPore3D - The volume of extracted blob is %d voxels.", vol_max);
        wr_log("Pore3D - Blob having maximum area extracted successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (lbl_rev != NULL) free(lbl_rev);
    if (volumes != NULL) free(volumes);

    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Free memory:
    if (lbl_rev != NULL) free(lbl_rev);
    if (volumes != NULL) free(volumes);

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Return OK:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}






