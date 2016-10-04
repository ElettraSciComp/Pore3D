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
#include <omp.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dUtils.h"

int p3dSkeletonLabeling(
        unsigned char* in_im, // IN: Input (binary) skeleton
        unsigned char* out_im, // OUT: Labeled skeleton
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...)
        ) {

    coords_list_t list; // FIFO data structure for coords storage
    coords_t coords; // Current coords

    int i, j, k;
    int a, b, c;
    int ct;

    // Dimensions for padded/cropped volumes:
    int a_dimx, a_dimy, a_dimz;
    int a_rad;

    // Padded and cropped temporary input and output:
    unsigned char* tmp_im;

    // Counters:
    int neigh, length;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dSkeletonLabeling");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Labeling the skeleton...");
    }

    // Create temporary input replicate-padded:
    a_rad = 1;

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;


    // Initialize input:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
    P3D_TRY(p3dZeroPadding3D_uchar2uchar(in_im, tmp_im, dimx, dimy, dimz, a_rad));


    // Volume scanning:
#pragma omp parallel for private(a, b, i, j, k, ct, coords, list, neigh, length) 
    for (c = a_rad; c < (a_dimz - a_rad); c++)
        for (b = a_rad; b < (a_dimy - a_rad); b++)
            for (a = a_rad; a < (a_dimx - a_rad); a++) {
                // If we're on a skeleton voxel:
                if (tmp_im[ I(a, b, c, a_dimx, a_dimy) ] == OBJECT) {
                    // Check Neighborhood:
                    ct = countNeighbors(tmp_im, a_dimx, a_dimy, a_dimz, a, b, c);

                    // Is an end point?
                    if (ct == 1) {
                        // Initialize variables:
                        coords_list_init(&list);

                        // End point found:
                        coords.x = a;
                        coords.y = b;
                        coords.z = c;

                        // Start tracking length:
                        length = 0;

                        do {
                            // Push the voxel into list:
                            coords_list_push(&list, coords);

                            // Get coords:
                            k = coords.z;
                            j = coords.y;
                            i = coords.x;

                            // Temporary delete voxel
                            tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = BACKGROUND;

                            // Determine number of neighbors:
                            neigh = findNeighbor(tmp_im, a_dimx, a_dimy, a_dimz, i, j, k, &coords);

                            // Increment counter of branch length:
                            length++;

                        } while (neigh == 1);

                        // At this point, we're on last voxel of node-to-end branch (ct > 1) or we 
                        // completely scanned a end-to-end branch (ct == 0). In the first case we 
                        // need to take care whether last voxel is a simple point or not:
                        if (neigh > 1) {
                            // Assign last voxel to NODE label:
                            tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = NODE_LABEL;

                            // Remove coordinates from the list:
                            coords = coords_list_pop(&list);

                            // Assign NODETOEND labels:
                            while (coords_list_isempty(list) == P3D_FALSE) {
                                // Get coordinates from the list:
                                coords = coords_list_pop(&list);

                                k = coords.z;
                                j = coords.y;
                                i = coords.x;

                                // Re-assign voxel:
                                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = NODETOEND_LABEL; // NODE-TO-END Label
                            }

                            // Set the endpoint as END_LABEL:
                            tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = END_LABEL;
                        } else {
                            tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = END_LABEL; // END-POINT Label

                            // Remove coordinates from the list:
                            coords = coords_list_pop(&list);

                            // Restore ENDTOEND branch:
                            while (coords_list_isempty(list) == P3D_FALSE) {
                                // Get coordinates from the list:
                                coords = coords_list_pop(&list);

                                k = coords.z;
                                j = coords.y;
                                i = coords.x;

                                // Re-assign voxel:
                                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = ENDTOEND_LABEL;
                            }

                            // Set the endpoint as END_LABEL:
                            tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = END_LABEL;
                        }
                    }// end of cycle on each endpoint				

                        // Is an isolated voxel?
                    else if (ct == 0) {
                        // Isolated voxel removal:
                        tmp_im[ I(a, b, c, a_dimx, a_dimy) ] = ISOLATED_LABEL;
                    }
                }
            }

    // At this point we need to label NODETONODE branches and NODES:	
#pragma omp parallel for private(i, j, ct) 
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                if (tmp_im[ I(i, j, k, a_dimx, a_dimy) ] == OBJECT) {
                    // Check Neighborhood:
                    ct = countNeighbors(tmp_im, a_dimx, a_dimy, a_dimz, i, j, k);

                    // Is a node voxel?
                    if (ct > 2) {
                        // Node point found:
                        tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = NODE_LABEL;
                    }
                    else {
                        // Node-to-Node voxel found:
                        tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = NODETONODE_LABEL;
                    }
                }
            }


    // Crop output:
    P3D_TRY(p3dCrop3D_uchar2uchar(tmp_im, out_im, a_dimx, a_dimy, a_dimz, a_rad));

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Skeleton labeled successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:	
    if (tmp_im != NULL) free(tmp_im);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:	
    if (tmp_im != NULL) free(tmp_im);

    // Return error:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}



