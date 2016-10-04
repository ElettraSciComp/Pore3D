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


int _p3dCheckForObjectVoxel(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        const int dimx,
        const int dimy,
        const int dimz,
        const int i,
        const int j,
        const int k,
        const int size
        ) {
    int x, y, z;

    // Check if there is an OBJECT pixel within box (excluding partial boxes):
    /*if ( ((k + size - 1) < dimz ) && 
             ((j + size - 1) < dimy ) &&
             ((i + size - 1) < dimx ) )
    {
            for (z = (int) (k + size - 1); z >= (int) k; z--)                     
                    for (y = (int) (j + size - 1); y >= (int) j; y--)
                            for (x = (int) (i + size - 1); x >= (int) i; x--)
                            {
                                    if ( in_im[ I(x,y,z,dimx,dimy)] == OBJECT )
                                    {
                                            return 1; // true
                                    }
                            }
    }*/

    // Check if there is an OBJECT pixel within box (including partial boxes):
    for (z = k; ((z < (k + size - 1)) && (z < dimz)); z++)
        for (y = j; ((y < (j + size - 1)) && (y < dimz)); y++)
            for (x = i; ((x < (i + size - 1)) && (x < dimz)); x++) {
                if (in_im[ I(x, y, z, dimx, dimy)] == OBJECT) {
                    return 1; // true
                }
            }

    return 0; // false;
}

int p3dBoxCountingFractalDimension(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        double* fd, // OUT: Box Counting Fractal Dimension
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...)
        ) {

    int i, j, k, ct;
    int size, boxct;

    int min_dim, max_box;

    double sizes_avg, boxes_avg, sizes_sum, boxes_sum;
    double cross, r;

    double* sizes;
    double* boxes;


    // Determine maximum box size in order to have more than one boxes for
    // the greater step size. For instance using a maximum step size of 128
    // for a 300x300x300 volume will result for sure in 4 full box e 5 
    // partial boxes for the 128 step size. This situation affects the 
    // accuracy of results.
    min_dim = MIN(MIN(dimx, dimy), dimz);
    max_box = (int) (log((double) min_dim) / log(2.0) + 0.5) - 1;

     if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tMinimum box side adopted: 2");    
        wr_log("\tMaximum box side adopted: %d. ",(int) ((pow(2, (double) max_box))));      
    }

    sizes = (double*) malloc(max_box * sizeof (double));
    boxes = (double*) malloc(max_box * sizeof (double));

    sizes_avg = 0;
    boxes_avg = 0;

    // For each size:	
    for (ct = 0; ct < max_box; ct++) {
        // Get current size:
        size = (int) (pow((double) 2.0, (double) ct + 1));

        // Reset box counter:
        boxct = 0;

        // Image scanning:  
#pragma omp parallel for private(i, j) reduction (+ : boxct)
        for (k = 0; k < dimz; k = k + size)
            for (j = 0; j < dimy; j = j + size)
                for (i = 0; i < dimx; i = i + size) {
                    // If there's a ON pixel within box increment counter:
                    if (_p3dCheckForObjectVoxel(in_im, (int) dimx, (int) dimy, (int) dimz, i, j, k, size) == 1)
                        boxct++;
                }

        // Set total number of boxes to vectors:
        sizes[ct] = -log((double) size);
        boxes[ct] = log((double) boxct);

        sizes_avg += sizes[ct];
        boxes_avg += boxes[ct];
    }

    sizes_avg /= (double) (max_box);
    boxes_avg /= (double) (max_box);

    sizes_sum = 0;
    boxes_sum = 0;
    cross = 0;

    // Compute sum of squares:
    for (ct = 0; ct < max_box; ct++) {
        cross += (boxes[ct] - boxes_avg)*(sizes[ct] - sizes_avg);
        boxes_sum += (boxes[ct] - boxes_avg)*(boxes[ct] - boxes_avg);
        sizes_sum += (sizes[ct] - sizes_avg)*(sizes[ct] - sizes_avg);
    }

    // Compute correlation coefficient and fractal dimension:
    r = cross / sqrt(boxes_sum * sizes_sum);
    *fd = r * sqrt(boxes_sum) / sqrt(sizes_sum);


    // Release resources:
    if (sizes != NULL) free(sizes);
    if (boxes != NULL) free(boxes);


    // Return OK:
    return P3D_SUCCESS;
}

int p3dTextureAnalysis(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        struct TextureStats* out_stats, // OUT: Texture statistics
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...)
        ) {

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dTextureAnalysis");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing texture analysis...");

    }

    p3dBoxCountingFractalDimension(in_im, &(out_stats->FD), dimx, dimy, dimz, wr_log);
  
    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tFractal dimension (box counting) (FD): %0.3f [-].", out_stats->FD);        
    }

   // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Texture analysis performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }


    // Return OK:
    return P3D_SUCCESS;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}





