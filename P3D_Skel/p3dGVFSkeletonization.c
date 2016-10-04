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
#include <omp.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "GVFSkeletonization/p3dGVF.h" 
#include "GVFSkeletonization/p3dCriticalPoints.h"
#include "GVFSkeletonization/p3dHighDivPointList.h"
#include "GVFSkeletonization/p3dGetHighDivPoints.h"
#include "GVFSkeletonization/p3dComputeCoreSkeleton.h"
#include "GVFSkeletonization/p3dComputeHierarchicalSkeleton.h"

#include "Common/p3dFCoordsList.h"
#include "Common/p3dThinning.h"
#include "Common/p3dUtils.h"

int p3dGVFSkeletonization(
        unsigned char* in_im,
        unsigned char* skl_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const double mu,
        const double eps,
        const double thresh,
        const double scale,
        int (*wr_log)(const char*, ...)
        ) {

    int a_rad = 1;
    int a_dimx, a_dimy, a_dimz;
    int i, j, k;

    crit_point_list_t crit_point_list = NULL;
    fcoords_list_t skel_point_list = NULL;
    highDiv_point_list_t highDiv_point_list = NULL;

    int ct = 0;

    double step = STEP_SIZE; // default value
    double close_dist = CLOSENESS; // default value

    fcoords_t tmp_coords;
    int a, b, c;

    // Matrixes:
    unsigned char* tmp_im;
    float* gvf_x;
    float* gvf_y;
    float* gvf_z;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dGVFSkeletonization");
    if (auth_code == '0') goto AUTH_ERROR;*/

    //
    // Initializations:
    //

    double t1 = 0.0;
    double t2 = 0.0;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Peforming GVF skeletonization...");
        wr_log("\tScaling factor: %0.3f.", scale);
        wr_log("\tAdditional branches factor: %0.3f.", thresh);
        wr_log("\tGradient Vector Flow (GVF) parameters: mu = %0.3f eps = %0.5f.", mu, eps);
        wr_log("\t----");
    }

    // Init critical points list and skeleton points list:
    crit_point_list_init(&crit_point_list);
    fcoords_list_init(&skel_point_list);
    highDiv_point_list_init(&highDiv_point_list);



    //
    // PHASE 1 - Perform resampling (if needed):
    //

    // Compute dimensions of resampled and padded volume:
    a_dimx = (int) ((double) (dimx) * scale + 0.5) + a_rad * 2;
    a_dimy = (int) ((double) (dimy) * scale + 0.5) + a_rad * 2;
    a_dimz = (int) ((double) (dimz) * scale + 0.5) + a_rad * 2;

    // Reset output image:
    memset(skl_im, 0, dimx * dimy * dimz * sizeof (unsigned char));

    // Init gradient vector flow:
    P3D_TRY(gvf_x = (float*) calloc(a_dimx * a_dimy*a_dimz, sizeof (float)));
    P3D_TRY(gvf_y = (float*) calloc(a_dimx * a_dimy*a_dimz, sizeof (float)));
    P3D_TRY(gvf_z = (float*) calloc(a_dimx * a_dimy*a_dimz, sizeof (float)));

    // Zero padding for temporary image (values 255 are replaced with 1):	
    P3D_TRY(tmp_im = (unsigned char*) calloc(a_dimx * a_dimy*a_dimz, sizeof (unsigned char)));

    //#pragma omp parallel for private(i, j)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                a = (int) ((i - a_rad) / scale);
                b = (int) ((j - a_rad) / scale);
                c = (int) ((k - a_rad) / scale);

                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = (unsigned char) (in_im [ I(a, b, c, dimx, dimy) ] / 255.0);
            }

    // Print elapsed time for this phase:	
    if (wr_log != NULL) {
        t2 = p3dGetElapsedTime();
        wr_log("\tInitializations and re-scaling performed in %0.3f sec.", t2 - t1);
        t1 = t2;
    }

    //
    // PHASE 2 - Compute gradient vector flow:
    //
    ct = p3dGVF(tmp_im, gvf_x, gvf_y, gvf_z, a_dimx, a_dimy, a_dimz, mu, eps);

    // Print elapsed time for this phase:
    if (wr_log != NULL) {
        t2 = p3dGetElapsedTime();
        wr_log("\tGradient vector flow (%d iterations) computed in %0.3f sec.", ct, t2 - t1);
        t1 = t2;
    }


    //
    // PHASE 3 - Assess critical points:
    //

    // Get and classify critical points:
    P3D_TRY(p3dGetCriticalPoints(tmp_im, gvf_x, gvf_y, gvf_z, a_dimx, a_dimy, a_dimz, &crit_point_list));
    P3D_TRY(p3dClassifyCriticalPoints(crit_point_list, gvf_x, gvf_y, gvf_z, a_dimx, a_dimy, a_dimz));

    // Print elapsed time for this phase:
    if (wr_log != NULL) {
        t2 = p3dGetElapsedTime();
        wr_log("\tCritical points determined and classified in %0.3f sec.", t2 - t1);
        t1 = t2;
    }


    //
    // PHASE 4 - Compute core skeleton following gradient vector flow from saddle points:
    //

    // Sub-voxel step-size and closeness criteria should be reduced in the case of downsampling:
    if (scale < 0.0) step *= scale;
    if (scale < 0.0) close_dist *= scale;

    P3D_TRY(p3dComputeCoreSkeleton(crit_point_list, &skel_point_list, gvf_x, gvf_y, gvf_z,
            a_dimx, a_dimy, a_dimz, step, close_dist));

    // Print elapsed time for this phase:
    if (wr_log != NULL) {
        t2 = p3dGetElapsedTime();
        wr_log("\tCore skeleton computed successfully in %0.3f sec.", t2 - t1);
        t1 = t2;
    }


    //
    // PHASE 5 - Add high divergence points:
    //
    if (fcoords_list_isempty(skel_point_list) == P3D_TRUE) {
        if (wr_log != NULL) {
            wr_log("\tWARNING: The core skeleton is empty!");
        }
    } else if (thresh > 0.0) {
        P3D_TRY(p3dGetHighDivPoints(tmp_im, gvf_x, gvf_y, gvf_z, a_dimx, a_dimy, a_dimz, thresh, &highDiv_point_list));
        P3D_TRY(p3dComputeHierarchicalSkeleton(&highDiv_point_list, &skel_point_list, gvf_x, gvf_y, gvf_z,
                a_dimx, a_dimy, a_dimz, step, close_dist));

        // Print elapsed time for this phase:
        if (wr_log != NULL) {
            t2 = p3dGetElapsedTime();
            wr_log("\tAdditional branches determined successfully in %0.3f sec.", t2 - t1);
            t1 = t2;
        }
    }


    //
    // PHASE 6 - Discretize skeleton points to output matrix (paying attention to previous 
    //           upsampling) and thinning:
    //	
    if (fcoords_list_isempty(skel_point_list) == P3D_FALSE) {
        // Perform voxelization:
        while (fcoords_list_isempty(skel_point_list) == P3D_FALSE) {
            tmp_coords = fcoords_list_pop(&skel_point_list);

            a = (int) ((tmp_coords.x + 0.5 - a_rad) / scale);
            b = (int) ((tmp_coords.y + 0.5 - a_rad) / scale);
            c = (int) ((tmp_coords.z + 0.5 - a_rad) / scale);

            skl_im[ I(a, b, c, dimx, dimy) ] = OBJECT;
        }

        // Perform a thinning procedure in order to get a one-voxel thick output skeleton:		
        p3dThinning(skl_im, dimx, dimy, dimz);

        // Print elapsed time for this phase:
        if (wr_log != NULL) {
            t2 = p3dGetElapsedTime();
            wr_log("\tRe-scaling and voxelization performed in %0.3f sec.", t2 - t1);
            t1 = t2;   

            // Print total time:
            wr_log("Pore3D - GVF skeleton computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
            //wr_log("Pore3D - GVF skeleton computed successfully in %0.3f sec.", p3dGetElapsedTime());
        }
    }


    //
    // Release resources:	
    //
    if (tmp_im != NULL) free(tmp_im);
    if (gvf_x != NULL) free(gvf_x);
    if (gvf_y != NULL) free(gvf_y);
    if (gvf_z != NULL) free(gvf_z);

    // Ensure lists are deleted:
    while (crit_point_list_isempty(crit_point_list) == P3D_FALSE) {
        crit_point_list_pop(&crit_point_list);
    }
    while (fcoords_list_isempty(skel_point_list) == P3D_FALSE) {
        fcoords_list_pop(&skel_point_list);
    }
    while (highDiv_point_list_isempty(highDiv_point_list) == P3D_FALSE) {
        highDiv_point_list_pop(&highDiv_point_list);
    }

    // Return OK:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:	
    if (tmp_im != NULL) free(tmp_im);
    if (gvf_x != NULL) free(gvf_x);
    if (gvf_y != NULL) free(gvf_y);
    if (gvf_z != NULL) free(gvf_z);

    // Ensure lists are deleted:
    while (crit_point_list_isempty(crit_point_list) == P3D_FALSE) {
        crit_point_list_pop(&crit_point_list);
    }
    while (fcoords_list_isempty(skel_point_list) == P3D_FALSE) {
        fcoords_list_pop(&skel_point_list);
    }
    while (highDiv_point_list_isempty(highDiv_point_list) == P3D_FALSE) {
        highDiv_point_list_pop(&highDiv_point_list);
    }

    return P3D_MEM_ERROR;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}