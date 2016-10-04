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
#include <limits.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"


// NOTE: Different implementation for the 8 bit case and the
// 16 bit case. 

#define P3DMEDIAN_FILTER_ELEM_SWAP(a,b) { register unsigned short t=(a);(a)=(b);(b)=t; }
#define P3DMEDIAN_FILTER_MEDIAN(a,n) _p3dMedianFilter_kth_smallest(a,n,((n)/2))

unsigned short _p3dMedianFilter_kth_smallest(unsigned short* a, int n, int k) {
    register int i, j, l, m;
    register unsigned short x;

    l = 0;
    m = n - 1;
    while (l < m) {
        x = a[k];
        i = l;
        j = m;
        do {
            while (a[i] < x) i++;
            while (x < a[j]) j--;
            if (i <= j) {
                P3DMEDIAN_FILTER_ELEM_SWAP(a[i], a[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (j < k) l = i;
        if (k < i) m = j;
    }
    return a[k];
}

int p3dMedianFilter2D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned char* tmp_im;

    int a_dimx, a_dimy;
    int i, j;
    int x, y;
    int pr = 0;

    // Variables for computing gaussian kernel:
    int a_rad, ct;

    // Temporary array:
    unsigned int* hist;
    unsigned int sum;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying median filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;

    // Get the replicate padded input:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * sizeof (unsigned char)));
    P3D_TRY(p3dReplicatePadding2D_8(in_im, tmp_im, dimx, dimy, a_rad, NULL, NULL));




    // Volume scanning:
#pragma omp parallel for private(i, x, y, ct, sum, hist) reduction( + : pr) 
    for (j = a_rad; j < (a_dimy - a_rad); j++) {
        // Allocate and initialize to zero kernel histogram:
        hist = (unsigned int*) calloc((UCHAR_MAX + 1), sizeof (unsigned int));

        // Compute histogram for first step:
        for (y = (j - a_rad); y <= (j + a_rad); y++)
            for (x = 0; x <= (2 * a_rad); x++) {
                hist[tmp_im[ I2(x, y, a_dimx) ]] = hist[tmp_im[ I2(x, y, a_dimx) ]] + 1;
            }


        // Compute median:
        ct = -1;
        sum = 0;
        while (sum <= ((size * size) / 2)) {
            sum += hist[++ct];
        }

        // Set out voxel with the median:
        out_im[ I2(0, j - a_rad, dimx) ] = ct;

        // Increment progress counter:
        pr++;


        // Scan along x dimension:
        for (i = (a_rad + 1); i < (a_dimx - a_rad); i++) {
            // Update "sliding" histogram:
            for (y = (j - a_rad); y <= (j + a_rad); y++) {
                hist[tmp_im[ I2(i - a_rad - 1, y, a_dimx) ]]--;
                hist[tmp_im[ I2(i + a_rad, y, a_dimx) ]]++;
            }

            // Compute median:
            ct = -1;
            sum = 0;
            while (sum <= ((size * size) / 2)) {
                sum += hist[++ct];
            }

            // Set out voxel with the median:
            out_im[ I2(i - a_rad, j - a_rad, dimx) ] = ct;

            // Increment progress counter:
            pr++;
        }

        // Clear histogram:
        if (hist != NULL) free(hist);

        // Update any progress bar:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy)*100 + 0.5));
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
         wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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
}

int p3dMedianFilter2D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned short* tmp_im;

    int a_dimx, a_dimy;
    int i, j;
    int x, y;

    // Variables for computing gaussian kernel:
    int a_rad, ct;

    // Temporary array:
    unsigned short* a_vtmp;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying median filter...");
        wr_log("\tKernel size: %d.", size);
    }

    // Init variables:
    a_rad = size / 2; // integer division

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;

    // Get the replicate padded input:
    tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * sizeof (unsigned short));
    p3dReplicatePadding2D_8(in_im, tmp_im, dimx, dimy, a_rad, NULL, NULL);

    // Allocate temp array:
    a_vtmp = (unsigned short*) malloc(size * size * sizeof (unsigned short));

    // Volume scanning:
#pragma omp parallel for private(i, x, y, ct)
    for (j = a_rad; j < (a_dimy - a_rad); j++)
        for (i = a_rad; i < (a_dimx - a_rad); i++) {
            ct = 0;
            // Fill temporary array:
            for (y = (j - a_rad); y <= (j + a_rad); y++)
                for (x = (i - a_rad); x <= (i + a_rad); x++)
                    a_vtmp[ct++] = tmp_im[ I2(x, y, a_dimx) ];

            out_im[ I2(i - a_rad, j - a_rad, dimx) ] = P3DMEDIAN_FILTER_MEDIAN(a_vtmp, size * size);
        }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
         wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (tmp_im != NULL) free(tmp_im);
    if (a_vtmp != NULL) free(a_vtmp);

    // Return OK:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);

    // Return error:
    return P3D_MEM_ERROR;
}

int p3dMedianFilter3D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned char* tmp_im;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int pr = 0;

    // Variables for computing kernel:
    int a_rad, ct;

    // Temporary array:
    unsigned int* hist;
    unsigned int sum;
    
    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dMedianFilter3D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying median filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Get the replicate padded input:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
    P3D_TRY(p3dReplicatePadding3D_8(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));




    // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, ct, sum, hist) reduction( + : pr)
    for (k = a_rad; k < (a_dimz - a_rad); k++) {
        for (j = a_rad; j < (a_dimy - a_rad); j++) {
            // Allocate and initialize to zero kernel histogram:
            hist = (unsigned int*) calloc((UCHAR_MAX + 1), sizeof (unsigned int));

            // Compute histogram for first step:
            for (z = (k - a_rad); z <= (k + a_rad); z++)
                for (y = (j - a_rad); y <= (j + a_rad); y++)
                    for (x = 0; x <= (2 * a_rad); x++) {
                        hist[tmp_im[ I(x, y, z, a_dimx, a_dimy) ]]++;
                    }


            // Compute median:
            ct = -1;
            sum = 0;
            while (sum <= ((size * size * size) / 2)) {
                sum += hist[++ct];
            }

            // Set out voxel with the median:
            out_im[ I(0, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

            // Increase progress counter:
            pr++;


            // Scan along x dimension:
            for (i = (a_rad + 1); i < (a_dimx - a_rad); i++) {
                // Update "sliding" histogram:
                for (z = (k - a_rad); z <= (k + a_rad); z++)
                    for (y = (j - a_rad); y <= (j + a_rad); y++) {
                        hist[tmp_im[ I(i - a_rad - 1, y, z, a_dimx, a_dimy) ]]--;
                        hist[tmp_im[ I(i + a_rad, y, z, a_dimx, a_dimy) ]]++;
                    }

                // Compute median:
                ct = -1;
                sum = 0;
                while (sum <= ((size * size * size) / 2)) {
                    sum += hist[++ct];
                }

                // Set out voxel with the median:
                out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

                // Increase progress counter:
                pr++;
            }

            // Clear histogram:
            if (hist != NULL) free(hist);
        }

        // Update any progress bar:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy * dimz)*100 + 0.5));
    }



    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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

int p3dMedianFilter3D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned short* tmp_im;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;

    // Variables for computing gaussian kernel:
    int a_rad, ct;

    // Temporary array:
    unsigned short* a_vtmp;
    
    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dMedianFilter3D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying median filter...");
        wr_log("\tKernel size: %d.", size);
    }

    // Init variables:
    a_rad = size / 2; // integer division

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Get the replicate padded input:
    tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned short));
    p3dReplicatePadding3D_16(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL);

    // Allocate temp array:
    a_vtmp = (unsigned short*) malloc(size * size * size * sizeof (unsigned short));

    // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                ct = 0;
                // Fill temporary array:
                for (z = (k - a_rad); z <= (k + a_rad); z++)
                    for (y = (j - a_rad); y <= (j + a_rad); y++)
                        for (x = (i - a_rad); x <= (i + a_rad); x++)
                            a_vtmp[ct++] = tmp_im[ I(x, y, z, a_dimx, a_dimy) ];

                out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = P3DMEDIAN_FILTER_MEDIAN(a_vtmp, size * size * size);
            }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (tmp_im != NULL) free(tmp_im);
    if (a_vtmp != NULL) free(a_vtmp);

    // Return OK:
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