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
#include "Common/p3dRingRemoverCommon.h"

#define P3DISIJBERSPOSTNOV_RINGREMOVER_ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
#define P3DISIJBERSPOSTNOV_RINGREMOVER_MEDIAN(a,n) _p3dSijbersPostnovRingRemover_kth_smallest(a,n,((n)/2))

double _p3dSijbersPostnovRingRemover_kth_smallest(double* a, int n, int k) {
    register int i, j, l, m;
    register double x;

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
                P3DISIJBERSPOSTNOV_RINGREMOVER_ELEM_SWAP(a[i], a[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (j < k) l = i;
        if (k < i) m = j;
    }
    return a[k];
}

// Function used to perform quick sort of the array for the computation
// of the median (see q_sort function of standard library):

/*int _cmp(const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}*/

// Should be thread-safe:

int _getSums_8(unsigned char* p_in_im,
        unsigned char* mask_im,
        unsigned char* p_mask_im,
        int p_dim,
        int i,
        int j,
        int winsize,
        double* sum,
        double* sqrsum
        ) {
    int k;

    // Init sums:
    *sum = 0.0;
    *sqrsum = 0.0;

    // Cycle for each element in line to compute mean and
    // variance of the line. This computation is meaningless
    // if line is not completely included into ROI.
    for (k = i; k < (i + winsize); k++) {

        // Check if element is outside ROI:
        if (mask_im != NULL) {
            if (p_mask_im[ I2(k, j, p_dim) ] == 0)
                // If element is outside ROI break computation to improve
                // performance:
                return 0; // false
        }

        // Compute sums for further use in mean and variance
        // computation:
        *sum += p_in_im[ I2(k, j, p_dim) ];
        *sqrsum += p_in_im[ I2(k, j, p_dim) ] * p_in_im[ I2(k, j, p_dim) ];
    }

    return 1; // true
}



// Should be thread-safe:

int _getSums_16(unsigned short* p_in_im,
        unsigned char* mask_im,
        unsigned char* p_mask_im,
        int p_dim,
        int i,
        int j,
        int winsize,
        double* sum,
        double* sqrsum
        ) {
    int k;

    // Init sums:
    *sum = 0.0;
    *sqrsum = 0.0;

    // Cycle for each element in line to compute mean and
    // variance of the line. This computation is meaningless
    // if line is not completely included into ROI.
    for (k = i; k < (i + winsize); k++) {

        // Check if element is outside ROI:
        if (mask_im != NULL) {
            if (p_mask_im[ I2(k, j, p_dim) ] == 0)
                // If element is outside ROI break computation to improve
                // performance:
                return 0; // false
        }

        // Compute sums for further use in mean and variance
        // computation:
        *sum += p_in_im[ I2(k, j, p_dim) ];
        *sqrsum += p_in_im[ I2(k, j, p_dim) ] * p_in_im[ I2(k, j, p_dim) ];
    }

    return 1; // true
}



// This procedure removes ring artifacts from CT images.

int p3dSijbersPostnovRingRemover2D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const double centerX,
        const double centerY,
        const int winsize,
        const double in_thresh,
        const int iterations,
        const double precision,
        unsigned char* mask_im,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int i, j, k, it_ct; // generic counters
    double sum, sqrsum;

    double mean, variance;
    double tmp_val;
    double thresh = in_thresh*256.0;

    // Polar images:
    unsigned char* p_in_im = NULL; // no need for malloc
    unsigned char* p_mask_im = NULL; // no need for malloc
    unsigned char* c_out_im = NULL; // no need for malloc

    int p_dim;

    // Matrix:
    double* matrix;
    int* matrix_mask;
    double* column;
    int ct;
    int row_ct;
    int prev_rowct;

    // Artifacts vectors:
    double* loc_art;
    double* glob_art;
    int* glob_mask;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dSijbersPostnovRingRemover2D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying Sijbers and Postnov ring remover...");
        wr_log("\tCenter of rings: [%0.1f,%0.1f].", centerX, centerY);
        wr_log("\tWinsize: %d.", winsize);
        wr_log("\tThreshold: %0.3f.", in_thresh);
        wr_log("\tIterations: %d.", iterations);
        if (mask_im != NULL)
            wr_log("\tMask adopted.");
        wr_log("\tPolar/Cartesian precision: %0.3f.", precision);
    }


    // STEP2: Transform in polar coordinates. Remember that memory allocation
    // for output images is performed within the procedure.

    p3dCartesian2polar_8(in_im, &p_in_im, dimx, dimy, centerX, centerY, precision, &p_dim);
    if (mask_im != NULL) {
        p3dCartesian2polar_8(mask_im, &p_mask_im, dimx, dimy, centerX, centerY, precision, &p_dim);
    }



    // STEP3: Artifact template selection.

    // Allocate dynamic memory:
    loc_art = (double*) malloc(winsize * sizeof (double));


    // Matrix should be a dynamic structure (i.e. a list of array) but for
    // performance reasons we adopt a matrix allocated with the maximum
    // dimensions (worst case):
    matrix = (double*) calloc(winsize*p_dim, sizeof (double));
    /**/matrix_mask = (int*) malloc(p_dim * sizeof (int));

    // Now that we know the dimensions of polar image we can allocate
    // the global artifacts vector:
    glob_art = (double*) malloc(p_dim * sizeof (double));
    glob_mask = (int*) malloc(p_dim * sizeof (int));

    for (it_ct = 0; it_ct < iterations; it_ct++) {
        // Initializations at every iterations:
        prev_rowct = 0;
        memset(loc_art, 0, winsize * sizeof (double));
        memset(glob_art, 0, p_dim * sizeof (double));
        memset(glob_mask, 0, p_dim * sizeof (int)); // init to 0 (false)

        // Within a sliding window:
        for (i = 0; i < (p_dim - winsize); i++) {
            // Init counters:
            ct = 0;
            row_ct = 0;

            // Initialization of matrix mask:
            /**/memset(matrix_mask, 0, p_dim * sizeof (int)); // init to 0 (false)

            // For each line of polar image:
#pragma omp parallel for private(sum, sqrsum, mean, variance, k) reduction (+ : row_ct)
            for (j = 0; j < p_dim; j++) {
                // If computation of sum and squared sum is meaningful:
                if (_getSums_8(p_in_im, mask_im, p_mask_im, p_dim, i, j, winsize, &sum, &sqrsum)) {
                    // Compute mean and variance of the line:
                    mean = sum / winsize;
                    variance = sqrsum / winsize - mean*mean;

                    // If variance is below threshold:
                    if ((variance < thresh) && (variance > 0)) {
                        // Set that current position in matrix is meaningful:
                        /**/matrix_mask[j] = 1; // true

                        // For each element subtract mean and add the line
                        // to the matrix:
                        for (k = i; k < (i + winsize); k++) {
                            // Subtract mean from the line and add
                            // it to a matrix
                            /**///matrix[ct++] = p_in_im[ I(k,j,p_dim) ] - mean;
                            matrix[ I2(k - i, j, winsize) ] = p_in_im[ I2(k, j, p_dim) ] - mean;
                        }

                        // Increment the number of lines that meets the
                        // homogeneity criterium:
                        row_ct++;
                    }
                }
            }


            // Compute median for each column and store the value in the
            // artifact vector for this sliding window:

            // Now that we know the number of rows that meets homogeneity
            // criterium we can allocate memory for a column:
            if (row_ct > 0) {
                // Allocate memory for the column (no OpenMP):
                //column = (double*) malloc( row_ct*sizeof(double) );

#pragma omp parallel for private(k, ct, column)
                for (j = 0; j < winsize; j++) {
                    // Allocate memory for the column:
                    column = (double*) malloc(row_ct * sizeof (double));

                    // Fill the column array:
                    /*for ( k = 0; k < row_ct; k++ )
                    {
                            column[k] = matrix[ I(j,k,winsize) ];
                    }*/
                    ct = 0;
                    for (k = 0; k < p_dim; k++) {
                        if (matrix_mask[k] == 1)
                            column[ct++] = matrix[ I2(j, k, winsize) ];
                    }

                    // Order the array:
                    /*qsort( column, row_ct, sizeof(double), _cmp);

                    // Get the median and put in "local" artifact vector:
                    loc_art[j] = column[ row_ct / 2 ];*/
                    loc_art[j] = P3DISIJBERSPOSTNOV_RINGREMOVER_MEDIAN(column, row_ct);

                    // Free memory for the column:
                    free(column);
                }

                // Free memory for the column (no OpenMP):
                //free(column);
            }
            /*else
            {
                    // Set to zero "local" artifact vector:
                    memset(loc_art,0,winsize*sizeof(double));
            }*/


            // Unwrap the "local" artifact vector of dimension W to the
            // "global" artifact vector of dimension P_DIM using the rule
            // based on number of rows that meet the homogeneity criterium:
            for (k = 0; k < winsize; k++) {
                if ((row_ct > prev_rowct) || (glob_mask[k] == 0)) {
                    glob_art[k + i] = loc_art[k];
                    glob_mask[k + i] = 1; // true
                }
            }

            // Set the previous equal to the current:
            prev_rowct = row_ct;
        }

        // The "global" artifact template vector is subtracted from each
        // row of the polar image P:
        if (mask_im != NULL) {
#pragma omp parallel for private(i, tmp_val)
            for (j = 0; j < p_dim; j++) {
                for (i = 0; i < p_dim; i++) {
                    if (p_mask_im[ I2(i, j, p_dim) ] != 0) {
                        // Take care of intensity bounds:
                        tmp_val = p_in_im[ I2(i, j, p_dim) ] - glob_art[i];

                        if (tmp_val < 0.0) tmp_val = 0.0;
                        if (tmp_val > UCHAR_MAX) tmp_val = UCHAR_MAX;

                        p_in_im[ I2(i, j, p_dim) ] = (unsigned char) tmp_val;
                    }
                }
            }
        } else {
#pragma omp parallel for private(i, tmp_val)
            for (j = 0; j < p_dim; j++) {
                for (i = 0; i < p_dim; i++) {
                    // Take care of intensity bounds:
                    tmp_val = p_in_im[ I2(i, j, p_dim) ] - glob_art[i];

                    if (tmp_val < 0.0) tmp_val = 0.0;
                    if (tmp_val > UCHAR_MAX) tmp_val = UCHAR_MAX;

                    p_in_im[ I2(i, j, p_dim) ] = (unsigned char) tmp_val;
                }
            }
        }
    }



    // Return in cartesian coordinates:
    p3dPolar2cartesian_8(p_in_im, &c_out_im, p_dim, centerX, centerY, dimx, dimy);

    // Copy C_OUT_IM to the output of the procedure:
    memcpy(out_im, c_out_im, dimx * dimy * sizeof (unsigned char));

    // Free memory:
    free(glob_art);
    free(glob_mask);
    free(loc_art);
    free(matrix);
    free(matrix_mask);

    free(p_in_im);
    free(c_out_im);

    if (mask_im != NULL) {
        free(p_mask_im);
    }
    
     // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Sijbers and Postnov ring remover applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dSijbersPostnovRingRemover2D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const double centerX,
        const double centerY,
        const int winsize,
        const double in_thresh,
        const int iterations,
        const double precision,
        const int bit12,
        unsigned char* mask_im,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int i, j, k, it_ct; // generic counters
    double sum, sqrsum;

    double mean, variance;
    double tmp_val;
    double thresh;

    // Polar images:
    unsigned short* p_in_im = NULL; // no need for malloc
    unsigned char* p_mask_im = NULL; // no need for malloc
    unsigned short* c_out_im = NULL; // no need for malloc

    int p_dim;

    // Matrix:
    double* matrix;
    int* matrix_mask;
    double* column;
    int ct;
    int row_ct;
    int prev_rowct = 0;

    // Artifacts vectors:
    double* loc_art;
    double* glob_art;
    int* glob_mask;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dSijbersPostnovRingRemover2D_16");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying Sijbers and Postnov ring remover...");
        wr_log("\tCenter of rings: [%0.1f,%0.1f].", centerX, centerY);
        wr_log("\tWinsize: %d.", winsize);
        wr_log("\tThreshold: %0.3f.", in_thresh);
        wr_log("\tIterations: %d.", iterations);
        if (bit12 == P3D_TRUE)
            wr_log("\t12-bit images flag: true.");
        else
            wr_log("\t12-bit images flag: false.");
        if (mask_im != NULL)
            wr_log("\tMask adopted.");        
        wr_log("\tPolar/Cartesian precision: %0.3f.", precision);
    }

    // Set the correct threshold value:
    if (bit12 == P3D_TRUE)
        thresh = in_thresh * 4096.0;
    else
        thresh = in_thresh * 65536.0;


    // STEP2: Transform in polar coordinates. Remember that memory allocation
    // for output images is performed within the procedure.

    p3dCartesian2polar_16(in_im, &p_in_im, dimx, dimy, centerX, centerY, precision, &p_dim);
    if (mask_im != NULL) {
        p3dCartesian2polar_8(mask_im, &p_mask_im, dimx, dimy, centerX, centerY, precision, &p_dim);
    }


    // STEP3: Artifact template selection.

    // Allocate dynamic memory:
    loc_art = (double*) calloc(winsize, sizeof (double));


    // Matrix should be a dynamic structure (i.e. a list of array) but for
    // performance reasons we adopt a matrix allocated with the maximum
    // dimensions (worst case):
    matrix = (double*) calloc(winsize*p_dim, sizeof (double));
    /**/matrix_mask = (int*) malloc(p_dim * sizeof (int));

    // Now that we know the dimensions of polar image we can allocate
    // the global artifacts vector:
    glob_art = (double*) malloc(p_dim * sizeof (double));
    glob_mask = (int*) malloc(p_dim * sizeof (int));

    for (it_ct = 0; it_ct < iterations; it_ct++) {
        // Initializations at every iterations:
        prev_rowct = 0;
        memset(loc_art, 0, winsize * sizeof (double));
        memset(glob_art, 0, p_dim * sizeof (double));
        memset(glob_mask, 0, p_dim * sizeof (int)); // init to 0 (false)

        // Within a sliding window:
        for (i = 0; i < (p_dim - winsize); i++) {
            // Init counters:
            ct = 0;
            row_ct = 0;

            // Initialization of matrix mask:
            /**/memset(matrix_mask, 0, p_dim * sizeof (int)); // init to 0 (false)

            // For each line of polar image:
#pragma omp parallel for private(sum, sqrsum, mean, variance, k) reduction (+ : row_ct)
            for (j = 0; j < p_dim; j++) {
                // If computation of sum and squared sum is meaningful:
                if (_getSums_16(p_in_im, mask_im, p_mask_im, p_dim, i, j, winsize, &sum, &sqrsum)) {
                    // Compute mean and variance of the line:
                    mean = sum / winsize;
                    variance = sqrsum / winsize - mean*mean;

                    // If variance is below threshold:
                    if ((variance < thresh) && (variance > 0)) {
                        // Set that current position in matrix is meaningful:
                        /**/matrix_mask[j] = 1; // true

                        // For each element subtract mean and add the line
                        // to the matrix:
                        for (k = i; k < (i + winsize); k++) {
                            // Subtract mean from the line and add
                            // it to a matrix
                            /**///matrix[ct++] = p_in_im[ I(k,j,p_dim) ] - mean;
                            matrix[ I2(k - i, j, winsize) ] = p_in_im[ I2(k, j, p_dim) ] - mean;
                        }

                        // Increment the number of lines that meets the
                        // homogeneity criterium:
                        row_ct++;
                    }
                }
            }


            // Compute median for each column and store the value in the
            // artifact vector for this sliding window:

            // Now that we know the number of rows that meets homogeneity
            // criterium we can allocate memory for a column:
            if (row_ct > 0) {
                // Allocate memory for the column (no OpenMP):
                //column = (double*) malloc( row_ct*sizeof(double) );

#pragma omp parallel for private(k, ct, column)
                for (j = 0; j < winsize; j++) {
                    // Allocate memory for the column:
                    column = (double*) malloc(row_ct * sizeof (double));

                    // Fill the column array:
                    /*for ( k = 0; k < row_ct; k++ )
                    {
                            column[k] = matrix[ I(j,k,winsize) ];
                    }*/
                    ct = 0;
                    for (k = 0; k < p_dim; k++) {
                        if (matrix_mask[k] == 1)
                            column[ct++] = matrix[ I2(j, k, winsize) ];
                    }

                    // Order the array:
                    /*qsort( column, row_ct, sizeof(double), _cmp);

                    // Get the median and put in "local" artifact vector:
                    loc_art[j] = column[ row_ct / 2 ];*/
                    loc_art[j] = P3DISIJBERSPOSTNOV_RINGREMOVER_MEDIAN(column, row_ct);

                    // Free memory for the column:
                    free(column);
                }

                // Free memory for the column (no OpenMP):
                //free(column);
            }
            /*else
            {
                    // Set to zero "local" artifact vector:
                    memset(loc_art,0,winsize*sizeof(double));
            }*/


            // Unwrap the "local" artifact vector of dimension W to the
            // "global" artifact vector of dimension P_DIM using the rule
            // based on number of rows that meets the homogeneity criterium:
            for (k = 0; k < winsize; k++) {
                if ((row_ct > prev_rowct) || (glob_mask[k] == 0)) {
                    glob_art[k + i] = loc_art[k];
                    glob_mask[k + i] = 1; // true
                }
            }

            // Set the previous equal to the current:
            prev_rowct = row_ct;
        }

        // The "global" artifact template vector is subtracted from each
        // row of the polar image P:
        if (mask_im != NULL) {
#pragma omp parallel for private(i, tmp_val)
            for (j = 0; j < p_dim; j++) {
                for (i = 0; i < p_dim; i++) {
                    if (p_mask_im[ I2(i, j, p_dim) ] != 0) {
                        // Take care of intensity bounds:
                        tmp_val = p_in_im[ I2(i, j, p_dim) ] - glob_art[i];

                        if (tmp_val < 0.0) tmp_val = 0.0;
                        if (tmp_val > USHRT_MAX) tmp_val = USHRT_MAX;

                        p_in_im[ I2(i, j, p_dim) ] = (unsigned short) tmp_val;
                    }
                }
            }
        } else {
#pragma omp parallel for private(i, tmp_val)
            for (j = 0; j < p_dim; j++) {
                for (i = 0; i < p_dim; i++) {
                    // Take care of intensity bounds:
                    tmp_val = p_in_im[ I2(i, j, p_dim) ] - glob_art[i];

                    if (tmp_val < 0.0) tmp_val = 0.0;
                    if (tmp_val > USHRT_MAX) tmp_val = USHRT_MAX;

                    p_in_im[ I2(i, j, p_dim) ] = (unsigned short) tmp_val;
                }
            }
        }
    }


    // Return in cartesian coordinates:
    p3dPolar2cartesian_16(p_in_im, &c_out_im, p_dim, centerX, centerY, dimx, dimy);

    // Copy C_OUT_IM to the output of the procedure:
    memcpy(out_im, c_out_im, dimx * dimy * sizeof (unsigned short));

    // Free memory:
    free(glob_art);
    free(glob_mask);
    free(loc_art);
    free(matrix);
    free(matrix_mask);

    free(p_in_im);
    free(c_out_im);

    if (mask_im != NULL) {
        free(p_mask_im);
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Sijbers and Postnov ring remover applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }


    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}