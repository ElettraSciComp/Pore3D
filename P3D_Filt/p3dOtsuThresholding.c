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

#include "p3dFilt.h"
#include "p3dTime.h"

double _p3dOtsuThresholding_w(double *p, int k) {
    int i;
    double x = 0.0;

    for (i = 0; i < k; i++)
        x += p[i];

    return x;
}

double _p3dOtsuThresholding_u(double *p, int k) {
    int i;
    double x = 0.0;

    for (i = 0; i < k; i++)
        x += (double) i * p[i];
    return x;
}

double _p3dOtsuThresholding_nu(double *p, int k, double ut, double vt) {
    double x, y;

    y = _p3dOtsuThresholding_w(p, k);
    x = ut * y - _p3dOtsuThresholding_u(p, k);
    x = x*x;
    y = y * (1.0 - y);
    if (y > 0)
        x = x / y;
    else
        x = 0.0;

    return x / vt;
}


int p3dOtsuThresholding_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned char* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double* prob;
    int i, j, k, m;
    double y, z;
    double ut, vt;
    int ct;

    /*char auth_code;
        
    //
    // Authenticate:
    //
    auth_code = authenticate("p3dOtsuThresholding_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Otsu's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    P3D_TRY(prob = (double*) calloc((UCHAR_MAX + 1), sizeof (double)));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= UCHAR_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute global mean: */
    ut = _p3dOtsuThresholding_u(prob, UCHAR_MAX);

    /* Copute global variance: */
    vt = 0.0;
    for (i = 0; i <= UCHAR_MAX; i++)
        vt += (i - ut)*(i - ut) * prob[i];

    j = -1;
    k = -1;
    for (i = 0; i <= UCHAR_MAX; i++) {
        if ((j < 0) && (prob[i] > 0.0))
            /* First index handling: */
            j = i;
        if (prob[i] > 0.0)
            /* Last index handling: */
            k = i;
    }
    z = -1.0;
    m = -1;
    for (i = j; i <= k; i++) {
        /* Compute NU: */
        y = _p3dOtsuThresholding_nu(prob, i, ut, vt);
        /* Check if it is the biggest and save value: */
        if (y >= z) {
            z = y;
            m = i;
        }
    }

    *thresh = (unsigned char) m;

    #pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (prob != NULL) free(prob);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Free memory:
    if (prob != NULL) free(prob);

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

int p3dOtsuThresholding_16(
        unsigned short* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned short* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double* prob;
    int i, j, k, m;
    double y, z;
    double ut, vt;
    int ct;

    /*char auth_code;
        
    //
    // Authenticate:
    //
    auth_code = authenticate("p3dOtsuThresholding_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Otsu's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    P3D_TRY(prob = (double*) calloc((USHRT_MAX + 1), sizeof (double)));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= USHRT_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute global mean: */
    ut = _p3dOtsuThresholding_u(prob, USHRT_MAX);

    /* Copute global variance: */
    vt = 0.0;
    for (i = 0; i <= USHRT_MAX; i++)
        vt += (i - ut)*(i - ut) * prob[i];

    j = -1;
    k = -1;
    for (i = 0; i <= USHRT_MAX; i++) {
        if ((j < 0) && (prob[i] > 0.0))
            /* First index handling: */
            j = i;
        if (prob[i] > 0.0)
            /* Last index handling: */
            k = i;
    }
    z = -1.0;
    m = -1;
    for (i = j; i <= k; i++) {
        /* Compute NU: */
        y = _p3dOtsuThresholding_nu(prob, i, ut, vt);
        /* Check if it is the biggest and save value: */
        if (y >= z) {
            z = y;
            m = i;
        }
    }

    *thresh = (unsigned short) m;

    #pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (prob != NULL) free(prob);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Free memory:
    if (prob != NULL) free(prob);

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}
