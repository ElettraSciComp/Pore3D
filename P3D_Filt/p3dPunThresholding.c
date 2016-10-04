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
#include <math.h>

#include "p3dFilt.h"
#include "p3dTime.h"

double _p3dPunThresholding_entropy(double *h, int a) {
    if (h[a] > 0.0)
        return -(h[a] * log((double) h[a]));
    return 0.0;
}

double _p3dPunThresholding_flog(double x) {
    if (x <= 0.0)
        return 0.0;
    return log((double) x);
}

double _p3dPunThresholding_maxtot(double *h, int i) {
    double x;
    int j;

    x = h[0];
    for (j = 1; j <= i; j++)
        if (x < h[j])
            x = h[j];
    return x;
}

double _p3dPunThresholding_maxfromt_8(double *h, int i) {
    int j;
    double x;

    x = h[i + 1];
    for (j = i + 2; j <= UCHAR_MAX; j++)
        if (x < h[j])
            x = h[j];
    return x;
}

double _p3dPunThresholding_maxfromt_16(double *h, int i) {
    int j;
    double x;

    x = h[i + 1];
    for (j = i + 2; j <= USHRT_MAX; j++)
        if (x < h[j])
            x = h[j];
    return x;
}


int p3dPunThresholding_8(
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
    double *Ht, *Pt, *F;
    double HT, x, y, z, to, from;
    int i, ct, t;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dPunThresholding_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Pun's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    Ht = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));
    Pt = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));
    F = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));
    prob = (double*) calloc((UCHAR_MAX + 1), sizeof (double));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= UCHAR_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute the factors */
    HT = Ht[0] = _p3dPunThresholding_entropy(prob, 0);
    Pt[0] = prob[0];
    for (i = 1; i <= UCHAR_MAX; i++) {
        Pt[i] = Pt[i - 1] + prob[i];
        x = _p3dPunThresholding_entropy(prob, i);
        Ht[i] = Ht[i - 1] + x;
        HT += x;
    }

    /* Calculate the function to be maximized at all levels */
    t = 0;
    for (i = 0; i <= UCHAR_MAX; i++) {
        to = (_p3dPunThresholding_maxtot(prob, i));
        from = _p3dPunThresholding_maxfromt_8(prob, i);
        if (to > 0.0 && from > 0.0) {
            x = (Ht[i] / HT) * _p3dPunThresholding_flog(Pt[i]) / _p3dPunThresholding_flog(to);
            y = 1.0 - (Ht[i] / HT);
            z = _p3dPunThresholding_flog(1 - Pt[i]) / _p3dPunThresholding_flog(from);
        } else x = y = z = 0.0;
        F[i] = x + y*z;
        if (i > 0 && F[i] > F[t]) t = i;
    }

    *thresh = (unsigned char) t;



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
    if (Ht != NULL) free(Ht);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);

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

int p3dPunThresholding_16(
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
    double *Ht, *Pt, *F;
    double HT, x, y, z, to, from;
    int i, ct, t;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dPunThresholding_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Pun's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    Ht = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));
    Pt = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));
    F = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));
    prob = (double*) calloc((USHRT_MAX + 1), sizeof (double));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= USHRT_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute the factors */
    HT = Ht[0] = _p3dPunThresholding_entropy(prob, 0);
    Pt[0] = prob[0];
    for (i = 1; i <= USHRT_MAX; i++) {
        Pt[i] = Pt[i - 1] + prob[i];
        x = _p3dPunThresholding_entropy(prob, i);
        Ht[i] = Ht[i - 1] + x;
        HT += x;
    }

    /* Calculate the function to be maximized at all levels */
    t = 0;
    for (i = 0; i <= USHRT_MAX; i++) {
        to = (_p3dPunThresholding_maxtot(prob, i));
        from = _p3dPunThresholding_maxfromt_16(prob, i);
        if (to > 0.0 && from > 0.0) {
            x = (Ht[i] / HT) * _p3dPunThresholding_flog(Pt[i]) / _p3dPunThresholding_flog(to);
            y = 1.0 - (Ht[i] / HT);
            z = _p3dPunThresholding_flog(1 - Pt[i]) / _p3dPunThresholding_flog(from);
        } else x = y = z = 0.0;
        F[i] = x + y*z;
        if (i > 0 && F[i] > F[t]) t = i;
    }

    *thresh = (unsigned short) t;



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
    if (Ht != NULL) free(Ht);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);

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
