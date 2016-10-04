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

double _p3dJohannsenThresholding_entropy(double h) {
    if (h > 0.0)
        return (-h * log(h));
    else return 0.0;
}

double _p3dJohannsenThresholding_dlog(double x) {
    if (x > 0.0) return log(x);
    else return 0.0;
}

int p3dJohannsenThresholding_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned char* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double *prob = NULL;
    double *Pt = NULL;
    double *F = NULL;
    double *Pq = NULL;
    int i, t, start, end;
    double Sb, Sw;
    int ct;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dJohannsenThresholding_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Johannsen's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    prob = (double*) calloc((UCHAR_MAX + 1), sizeof (double));
    Pt = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));
    F = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));
    Pq = (double *) malloc(sizeof (double) *(UCHAR_MAX + 1));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= UCHAR_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute the factors */
    Pt[0] = prob[0];
    Pq[0] = 1.0 - Pt[0];
    for (i = 1; i <= UCHAR_MAX; i++) {
        Pt[i] = Pt[i - 1] + prob[i];
        Pq[i] = 1.0 - Pt[i - 1];
    }

    start = 0;
    while (prob[start++] <= 0.0);
    end = UCHAR_MAX;
    while (prob[end--] <= 0.0);

    /* Calculate the function to be minimized at all levels */
    t = -1;
    for (i = start; i <= end; i++) {
        if (prob[i] <= 0.0) continue;
        Sb = _p3dJohannsenThresholding_dlog(Pt[i]) + (1.0 / Pt[i])*
                (_p3dJohannsenThresholding_entropy(prob[i]) + _p3dJohannsenThresholding_entropy(Pt[i - 1]));
        Sw = _p3dJohannsenThresholding_dlog(Pq[i]) + (1.0 / Pq[i])*
                (_p3dJohannsenThresholding_entropy(prob[i]) + _p3dJohannsenThresholding_entropy(Pq[i + 1]));
        F[i] = Sb + Sw;
        if (t < 0) t = i;
        else if (F[i] < F[t]) t = i;
    }


    *thresh = (unsigned char) t;



#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;


    // Free memory:
    if (prob != NULL) free(prob);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);
    if (Pq != NULL) free(Pq);

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Free memory:
    if (prob != NULL) free(prob);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

int p3dJohannsenThresholding_16(
        unsigned short* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned short* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double *prob = NULL;
    double *Pt = NULL;
    double *F = NULL;
    double *Pq = NULL;
    int i, t, start, end;
    double Sb, Sw;
    int ct;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dJohannsenThresholding_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Johannsen's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    prob = (double*) calloc((USHRT_MAX + 1), sizeof (double));
    Pt = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));
    F = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));
    Pq = (double *) malloc(sizeof (double) *(USHRT_MAX + 1));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;

    /* Compute probabilities: */
    for (ct = 0; ct <= USHRT_MAX; ct++)
        prob[ct] = prob[ct] / (double) (dimx * dimy * dimz);

    /* Compute the factors */
    Pt[0] = prob[0];
    Pq[0] = 1.0 - Pt[0];
    for (i = 1; i <= USHRT_MAX; i++) {
        Pt[i] = Pt[i - 1] + prob[i];
        Pq[i] = 1.0 - Pt[i - 1];
    }

    start = 0;
    while (prob[start++] <= 0.0);
    end = USHRT_MAX;
    while (prob[end--] <= 0.0);

    /* Calculate the function to be minimized at all levels */
    t = -1;
    for (i = start; i <= end; i++) {
        if (prob[i] <= 0.0) continue;
        Sb = _p3dJohannsenThresholding_dlog(Pt[i]) + (1.0 / Pt[i])*
                (_p3dJohannsenThresholding_entropy(prob[i]) + _p3dJohannsenThresholding_entropy(Pt[i - 1]));
        Sw = _p3dJohannsenThresholding_dlog(Pq[i]) + (1.0 / Pq[i])*
                (_p3dJohannsenThresholding_entropy(prob[i]) + _p3dJohannsenThresholding_entropy(Pq[i + 1]));
        F[i] = Sb + Sw;
        if (t < 0) t = i;
        else if (F[i] < F[t]) t = i;
    }


    *thresh = (unsigned short) t;



#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;


    // Free memory:
    if (prob != NULL) free(prob);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);
    if (Pq != NULL) free(Pq);

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Free memory:
    if (prob != NULL) free(prob);
    if (Pt != NULL) free(Pt);
    if (F != NULL) free(F);
    if (Pq != NULL) free(Pq);

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}
