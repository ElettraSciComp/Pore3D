#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

int p3dRidlerThresholding_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned char* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double tt, tb, to, t2;
    int t;
    long N, no, nb;

    int ct;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dRidlerThresholding_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Ridler's method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    N = dimx * dimy*dimz;
    tb = 0.0;
    to = 0.0;
    no = 0;
    for (ct = 0; ct < dimx * dimy * dimz; ct++)
        to = to + (in_im[ct]);
    tt = (to / (double) N);

    while (N) {
        no = 0;
        nb = 0;
        tb = 0.0;
        to = 0.0;
        for (ct = 0; ct < dimx * dimy * dimz; ct++) {
            if ((double) (in_im[ct]) >= tt) {
                to = to + (double) (in_im[ct]);
                no++;
            }
            else {
                tb = tb + (double) (in_im[ct]);
                nb++;
            }
        }

        if (no == 0) no = 1;
        if (nb == 0) nb = 1;
        t2 = (tb / (double) nb + to / (double) no) / 2.0;
        if (t2 == tt) N = 0;
        tt = t2;
    }
    t = (int) tt;

    *thresh = (unsigned char) t;



#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;

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

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

int p3dRidlerThresholding_16(
        unsigned short* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned short* thresh,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {

    double tt, tb, to, t2;
    int t;
    long N, no, nb;

    int ct;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dRidlerThresholding_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Ridler's iterative method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    N = dimx * dimy*dimz;
    tb = 0.0;
    to = 0.0;
    no = 0;
    for (ct = 0; ct < dimx * dimy * dimz; ct++)
        to = to + (in_im[ct]);
    tt = (to / (double) N);

    while (N) {
        no = 0;
        nb = 0;
        tb = 0.0;
        to = 0.0;
        for (ct = 0; ct < dimx * dimy * dimz; ct++) {
            if ((double) (in_im[ct]) >= tt) {
                to = to + (double) (in_im[ct]);
                no++;
            }
            else {
                tb = tb + (double) (in_im[ct]);
                nb++;
            }
        }

        if (no == 0) no = 1;
        if (nb == 0) nb = 1;
        t2 = (tb / (double) nb + to / (double) no) / 2.0;
        if (t2 == tt) N = 0;
        tt = t2;
    }
    t = (int) tt;

    *thresh = (unsigned short) t;



#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;

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

    // Return error:
    return (int) P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}
