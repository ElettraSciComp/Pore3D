#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

double _p3dKittlerTresholding_log(double x) {
    if (x > 0.0)
        return log(x);
    else
        return 0.0;
}

double _p3dKittlerTresholding_P1(int t, double* h) {
    long i;
    double sum;

    sum = 0.0;
    for (i = 0; i <= t; i++)
        sum += h[i];

    return sum;
}

double _p3dKittlerTresholding_P2(int t, double* h, int hist_length) {
    long i;
    double sum;

    sum = 0.0;
    for (i = t + 1; i <= hist_length; i++)
        sum += h[i];

    return sum;
}

double _p3dKittlerTresholding_u1(int t, double* h) {
    int i;
    double sum, p;

    sum = 0.0;
    p = _p3dKittlerTresholding_P1(t, h);
    if (p <= 0.0)
        return 0.0;

    for (i = 0; i <= t; i++)
        sum += h[i] * i;

    return sum / p;
}

double _p3dKittlerTresholding_u2(int t, double* h, int hist_length) {
    int i;
    double sum, p;

    sum = 0.0;
    p = _p3dKittlerTresholding_P2(t, h, hist_length);
    if (p <= 0.0)
        return 0.0;

    for (i = t + 1; i <= hist_length; i++)
        sum += h[i] * i;

    return sum / p;
}

double _p3dKittlerTresholding_s1(int t, double* h) {
    int i;
    double sum, p, u, x;

    sum = 0.0;
    p = _p3dKittlerTresholding_P1(t, h);

    if (p <= 0.0)
        return 0.0;

    u = _p3dKittlerTresholding_u1(t, h);
    for (i = 0; i <= t; i++) {
        x = (i - u)*(i - u);
        sum += x * h[i];
    }

    return sum / p;
}

double _p3dKittlerTresholding_s2(int t, double* h, int hist_length) {
    int i;
    double sum, p, u, x;

    sum = 0.0;
    p = _p3dKittlerTresholding_P2(t, h, hist_length);

    if (p <= 0.0)
        return 0.0;

    u = _p3dKittlerTresholding_u2(t, h, hist_length);
    for (i = t + 1; i <= hist_length; i++) {
        x = (i - u)*(i - u);
        sum += x * h[i];
    }

    return sum / p;
}

double _p3dKittlerTresholding_J(int t, double* h, int hist_length) {
    double a, b, c, d, x1;

    a = _p3dKittlerTresholding_P1(t, h);
    b = _p3dKittlerTresholding_s1(t, h);
    c = _p3dKittlerTresholding_P2(t, h, hist_length);
    d = _p3dKittlerTresholding_s2(t, h, hist_length);

    x1 = 1.0 + 2.0 * (a * _p3dKittlerTresholding_log(b) + c * _p3dKittlerTresholding_log(d)) -
            2.0 * (a * _p3dKittlerTresholding_log(a) + c * _p3dKittlerTresholding_log(c));

    return x1;
}

int p3dKittlerThresholding_8(
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
    double* F;
    int i, ct;    
    int tbest = 0;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dKittlerThresholding_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Kittler and Illingworth method...");
    }

    /* Allocate and initialize to zero kernel histogram: */
    P3D_TRY(prob = (double*) calloc((UCHAR_MAX + 1), sizeof (double)));
    P3D_TRY(F = (double*) calloc((UCHAR_MAX + 1), sizeof (double)));
    
    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;


    /* Compute thresholding according to Kittler's minimum error thresholding: */
    for (i = 1; i <= UCHAR_MAX; i++) {
        F[i] = _p3dKittlerTresholding_J(i, prob, UCHAR_MAX);
        if (F[i] < F[tbest])
            tbest = i;
    }

    *thresh = (unsigned char) tbest;        


    /* Threshold image: */
    //#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;
    
    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    } 


    // Free memory:
    if (prob != NULL) free(prob);
    if (F != NULL) free(F);

        // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

 // Free memory:
    if (prob != NULL) free(prob);
    if (F != NULL) free(F);

    // Return error:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

int p3dKittlerThresholding_16(
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
    double* F;
    int i, ct;
    int tbest = 0;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dKittlerThresholding_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Thresholding image according to Kittler and Illingworth method...");
    }


    /* Allocate and initialize to zero kernel histogram: */
    P3D_TRY(prob = (double*) calloc((USHRT_MAX + 1), sizeof (double)));
    P3D_TRY(F = (double*) calloc((USHRT_MAX + 1), sizeof (double)));

    /* Compute image histogram: */
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        prob[in_im[ ct ]] = prob[in_im[ ct ]] + 1.0;


    /* Compute thresholding according to Kittler's minimum error thresholding: */
    for (i = 1; i <= USHRT_MAX; i++) {
        F[i] = _p3dKittlerTresholding_J(i, prob, USHRT_MAX);
        if (F[i] < F[tbest])
            tbest = i;
    }

    *thresh = (unsigned short) tbest;

    /* Threshold image: */
    //#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        out_im[ ct ] = (in_im[ ct ] > (*thresh)) ? OBJECT : BACKGROUND;

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tDetermined threshold: %d.", *thresh);
        wr_log("Pore3D - Image thresholded successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (prob != NULL) free(prob);
    if (F != NULL) free(F);

        // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

 // Free memory:
    if (prob != NULL) free(prob);
    if (F != NULL) free(F);

    // Return error:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }
	
    return P3D_AUTH_ERROR;*/

}
