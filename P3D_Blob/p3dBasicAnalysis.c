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

/* given a 3-dimensional binary image bin_image of size n[0..2],
   the image is convolved with a filter mask F and
   the grey-tone histogram h[0..255] of the convolved image is returned 
 */
void ghist(
        unsigned char* im,
        double* h, //	h = (double *) calloc ( (UCHAR_MAX+1)*sizeof(double) );
        int dimx,
        int dimy,
        int dimz
        ) {
    int i, j, k;
    int l;
    int ct;

    // Convert image to int:
    unsigned char* tmp_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char));

#pragma omp parallel for private( ct )
    for (ct = 0; ct < (dimx * dimy * dimz); ct++)
        tmp_im [ ct ] = (im [ ct ] == OBJECT) ? 1 : 0;

    // Compute histogram:
    for (i = 0; i < (dimx - 1); i++)
        for (j = 0; j < (dimy - 1); j++) {
            l = tmp_im[ I(i, j, 0, dimx, dimy) ] + (tmp_im[ I(i + 1, j, 0, dimx, dimy) ] << 1)
                    + (tmp_im[ I(i, j + 1, 0, dimx, dimy) ] << 2) + (tmp_im[ I(i + 1, j + 1, 0, dimx, dimy) ] << 3);

            for (k = 0; k < (dimz - 1); k++) {
                l += (tmp_im[ I(i, j, k + 1, dimx, dimy) ] << 4) + (tmp_im[ I(i + 1, j, k + 1, dimx, dimy) ] << 5)
                        + (tmp_im[ I(i, j + 1, k + 1, dimx, dimy) ] << 6) + (tmp_im[ I(i + 1, j + 1, k + 1, dimx, dimy) ] << 7);

                h[l] = h[l] + 1.0;
                l >>= 4;
            }
        }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);

}

/* returns an estimate of the volume fraction V_V from the vector h[0..255] 
   of absolute frequencies of neighborhood configurations of a binary image
 */
double volfrac(double *h) {
    int l;
    double iVol = 0.0;
    double iVol1 = 0.0;

    //#pragma omp parallel reduction (+ : iVol, iVol1)
    for (l = 0; l <= UCHAR_MAX; l++) {
        iVol += h[l];
        if (l == (l | 1))
            iVol1 += h[l];
    }

    return iVol1 / iVol;
}

/* returns the specific surface area S_V from the gray-tone histogram h[0..255], 
    the grid spacing Delta[0..2] is input
 */
double specsurf(double* h, double* Delta) {
    int kl[13][2] = {
        {1, 2},
        {1, 4},
        {1, 16},
        {1, 8},
        {2, 4},
        {1, 32},
        {2, 16},
        {1, 64},
        {4, 16},
        {1, 128},
        {2, 64},
        {4, 32},
        {8, 16}
    };
    double c[13] = {0.045778, 0.045778, 0.045778, 0.036981, 0.036981, 0.036981,
        0.036981, 0.036981, 0.036981, 0.035196, 0.035196, 0.035196,
        0.035196};
    double S_V = 0.0;
    int l, ny;
    double iVol = 0.0;
    double r[13];

    r[0] = Delta[0];
    r[1] = Delta[1];
    r[2] = Delta[2];
    r[3] = r[4] = sqrt(Delta[0] * Delta[0] + Delta[1] * Delta[1]);
    r[5] = r[6] = sqrt(Delta[0] * Delta[0] + Delta[2] * Delta[2]);
    r[7] = r[8] = sqrt(Delta[1] * Delta[1] + Delta[2] * Delta[2]);
    r[9] = r[10] = r[11] = r[12] = sqrt(Delta[0] * Delta[0] + Delta[1] * Delta[1] + Delta[2] * Delta[2]);

    //#pragma omp parallel for private ( ny ) reduction (+ : iVol, S_V)
    for (l = 0; l <= UCHAR_MAX; l++) {
        iVol += h[l];

        for (ny = 0; ny < 13; ny++)
            S_V += h[l] * c[ny] / r[ny]*((l == (l | kl[ny][0]))*(0 == (l & kl[ny][1]))
                +(l == (l | kl[ny][1]))*(0 == (l & kl[ny][0])));
    }

    return 4.0 * S_V / iVol;
}

/* returns the specific integral of mean curvature M_V from the gray-tone 
   histogram h[0..255], the grid spacing Delta[0..2] is input
 */
double specimc(double *h, double *Delta) {
    int kr[9][4] = {
        {1, 2, 4, 8},
        {1, 2, 16, 32},
        {1, 4, 16, 64},
        {1, 2, 64, 128},
        {4, 16, 8, 32},
        {1, 32, 4, 128},
        {2, 8, 16, 64},
        {2, 4, 32, 64},
        {1, 16, 8, 128}
    };
    int kt[8][3] = {
        {1, 64, 32},
        {2, 16, 128},
        {8, 64, 32},
        {4, 16, 128},
        {2, 4, 128},
        {8, 1, 64},
        {2, 4, 16},
        {8, 1, 32}
    };
    double c[13] = {0.045778, 0.045778, 0.045778, 0.036981, 0.036981, 0.036981,
        0.036981, 0.036981, 0.036981, 0.035196, 0.035196, 0.035196,
        0.035196};
    double delta01 = sqrt(Delta[0] * Delta[0] + Delta[1] * Delta[1]);
    double delta02 = sqrt(Delta[0] * Delta[0] + Delta[2] * Delta[2]);
    double delta12 = sqrt(Delta[1] * Delta[1] + Delta[2] * Delta[2]);
    double s = (delta01 + delta02 + delta12) / 2;
    double a[13];
    double M_V = 0.0;
    int ir, l, ny;
    double iVol = 0.0;

    a[0] = Delta[0] * Delta[1];
    a[1] = Delta[0] * Delta[2];
    a[2] = Delta[1] * Delta[2];
    a[3] = a[4] = Delta[2] * delta01;
    a[5] = a[6] = Delta[1] * delta02;
    a[7] = a[8] = Delta[0] * delta12;
    a[9] = a[10] = a[11] = a[12] = 2 * sqrt(s * (s - delta01)*(s - delta02)*(s - delta12));

    //#pragma omp parallel for private ( ny, ir ) reduction (+ : iVol, M_V)
    for (l = 0; l <= UCHAR_MAX; l++) {
        iVol += h[l];

        for (ny = 0; ny < 9; ny++)
            for (ir = 0; ir < 4; ir++)
                M_V += (double) h[l] * c[ny] / (4.0 * a[ny])
                *((l == (l | kr[ny][ir]))*(0 == (l & kr[ny][(ir + 1) % 4]))
                    *(0 == (l & kr[ny][(ir + 2) % 4]))*(0 == (l & kr[ny][(ir + 3) % 4]))
                    -(l == (l | kr[ny][ir]))*(l == (l | kr[ny][(ir + 1) % 4]))
                    *(l == (l | kr[ny][(ir + 2) % 4]))*(0 == (l & kr[ny][(ir + 3) % 4])));

        for (ny = 9; ny < 13; ny++)
            for (ir = 0; ir < 3; ir++)
                M_V += (double) h[l] * c[ny] / (3.0 * a[ny])
                *((l == (l | kt[ny - 9][ir]))*(0 == (l & kt[ny - 9][(ir + 1) % 3]))
                    *(0 == (l & kt[ny - 9][(ir + 2) % 3]))
                    -(l == (l | kt[ny - 5][ir]))*(l == (l | kt[ny - 5][(ir + 1) % 3]))
                    *(0 == (l & kt[ny - 5][(ir + 2) % 3])));
    }

    return 4.0 * M_PI * M_V / iVol;
}

/* returns the specific Euler Characteristic C_V from the vector h[0..255] 
   of absolute freqencies, the lattice distances Delta[0..2] are input
 */
double euler(double *h, double *Delta) {
    int l;
    double iChi = 0.0;
    double iVol = 0.0;

    int iu[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0.. 15
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, // 16.. 31
        0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, // 32.. 47
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, // 48.. 63
        0, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, // 64.. 79
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, // 80.. 95
        -1, 0, -1, 0, -1, 0, 0, 0, -2, 0, -1, 0, -1, 0, 0, 0, // 96..111
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, // 112..127
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128..143
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, // 144..159
        0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, // 160..175
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, // 176..191
        0, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, // 192..207
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, // 208..223
        -1, 0, -1, 0, -1, 0, 0, 0, -2, 0, -1, 0, -1, 0, 0, 0, // 224..239
        0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0 // 240..255
    };

    //#pragma omp parallel for reduction (+ : iChi, iVol)
    for (l = 0; l <= UCHAR_MAX; l++) {
        iChi += iu[l] * h[l];
        iVol += h[l];
    }

    return iChi / (iVol * Delta[0] * Delta[1] * Delta[2]);
}

/* returns the specific integral of total curvature K_V from the vector h[0..255]
of absolute freqencies, the lattice distances Delta[0..2] are input*/
double euler_plus(double *h, double *Delta) {
    int l;
    double iChi = 0.0;
    double iVol = 0.0;

    int iv[256] = {0, 3, 3, 0, 3, 0, -6, -3, 3, -6, 0, -3, 0, -3, -3, 0, // 0.. 15
        3, 0, -6, -3, -6, -3, -3, -6, -12, -8, -8, -6, -8, -6, 0, -3, // 16.. 31
        3, -6, 0, -3, -12, -8, -8, -6, -6, -3, -3, -6, -8, 0, -6, -3, // 32.. 47
        0, -3, -3, 0, -8, -6, 0, -3, -8, 0, -6, -3, 0, 3, 3, 0, // 48.. 63
        3, -6, -12, -8, 0, -3, -8, -6, -6, -3, -8, 0, -3, -6, -6, -3, // 64.. 79
        0, -3, -8, -6, -3, 0, 0, -3, -8, 0, 0, 3, -6, -3, 3, 0, // 80.. 95
        -6, -3, -8, 0, -8, 0, 0, 3, -3, 12, 0, 9, 0, 9, 3, 6, // 96..111
        -3, -6, -6, -3, -6, -3, 3, 0, 0, 9, 3, 6, 3, 6, 6, 3, // 112..127
        3, -12, -6, -8, -6, -8, -3, 0, 0, -8, -3, -6, -3, -6, -6, -3, // 128..143
        -6, -8, -3, 0, -3, 0, 12, 9, -8, 0, 0, 3, 0, 3, 9, 6, // 144..159
        0, -8, -3, -6, -8, 0, 0, 3, -3, 0, 0, -3, -6, 3, -3, 0, // 160..175
        -3, -6, -6, -3, 0, 3, 9, 6, -6, 3, -3, 0, 3, 6, 6, 3, // 176..191
        0, -8, -8, 0, -3, -6, 0, 3, -3, 0, -6, 3, 0, -3, -3, 0, // 192..207
        -3, -6, 0, 3, -6, -3, 9, 6, -6, 3, 3, 6, -3, 0, 6, 3, // 208..223
        -3, 0, -6, 3, -6, 3, 3, 6, -6, 9, -3, 6, -3, 6, 0, 3, // 224..239
        0, -3, -3, 0, -3, 0, 6, 3, -3, 6, 0, 3, 0, 3, 3, 0 // 240..255
    };

    for (l = 0; l <= UCHAR_MAX; l++) {
        iChi += iv[l] * h[l];
        iVol += h[l];
    }

    return M_PI / 6 * iChi / (iVol * Delta[0] * Delta[1] * Delta[2]);
}

int p3dBasicAnalysis(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        struct BasicStats* out_stats, // OUT: Basic characteristics
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize, // IN: voxel resolution
        int (*wr_log)(const char*, ...)
        ) {

    // Allocate memory for temporary variables:
    double* h = (double *) calloc((UCHAR_MAX + 1), sizeof (double));
    double* delta = (double *) malloc(3 * sizeof (double));
    delta[0] = voxelsize;
    delta[1] = voxelsize;
    delta[2] = voxelsize;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dBasicAnalysis");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing basic analysis...");
        wr_log("\tAdopted voxelsize: %0.6f mm.", voxelsize);
    }


    // Compute histogram:
    ghist(in_im, h, dimx, dimy, dimz);

    // Compute density:
    out_stats->Vv = volfrac(h);

    // Compute specific surface volume:
    out_stats->Sv = specsurf(h, delta);

    // Compute the specific integral of mean curvature:
    out_stats->Mv = specimc(h, delta);

    // Compute the euler characteristic (related to specific 
    // integral of total curvature):   
    out_stats->Cv = euler(h, delta);


    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tDensity (Vv): %0.3f [-].", out_stats->Vv);
        wr_log("\tSpecific Surface Area (Sv): %0.3f [mm^-1].", out_stats->Sv);
        wr_log("\tIntegral of Mean Curvature (Mv): %0.3f [mm^-2].", out_stats->Mv);
        wr_log("\tEuler characteristic (Cv): %0.3f [mm^-3].", out_stats->Cv);
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Basic analysis computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:
    if (h != NULL) free(h);
    if (delta != NULL) free(delta);

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    // Release resources:
    if (h != NULL) free(h);
    if (delta != NULL) free(delta);

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}
