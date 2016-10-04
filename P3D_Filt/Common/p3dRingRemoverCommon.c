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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "../p3dFilt.h"

#define PI 3.1415926535897932384626433

/*int _isPath(char* path)
{
    struct stat statbuf;
    
    if ( stat(path, &statbuf) == -1 )
        return -1;
    else
        return S_ISDIR(statbuf.st_mode);
}*/

double _cubic(double x) {
    double z = 0;
    double a = 1; // a = 0.5 for less sharpening

    // Absolute value:
    if (x < 0)
        x = -x;

    // Apply polinomial definition:
    if (x < 1)
        z = (-a + 2) * x * x * x + (a - 3) * x * x + 1;
    else if (x < 2)
        z = -a * x * x * x + 5 * a * x * x - 8 * a * x + 4 * a;

    return z;
}
// Should be thread-safe

unsigned char p3dBicubicInterpolation_8(
        unsigned char* in_im,
        int dimx,
        int dimy,
        double x0,
        double y0
        ) {
    int u0 = (int) x0;
    int v0 = (int) y0;
    int i, j;
    int v;
    double p;
    int u, tmpu, tmpv;

    double q = 0;

    for (j = -1; j <= 2; j++) {
        v = v0 + j;
        p = 0;
        for (i = -1; i <= 2; i++) {
            u = u0 + i;
            tmpu = u;
            tmpv = v;

            // Replicate padding for X and Y direction:
            if (u < 0)
                tmpu = 0;
            if (u > ((int) dimx - 1))
                tmpu = (int) dimx - 1;
            if (v < 0)
                tmpv = 0;
            if (v > ((int) dimy - 1))
                tmpv = (int) dimy - 1;

            p = p + in_im[ I2(tmpu, tmpv, dimx) ] * _cubic(x0 - u);
        }
        q = q + p * _cubic(y0 - v);
    }

    q = q + 0.5;

    // Check intensity boundary:
    if (q < 0.0) q = 0.0;
    if (q >= UCHAR_MAX) q = UCHAR_MAX;

    return (unsigned char) q;
}



// Should be thread-safe

unsigned short p3dBicubicInterpolation_16(
        unsigned short* in_im,
        int dimx,
        int dimy,
        double x0,
        double y0
        ) {
    int u0 = (int) x0;
    int v0 = (int) y0;

    int i, j;
    int u, tmpu, tmpv;
    int v;
    double p;

    double q = 0;

    for (j = -1; j <= 2; j++) {
        v = v0 + j;
        p = 0;
        for (i = -1; i <= 2; i++) {
            u = u0 + i;
            tmpu = u;
            tmpv = v;

            // Replicate padding for X and Y direction:
            if (u < 0)
                tmpu = 0;
            if (u > ((int) dimx - 1))
                tmpu = (int) dimx - 1;
            if (v < 0)
                tmpv = 0;
            if (v > ((int) dimy - 1))
                tmpv = (int) dimy - 1;

            p = p + in_im[ I2(tmpu, tmpv, dimx) ] * _cubic(x0 - u);
        }
        q = q + p * _cubic(y0 - v);
    }

    q = q + 0.5;

    // Check intensity boundary:
    if (q < 0.0) q = 0.0;
    if (q >= USHRT_MAX) q = USHRT_MAX;

    return (unsigned short) q;
}



// Should be thread-safe

unsigned char p3dBicubicInterpolation_p2c_8(
        unsigned char* in_im,
        int dimx,
        int dimy,
        double x0,
        double y0
        ) {
    int u0 = (int) x0;
    int v0 = (int) y0;
    int u, tmpu, tmpv;
    int v;
    double p;
    int i, j;

    double q = 0;

    for (j = -1; j <= 2; j++) {
        v = v0 + j;
        p = 0;
        for (i = -1; i <= 2; i++) {
            u = u0 + i;
            tmpu = u;
            tmpv = v;

            // Replicate padding for X direction and circular padding for Y
            // direction:
            if (u < 0)
                tmpu = 0;
            if (u > ((int) dimx - 1))
                tmpu = (int) dimx - 1;
            if (v < 0)
                tmpv = v + (int) dimy;
            if (v > ((int) dimy - 1))
                tmpv = v - (int) dimy;

            p = p + in_im[ I2(tmpu, tmpv, dimx) ] * _cubic(x0 - u);
        }
        q = q + p * _cubic(y0 - v);
    }

    q = q + 0.5;

    // Check intensity boundary:
    if (q < 0.0) q = 0.0;
    if (q >= UCHAR_MAX) q = UCHAR_MAX;

    return (unsigned char) q;
}


// Should be thread-safe

unsigned short p3dBicubicInterpolation_p2c_16(
        unsigned short* in_im,
        int dimx,
        int dimy,
        double x0,
        double y0
        ) {
    int u0 = (int) x0;
    int v0 = (int) y0;
    int i, j;
    int u, v, tmpu, tmpv;
    double p;

    double q = 0;

    for (j = -1; j <= 2; j++) {
        v = v0 + j;
        p = 0;
        for (i = -1; i <= 2; i++) {
            u = u0 + i;
            tmpu = u;
            tmpv = v;

            // Replicate padding for X direction and circular padding for Y
            // direction:
            if (u < 0)
                tmpu = 0;
            if (u > ((int) dimx - 1))
                tmpu = (int) dimx - 1;
            if (v < 0)
                tmpv = v + (int) dimy;
            if (v > ((int) dimy - 1))
                tmpv = v - (int) dimy;

            p = p + in_im[ I2(tmpu, tmpv, dimx) ] * _cubic(x0 - u);
        }
        q = q + p * _cubic(y0 - v);
    }

    q = q + 0.5;

    // Check intensity boundary:
    if (q < 0.0) q = 0.0;
    if (q >= USHRT_MAX) q = USHRT_MAX;

    return (unsigned short) q;
}


// A square image with dimensions POLARX x POLARX is returned. Due to the fact
// the POLARX is not known a priori, memory image is allocated within this
// procedure so parameter OUT_IM should be passed as reference. The out parameter
// POLARX should be used for further handling of OUT_IM.

void p3dCartesian2polar_8(
        unsigned char* in_im, // IN: input image in cartesian coordinates
        unsigned char** out_im, // OUT: output image in polar coordinates
        const int dimx, // IN: width of input image
        const int dimy, // IN: heigth of input image
        const double centerX, // IN: X coordinate for center of polar transform
        const double centerY, // IN: Y coordinate for center of polar transform
        const double precision,
        int* polarX // OUT: side of the square image returned as output
        ) {
    int i, j;
    double r, phi;
    double x, y;

    double r1, r2, r3, r4;


    // Find the greatest radius:

    // Top-Left Corner (0,0):
    r1 = sqrt((centerX - 0)*(centerX - 0) + (centerY - 0)*(centerY - 0));
    // Top-Right Corner (widthInitial, 0):
    r2 = sqrt((centerX - dimx)*(centerX - dimx) + (centerY - 0)*(centerY - 0));
    // Bottom-Left Corner (0, heightInitial):
    r3 = sqrt((centerX - 0)*(centerX - 0) + (centerY - dimy)*(centerY - dimy));
    // Bottom-Right Corner (widthInitial , heightInitial):
    r4 = sqrt((centerX - dimx)*(centerX - dimx) + (centerY - dimy)*(centerY - dimy));


    // The greatest radius is the semi-width of the polar image:
    *polarX = (int) (precision * (MAX(MAX(MAX(r1, r2), r3), r4) + 0.5));

    // Now that we know dimensions, allocate memory:
    (*out_im) = (unsigned char*) calloc((*polarX)*(*polarX), sizeof (unsigned char));


    // Fill the polar grid:
#pragma omp parallel for private(i, r, phi, x, y)
    for (j = 0; j < (*polarX); j++)
        for (i = 0; i < (*polarX); i++) {
            // Get [r,phi]:
            r = (double) (i);
            phi = ((double) (j) / (*polarX)) * PI * 2;

            // Get [x,y] from [r,phi]:
            x = r * cos(phi) + centerX;
            y = r * sin(phi) + centerY;


            // Perform bicubic interpolation:
            (*out_im)[ I2(i, j, (*polarX)) ] = p3dBicubicInterpolation_8(in_im, dimx, dimy, x, y);

        }

}

void p3dCartesian2polar_16(
        unsigned short* in_im, // IN: input image in cartesian coordinates
        unsigned short** out_im, // OUT: output image in polar coordinates
        const int dimx, // IN: width of input image
        const int dimy, // IN: heigth of input image
        const double centerX, // IN: X coordinate for center of polar transform
        const double centerY, // IN: Y coordinate for center of polar transform
        const double precision,
        int* polarX // OUT: side of the square image returned as output
        ) {
    int i, j;
    double r, phi;
    double x, y;

    double r1, r2, r3, r4;


    // Find the greatest radius:

    // Top-Left Corner (0,0):
    r1 = sqrt((centerX - 0)*(centerX - 0) + (centerY - 0)*(centerY - 0));
    // Top-Right Corner (widthInitial, 0):
    r2 = sqrt((centerX - dimx)*(centerX - dimx) + (centerY - 0)*(centerY - 0));
    // Bottom-Left Corner (0, heightInitial):
    r3 = sqrt((centerX - 0)*(centerX - 0) + (centerY - dimy)*(centerY - dimy));
    // Bottom-Right Corner (widthInitial , heightInitial):
    r4 = sqrt((centerX - dimx)*(centerX - dimx) + (centerY - dimy)*(centerY - dimy));


    // The greatest radius is the semi-width of the polar image:
    *polarX = (int) (precision * (MAX(MAX(MAX(r1, r2), r3), r4) + 0.5));

    // Now that we know dimensions, allocate memory:
    (*out_im) = (unsigned short*) calloc((*polarX)*(*polarX), sizeof (unsigned short));


    // Fill the polar grid:
#pragma omp parallel for private(i, r, phi, x, y)
    for (j = 0; j < (*polarX); j++)
        for (i = 0; i < (*polarX); i++) {
            // Get [r,phi]:
            r = (double) (i);
            phi = ((double) (j) / (*polarX)) * PI * 2;

            // Get [x,y] from [r,phi]:
            x = r * cos(phi) + centerX;
            y = r * sin(phi) + centerY;


            // Perform bicubic interpolation:
            (*out_im)[ I2(i, j, (*polarX)) ] = p3dBicubicInterpolation_16(in_im, dimx, dimy, x, y);

        }

}



// For coherence with dual procedure CARTESIAN2POLAR memory image is allocated
// within this procedure so parameter OUT_IM should be passed as reference.

void p3dPolar2cartesian_8(
        unsigned char* in_im, // IN: input image in polar coordinates
        unsigned char** out_im, // OUT: output image in cartesian coordinates
        const int polarX, // IN: side of the square of polar image
        const double centerX, // IN: X coordinate for center of polar transform
        const double centerY, // IN: Y coordinate for center of polar transform
        const int original_dimx, // IN: width of the output cartesian image
        const int original_dimy // IN: heigth of the output cartesian image
        ) {
    double r, phi;
    int i, j;
    double x, y;

    // Allocate memory:
    (*out_im) = (unsigned char*) calloc(original_dimx*original_dimy, sizeof (unsigned char));

    // Fill the cartesian grid:
#pragma omp parallel for private(i, r, phi, x, y)
    for (j = 0; j < original_dimy; j++)
        for (i = 0; i < original_dimx; i++) {
            x = (double) (i) - centerX;
            y = (double) (j) - centerY;

            // Get [r,phi] from [x,y]:
            r = sqrt(x * x + y * y);
            phi = (atan2(y, x) * polarX) / (PI * 2); // Degree unit
            if (phi < 0)
                phi += polarX;

            // Perform bicubic interpolation:
            (*out_im)[ I2(i, j, original_dimx) ] = p3dBicubicInterpolation_p2c_8(in_im, polarX, polarX, r, phi);

        }
}

void p3dPolar2cartesian_16(
        unsigned short* in_im, // IN: input image in polar coordinates
        unsigned short** out_im, // OUT: output image in cartesian coordinates
        const int polarX, // IN: side of the square of polar image
        const double centerX, // IN: X coordinate for center of polar transform
        const double centerY, // IN: Y coordinate for center of polar transform
        const int original_dimx, // IN: width of the output cartesian image
        const int original_dimy // IN: heigth of the output cartesian image
        ) {
    double r, phi;
    int i, j;
    double x, y;

    // Allocate memory:
    (*out_im) = (unsigned short*) calloc(original_dimx*original_dimy, sizeof (unsigned short));

    // Fill the cartesian grid:
#pragma omp parallel for private(i, r, phi, x, y)
    for (j = 0; j < original_dimy; j++)
        for (i = 0; i < original_dimx; i++) {
            x = (double) (i) - centerX;
            y = (double) (j) - centerY;

            // Get [r,phi] from [x,y]:
            r = sqrt(x * x + y * y);
            phi = (atan2(y, x) * polarX) / (PI * 2); // Degree unit
            if (phi < 0)
                phi += polarX;

            // Perform bicubic interpolation:
            (*out_im)[ I2(i, j, original_dimx) ] = p3dBicubicInterpolation_p2c_16(in_im, polarX, polarX, r, phi);

        }
}