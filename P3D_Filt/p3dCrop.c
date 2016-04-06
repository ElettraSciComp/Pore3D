#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

int p3dCrop2D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx, // ncols
        const int dimy, // nrows
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int a_dimx, a_dimy;
    int i, j;


    // Compute dimensions of padded REV:
    a_dimx = dimx - size * 2;
    a_dimy = dimy - size * 2;

    // Copy original (internal) values:
    for (j = 0; j < a_dimy; j++)
        for (i = 0; i < a_dimx; i++)
            out_im[ I2(i, j, a_dimx) ] = in_im[ I2(i + size, j + size, dimx) ];

    return P3D_SUCCESS;

}

int p3dCrop2D_16(
        unsigned short* in_rev,
        unsigned short* out_rev,
        const int dimx, // ncols
        const int dimy, // nrows
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int a_dimx, a_dimy;
    int i, j;


    // Compute dimensions of padded REV:
    a_dimx = dimx - size * 2;
    a_dimy = dimy - size * 2;

    // Copy original (internal) values:
    for (j = 0; j < a_dimy; j++)
        for (i = 0; i < a_dimx; i++)
            out_rev[ I2(i, j, a_dimx) ] = in_rev[ I2(i + size, j + size, dimx) ];

    return P3D_SUCCESS;
}

int p3dCrop3D_8(
        unsigned char* in_rev,
        unsigned char* out_rev,
        const int dimx, // ncols
        const int dimy, // nrows
        const int dimz, // nplanes
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int a_dimx, a_dimy, a_dimz;
    int i, j, k;


    // Compute dimensions of padded REV:
    a_dimx = dimx - size * 2;
    a_dimy = dimy - size * 2;
    a_dimz = dimz - size * 2;

    // Copy original (internal) values:
    for (k = 0; k < a_dimz; k++)
        for (j = 0; j < a_dimy; j++)
            for (i = 0; i < a_dimx; i++)
                out_rev[ I(i, j, k, a_dimx, a_dimy) ] =
                    in_rev[ I(i + size, j + size, k + size, dimx, dimy) ];

    // Return OK:
    return P3D_SUCCESS;
}

int p3dCrop3D_16(
        unsigned short* in_rev,
        unsigned short* out_rev,
        const int dimx, // ncols
        const int dimy, // nrows
        const int dimz, // nplanes
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    int a_dimx, a_dimy, a_dimz;
    int i, j, k;


    // Compute dimensions of padded REV:
    a_dimx = dimx - size * 2;
    a_dimy = dimy - size * 2;
    a_dimz = dimz - size * 2;

    // Copy original (internal) values:
    for (k = 0; k < a_dimz; k++)
        for (j = 0; j < a_dimy; j++)
            for (i = 0; i < a_dimx; i++)
                out_rev[ I(i, j, k, a_dimx, a_dimy) ] =
                    in_rev[ I(i + size, j + size, k + size, dimx, dimy) ];

    // Return OK:
    return P3D_SUCCESS;
}

