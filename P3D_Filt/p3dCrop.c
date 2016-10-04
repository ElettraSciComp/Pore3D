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

#include "p3dFilt.h"
#include "p3dTime.h"

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

