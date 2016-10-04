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

int p3dFrom16To8(
        unsigned short* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        unsigned short min,
        unsigned short max,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    double new_min = 0.0;
    double new_max = UCHAR_MAX * 1.0;
    double tmpval;
    int i;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Rescaling and converting image from 16-bit to 8-bit format...");
        wr_log("\tMin/max values to rescale into [0,255] range: [%d.%d].", min, max);
    }

#pragma omp parallel for private (tmpval)
    for (i = 0; i < (dimx * dimy * dimz); i++) {
        tmpval = (in_im[i] - min) / ((max - min) * 1.0)*(new_max - new_min);
        if (tmpval > new_max) tmpval = new_max;
        if (tmpval < new_min) tmpval = new_min;
        out_im[i] = (unsigned char) (tmpval);
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Image rescaled and converted successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());

    }

    // Return OK:
    return P3D_SUCCESS;
}