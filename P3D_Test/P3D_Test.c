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
#include <stdarg.h>

#include "p3dBlob.h"
#include "p3dFilt.h"
#include "p3dSkel.h"

int customPrint_nolf(const char *msg, ...) {
	va_list fmtargs;
	char buffer[1024];	

	va_start(fmtargs, msg);
	vsnprintf(buffer, sizeof (buffer) - 1, msg, fmtargs);
	va_end(fmtargs);

	return printf("%s", buffer);
}

int customPrint(const char *msg, ...) {
	va_list fmtargs;
	char buffer[1024];

	va_start(fmtargs, msg);
	vsnprintf(buffer, sizeof (buffer) - 1, msg, fmtargs);
	va_end(fmtargs);

	return printf("%s\n", buffer);
}

int customProgress(const int msg, ...) {

	return customPrint_nolf("\tProgress: %d.\r", msg);

}

int main(int argc, char* argv[]) {
	const int dimx = atoi(argv[3]);
	const int dimy = atoi(argv[4]);
	const int dimz = atoi(argv[5]);	


	// Allocate memory for input and output images:
	unsigned char* in_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char));
	unsigned char* out_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char));
	
	// Read the slices:
	p3dReadRaw8 ( argv[1], in_im, dimx, dimy, dimz, customPrint, NULL);

	// Filtering:
	p3dBilateralFilter3D_8(in_im, out_im, dimx, dimy, dimz, 3, 1.0, 10.0, 50, customPrint, customProgress );
	
	// Output:
	p3dWriteRaw8(out_im, argv[2], dimx, dimy, dimz, customPrint, NULL);


	// Release resources:	
	if (in_im != NULL) free(in_im);
	if (out_im != NULL) free(out_im);
}