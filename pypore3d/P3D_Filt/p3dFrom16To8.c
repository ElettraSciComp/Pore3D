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
	) 
{
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