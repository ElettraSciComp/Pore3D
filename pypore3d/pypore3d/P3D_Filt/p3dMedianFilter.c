#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"




int p3dMedianFilter3D_8(
	unsigned char* in_im,
	unsigned char* out_im,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	// Padded input and related dims:
	unsigned char* tmp_im;

	int a_dimx, a_dimy, a_dimz;
	int i, j, k;
	int x, y, z;
	int pr = 0;

	// Variables for computing kernel:
	int a_rad, ct;

	// Temporary array:
	unsigned int* hist;
	unsigned int sum;


	// Start tracking computational time:
	if (wr_log != NULL) {
		p3dResetStartTime();
		wr_log("Pore3D - Applying median filter...");
		wr_log("\tKernel size: %d.", size);
	}


	// Init variables:
	a_rad = size / 2; // integer division   

	// Compute dimensions of padded REV:
	a_dimx = dimx + a_rad * 2;
	a_dimy = dimy + a_rad * 2;
	a_dimz = dimz + a_rad * 2;

	// Get the replicate padded input:
	P3D_MEM_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
	P3D_TRY(p3dReplicatePadding3D_8(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));


	// Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, ct, sum, hist) reduction( + : pr)
	for (k = a_rad; k < (a_dimz - a_rad); k++) {
		for (j = a_rad; j < (a_dimy - a_rad); j++) {
			// Allocate and initialize to zero kernel histogram:
			hist = (unsigned int*) calloc((UCHAR_MAX + 1), sizeof (unsigned int));

			// Compute histogram for first step:
			for (z = (k - a_rad); z <= (k + a_rad); z++) {
				for (y = (j - a_rad); y <= (j + a_rad); y++) {
					for (x = 0; x <= (2 * a_rad); x++) {
						hist[tmp_im[ I(x, y, z, a_dimx, a_dimy) ]]++;
					}
				}
			}


			// Compute median:
			ct = -1;
			sum = 0;
			while (sum <= ((unsigned int) ((size * size * size) / 2))) {
				sum += hist[++ct];
			}

			// Set out voxel with the median:
			out_im[ I(0, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

			// Increase progress counter:
			pr++;


			// Scan along x dimension:
			for (i = (a_rad + 1); i < (a_dimx - a_rad); i++) {
				// Update "sliding" histogram:
				for (z = (k - a_rad); z <= (k + a_rad); z++) {
					for (y = (j - a_rad); y <= (j + a_rad); y++) {
						hist[tmp_im[ I(i - a_rad - 1, y, z, a_dimx, a_dimy) ]]--;
						hist[tmp_im[ I(i + a_rad, y, z, a_dimx, a_dimy) ]]++;
					}
				}

				// Compute median:
				ct = -1;
				sum = 0;
				while (sum <= ( (unsigned int) (((size * size * size) / 2)))) {
					sum += hist[++ct];
				}

				// Set out voxel with the median:
				out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

				// Increase progress counter:
				pr++;
			}

			// Clear histogram:
			if (hist != NULL) free(hist);
		}

		// Update any progress bar:
		if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy * dimz)*100 + 0.5));
	}


	// Print elapsed time (if required):
	if (wr_log != NULL) {
		wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}

	// Release resources:
	if (tmp_im != NULL) free(tmp_im);

	// Return success:
	return P3D_SUCCESS;


MEM_ERROR:

	if (wr_log != NULL) {
		wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
	}

	// Release resources:	
	if (tmp_im != NULL) free(tmp_im);

	// Return error:
	return P3D_ERROR;  

}


int p3dMedianFilter3D_16(
	unsigned short* in_im,
	unsigned short* out_im,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	)
{
	// Padded input and related dims:
	unsigned short* tmp_im;

	int a_dimx, a_dimy, a_dimz;
	int i, j, k;
	int x, y, z;
	int pr = 0;

	// Variables for computing kernel:
	int a_rad, ct;

	// Temporary array:
	unsigned int* hist;
	unsigned int sum;


	// Start tracking computational time:
	if (wr_log != NULL) {
		p3dResetStartTime();
		wr_log("Pore3D - Applying median filter...");
		wr_log("\tKernel size: %d.", size);
	}


	// Init variables:
	a_rad = size / 2; // integer division   

	// Compute dimensions of padded REV:
	a_dimx = dimx + a_rad * 2;
	a_dimy = dimy + a_rad * 2;
	a_dimz = dimz + a_rad * 2;

	// Get the replicate padded input:
	P3D_MEM_TRY(tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned short)));
	P3D_TRY(p3dReplicatePadding3D_16(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));


	// Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, ct, sum, hist) reduction( + : pr)
	for (k = a_rad; k < (a_dimz - a_rad); k++) {
		for (j = a_rad; j < (a_dimy - a_rad); j++) {
			// Allocate and initialize to zero kernel histogram:
			hist = (unsigned int*) calloc((USHRT_MAX + 1), sizeof (unsigned int));

			// Compute histogram for first step:
			for (z = (k - a_rad); z <= (k + a_rad); z++) {
				for (y = (j - a_rad); y <= (j + a_rad); y++) {
					for (x = 0; x <= (2 * a_rad); x++) {
						hist[tmp_im[ I(x, y, z, a_dimx, a_dimy) ]]++;
					}
				}
			}


			// Compute median:
			ct = -1;
			sum = 0;
			while (sum <= ((unsigned int) ((size * size * size) / 2))) {
				sum += hist[++ct];
			}

			// Set out voxel with the median:
			out_im[ I(0, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

			// Increase progress counter:
			pr++;


			// Scan along x dimension:
			for (i = (a_rad + 1); i < (a_dimx - a_rad); i++) {
				// Update "sliding" histogram:
				for (z = (k - a_rad); z <= (k + a_rad); z++) {
					for (y = (j - a_rad); y <= (j + a_rad); y++) {
						hist[tmp_im[ I(i - a_rad - 1, y, z, a_dimx, a_dimy) ]]--;
						hist[tmp_im[ I(i + a_rad, y, z, a_dimx, a_dimy) ]]++;
					}
				}

				// Compute median:
				ct = -1;
				sum = 0;
				while (sum <= ( (unsigned int) (((size * size * size) / 2)))) {
					sum += hist[++ct];
				}

				// Set out voxel with the median:
				out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = ct;

				// Increase progress counter:
				pr++;
			}

			// Clear histogram:
			if (hist != NULL) free(hist);
		}

		// Update any progress bar:
		if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy * dimz)*100 + 0.5));
	}


	// Print elapsed time (if required):
	if (wr_log != NULL) {
		wr_log("Pore3D - Median filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}

	// Release resources:
	if (tmp_im != NULL) free(tmp_im);

	// Return success:
	return P3D_SUCCESS;


MEM_ERROR:

	if (wr_log != NULL) {
		wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
	}

	// Release resources:	
	if (tmp_im != NULL) free(tmp_im);

	// Return error:
	return P3D_ERROR;  

}