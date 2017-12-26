#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"

int p3dZeroPadding2D_8(
	unsigned char* in_im,
	unsigned char* out_im,
	const int dimx, // ncols
	const int dimy, // nrows
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy;
	int i, j;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;

	// Set to zero all values:
	memset(out_im, 0, a_dimx * a_dimy * sizeof (unsigned char));

	// Copy original (internal) values:
	for (j = 0; j < dimy; j++) {
		for (i = 0; i < dimx; i++) {
			out_im[ I2(i + size, j + size, a_dimx) ] = in_im[ I2(i, j, dimx) ];
		}
	}

	return P3D_SUCCESS;
}

int p3dZeroPadding2D_16(
	unsigned short* in_im,
	unsigned short* out_im,
	const int dimx, // ncols
	const int dimy, // nrows
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy;
	int i, j;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;

	// Set to zero all values:
	memset(out_im, 0, a_dimx * a_dimy * sizeof (unsigned short));

	// Copy original (internal) values:
	for (j = 0; j < dimy; j++) {
		for (i = 0; i < dimx; i++) {
			out_im[ I2(i + size, j + size, a_dimx) ] = in_im[ I2(i, j, dimx) ];
		}
	}

	return P3D_SUCCESS;
}

int p3dZeroPadding3D_8(
	unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;

	// Set to zero all values:
	memset(out_rev, 0, a_dimx * a_dimy * a_dimz * sizeof (unsigned char));


	// Copy original (internal) values:
	for (k = 0; k < dimz; k++) {
		for (j = 0; j < dimy; j++) {
			for (i = 0; i < dimx; i++) {
				out_rev[ I(i + size, j + size, k + size, a_dimx, a_dimy) ] =
					in_rev[ I(i, j, k, dimx, dimy) ];
			}
		}
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dZeroPadding3D_16(
	unsigned short* in_rev,
	unsigned short* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;

	// Set to zero all values:
	memset(out_rev, 0, a_dimx * a_dimy * a_dimz * sizeof (unsigned short));


	// Copy original (internal) values:
	for (k = 0; k < dimz; k++) {
		for (j = 0; j < dimy; j++) {
			for (i = 0; i < dimx; i++) {
				out_rev[ I(i + size, j + size, k + size, a_dimx, a_dimy) ] =
					in_rev[ I(i, j, k, dimx, dimy) ];
			}
		}
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dReplicatePadding2D_8(
	unsigned char* in_im,
	unsigned char* out_im,
	const int dimx, // ncols
	const int dimy, // nrows
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy;
	int i, j, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;


	// Perform first zero padding:
	p3dZeroPadding2D_8(in_im, out_im, dimx, dimy, size, NULL, NULL);


	// Replicate border values:
	for (ct = size; ct > 0; ct--) {
		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_im[ I2(i, ct - 1, a_dimx) ] =
				out_im[ I2(i, ct, a_dimx) ];
			out_im[ I2(i, a_dimy - ct, a_dimx) ] =
				out_im[ I2(i, a_dimy - 1 - ct, a_dimx) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_im[ I2(ct - 1, j, a_dimx) ] =
				out_im[ I2(ct, j, a_dimx) ];
			out_im[ I2(a_dimx - ct, j, a_dimx) ] =
				out_im[ I2(a_dimx - 1 - ct, j, a_dimx) ];
		}

		// Corners:

		out_im[ I2(ct - 1, ct - 1, a_dimx)] =
			out_im[ I2(ct, ct, a_dimx) ];

		out_im[ I2(a_dimx - ct, ct - 1, a_dimx)] =
			out_im[ I2(a_dimx - 1 - ct, ct, a_dimx) ];

		out_im[ I2(ct - 1, a_dimy - ct, a_dimx)] =
			out_im[ I2(ct, a_dimy - 1 - ct, a_dimx) ];

		out_im[ I2(a_dimx - ct, a_dimy - ct, a_dimx)] =
			out_im[ I2(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimx) ];

	}

	return P3D_SUCCESS;
}

int p3dReplicatePadding2D_16(
	unsigned short* in_im,
	unsigned short* out_im,
	const int dimx, // ncols
	const int dimy, // nrows
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy;
	int i, j, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;


	// Perform first zero padding:
	p3dZeroPadding2D_16(in_im, out_im, dimx, dimy, size, NULL, NULL);


	// Replicate border values:
	for (ct = size; ct > 0; ct--) {
		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_im[ I2(i, ct - 1, a_dimx) ] =
				out_im[ I2(i, ct, a_dimx) ];
			out_im[ I2(i, a_dimy - ct, a_dimx) ] =
				out_im[ I2(i, a_dimy - 1 - ct, a_dimx) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_im[ I2(ct - 1, j, a_dimx) ] =
				out_im[ I2(ct, j, a_dimx) ];
			out_im[ I2(a_dimx - ct, j, a_dimx) ] =
				out_im[ I2(a_dimx - 1 - ct, j, a_dimx) ];
		}

		// Corners:

		out_im[ I2(ct - 1, ct - 1, a_dimx)] =
			out_im[ I2(ct, ct, a_dimx) ];

		out_im[ I2(a_dimx - ct, ct - 1, a_dimx)] =
			out_im[ I2(a_dimx - 1 - ct, ct, a_dimx) ];

		out_im[ I2(ct - 1, a_dimy - ct, a_dimx)] =
			out_im[ I2(ct, a_dimy - 1 - ct, a_dimx) ];

		out_im[ I2(a_dimx - ct, a_dimy - ct, a_dimx)] =
			out_im[ I2(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimx) ];

	}

	return P3D_SUCCESS;

}

int p3dReplicatePadding3D_8(
	unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;


	// Perform first zero padding:
	p3dZeroPadding3D_8(in_rev, out_rev, dimx, dimy, dimz, size, NULL, NULL);


	// Replicate border values:
	for (ct = size; ct > 0; ct--) {
		// Faces:

		for (i = ct; i < (a_dimx - ct); i++) {
			for (j = ct; j < (a_dimy - ct); j++) {
				out_rev[ I(i, j, ct - 1, a_dimx, a_dimy) ] =
					out_rev[ I(i, j, ct, a_dimx, a_dimy) ];
				out_rev[ I(i, j, a_dimz - ct, a_dimx, a_dimy)] =
					out_rev[ I(i, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			}
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(ct - 1, j, k, a_dimx, a_dimy) ] =
					out_rev[ I(ct, j, k, a_dimx, a_dimy) ];
				out_rev[ I(a_dimx - ct, j, k, a_dimx, a_dimy)] =
					out_rev[ I(a_dimx - 1 - ct, j, k, a_dimx, a_dimy) ];
			}
		}

		for (i = ct; i < (a_dimx - ct); i++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(i, ct - 1, k, a_dimx, a_dimy) ] =
					out_rev[ I(i, ct, k, a_dimx, a_dimy) ];
				out_rev[ I(i, a_dimy - ct, k, a_dimx, a_dimy)] =
					out_rev[ I(i, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			}
		}

		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_rev[ I(i, ct - 1, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_rev[ I(ct - 1, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(a_dimx - 1 - ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (k = ct; k < (a_dimz - ct); k++) {
			out_rev[ I(ct - 1, ct - 1, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
		}

		// Corners:

		out_rev[ I(ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


		out_rev[ I(a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dReplicatePadding3D_16(
	unsigned short* in_rev,
	unsigned short* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;


	// Perform first zero padding:
	p3dZeroPadding3D_16(in_rev, out_rev, dimx, dimy, dimz, size, NULL, NULL);


	// Replicate border values:
	for (ct = size; ct > 0; ct--) {
		// Faces:

		for (i = ct; i < (a_dimx - ct); i++) {
			for (j = ct; j < (a_dimy - ct); j++) {
				out_rev[ I(i, j, ct - 1, a_dimx, a_dimy) ] =
					out_rev[ I(i, j, ct, a_dimx, a_dimy) ];
				out_rev[ I(i, j, a_dimz - ct, a_dimx, a_dimy)] =
					out_rev[ I(i, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			}
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(ct - 1, j, k, a_dimx, a_dimy) ] =
					out_rev[ I(ct, j, k, a_dimx, a_dimy) ];
				out_rev[ I(a_dimx - ct, j, k, a_dimx, a_dimy)] =
					out_rev[ I(a_dimx - 1 - ct, j, k, a_dimx, a_dimy) ];
			}
		}

		for (i = ct; i < (a_dimx - ct); i++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(i, ct - 1, k, a_dimx, a_dimy) ] =
					out_rev[ I(i, ct, k, a_dimx, a_dimy) ];
				out_rev[ I(i, a_dimy - ct, k, a_dimx, a_dimy)] =
					out_rev[ I(i, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			}
		}

		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_rev[ I(i, ct - 1, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_rev[ I(ct - 1, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(a_dimx - 1 - ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (k = ct; k < (a_dimz - ct); k++) {
			out_rev[ I(ct - 1, ct - 1, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
		}

		// Corners:
		out_rev[ I(ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


		out_rev[ I(a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


	}

	// Return OK:
	return P3D_SUCCESS;
}


int _p3dZeroPadding3D_float(
	float* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;

	// Set to zero all values:
	memset(out_rev, 0, a_dimx * a_dimy * a_dimz * sizeof (float));


	// Copy original (internal) values:
	#pragma omp parallel for private(i, j)
	for (k = 0; k < dimz; k++) {
		for (j = 0; j < dimy; j++) {
			for (i = 0; i < dimx; i++) {
				out_rev[ I(i + size, j + size, k + size, a_dimx, a_dimy) ] =
					(float) (in_rev[ I(i, j, k, dimx, dimy) ]);
			}
		}
	}

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dReplicatePadding3D_float(
	float* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;


	// Perform first zero padding:
	_p3dZeroPadding3D_float(in_rev, out_rev, dimx, dimy, dimz, size);


	// Replicate border values:
	#pragma omp parallel for private(i, j, k)
	for (ct = size; ct > 0; ct--) {
		// Faces:

		for (i = ct; i < (a_dimx - ct); i++) {
			for (j = ct; j < (a_dimy - ct); j++) {
				out_rev[ I(i, j, ct - 1, a_dimx, a_dimy) ] =
					out_rev[ I(i, j, ct, a_dimx, a_dimy) ];
				out_rev[ I(i, j, a_dimz - ct, a_dimx, a_dimy)] =
					out_rev[ I(i, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			}
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(ct - 1, j, k, a_dimx, a_dimy) ] =
					out_rev[ I(ct, j, k, a_dimx, a_dimy) ];
				out_rev[ I(a_dimx - ct, j, k, a_dimx, a_dimy)] =
					out_rev[ I(a_dimx - 1 - ct, j, k, a_dimx, a_dimy) ];
			}
		}

		for (i = ct; i < (a_dimx - ct); i++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(i, ct - 1, k, a_dimx, a_dimy) ] =
					out_rev[ I(i, ct, k, a_dimx, a_dimy) ];
				out_rev[ I(i, a_dimy - ct, k, a_dimx, a_dimy)] =
					out_rev[ I(i, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			}
		}

		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_rev[ I(i, ct - 1, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(i, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_rev[ I(ct - 1, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] =
				out_rev[ I(a_dimx - 1 - ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (k = ct; k < (a_dimz - ct); k++) {
			out_rev[ I(ct - 1, ct - 1, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] =
				out_rev[ I(ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] =
				out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
		}

		// Corners:

		out_rev[ I(ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


		out_rev[ I(a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


	}

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dZeroPadding3D_uchar2float(
	unsigned char* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;

	// Set to zero all values:
	memset(out_rev, 0, a_dimx * a_dimy * a_dimz * sizeof (float));


	// Copy original (internal) values:
	#pragma omp parallel for private(i, j)
	for (k = 0; k < dimz; k++) {
		for (j = 0; j < dimy; j++) {
			for (i = 0; i < dimx; i++) {
				out_rev[ I(i + size, j + size, k + size, a_dimx, a_dimy) ] =
					(float) (in_rev[ I(i, j, k, dimx, dimy) ]);
			}
		}
	}

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dReplicatePadding3D_uchar2float(
	unsigned char* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;


	// Perform first zero padding:
	_p3dZeroPadding3D_uchar2float(in_rev, out_rev, dimx, dimy, dimz, size);


	// Replicate border values:
	#pragma omp parallel for private(i, j, k)
	for (ct = size; ct > 0; ct--) {
		// Faces:

		for (i = ct; i < (a_dimx - ct); i++) {
			for (j = ct; j < (a_dimy - ct); j++) {
				out_rev[ I(i, j, ct - 1, a_dimx, a_dimy) ] =
					(float) out_rev[ I(i, j, ct, a_dimx, a_dimy) ];
				out_rev[ I(i, j, a_dimz - ct, a_dimx, a_dimy)] =
					(float) out_rev[ I(i, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			}
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(ct - 1, j, k, a_dimx, a_dimy) ] =
					(float) out_rev[ I(ct, j, k, a_dimx, a_dimy) ];
				out_rev[ I(a_dimx - ct, j, k, a_dimx, a_dimy)] =
					(float) out_rev[ I(a_dimx - 1 - ct, j, k, a_dimx, a_dimy) ];
			}
		}

		for (i = ct; i < (a_dimx - ct); i++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(i, ct - 1, k, a_dimx, a_dimy) ] =
					(float) out_rev[ I(i, ct, k, a_dimx, a_dimy) ];
				out_rev[ I(i, a_dimy - ct, k, a_dimx, a_dimy)] =
					(float) out_rev[ I(i, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			}
		}

		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_rev[ I(i, ct - 1, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(i, ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(i, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(i, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_rev[ I(ct - 1, j, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(a_dimx - 1 - ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (k = ct; k < (a_dimz - ct); k++) {
			out_rev[ I(ct - 1, ct - 1, k, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
		}

		// Corners:

		out_rev[ I(ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


		out_rev[ I(a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


	}

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dZeroPadding3D_ushort2float(
	unsigned short* in_rev,
	float* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;

	// Set to zero all values:
	memset(out_rev, 0, a_dimx * a_dimy * a_dimz * sizeof (float));


	// Copy original (internal) values:
	#pragma omp parallel for private(i, j)
	for (k = 0; k < dimz; k++) {
		for (j = 0; j < dimy; j++) {
			for (i = 0; i < dimx; i++) {
				out_rev[ I(i + size, j + size, k + size, a_dimx, a_dimy) ] =
					(float) (in_rev[ I(i, j, k, dimx, dimy) ]);
			}
		}
	}

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dReplicatePadding3D_ushort2float(
	unsigned short* in_rev,
	float* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size
	) 
{
	int a_dimx, a_dimy, a_dimz;
	int i, j, k, ct;


	// Compute dimensions of padded REV:
	a_dimx = dimx + size * 2;
	a_dimy = dimy + size * 2;
	a_dimz = dimz + size * 2;


	// Perform first zero padding:
	_p3dZeroPadding3D_ushort2float(in_rev, out_rev, dimx, dimy, dimz, size);


	// Replicate border values:
	#pragma omp parallel for private(i, j, k)
	for (ct = size; ct > 0; ct--) {
		// Faces:

		for (i = ct; i < (a_dimx - ct); i++) {
			for (j = ct; j < (a_dimy - ct); j++) {
				out_rev[ I(i, j, ct - 1, a_dimx, a_dimy) ] =
					(float) out_rev[ I(i, j, ct, a_dimx, a_dimy) ];
				out_rev[ I(i, j, a_dimz - ct, a_dimx, a_dimy)] =
					(float) out_rev[ I(i, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			}
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(ct - 1, j, k, a_dimx, a_dimy) ] =
					(float) out_rev[ I(ct, j, k, a_dimx, a_dimy) ];
				out_rev[ I(a_dimx - ct, j, k, a_dimx, a_dimy)] =
					(float) out_rev[ I(a_dimx - 1 - ct, j, k, a_dimx, a_dimy) ];
			}
		}

		for (i = ct; i < (a_dimx - ct); i++) {
			for (k = ct; k < (a_dimz - ct); k++) {
				out_rev[ I(i, ct - 1, k, a_dimx, a_dimy) ] =
					(float) out_rev[ I(i, ct, k, a_dimx, a_dimy) ];
				out_rev[ I(i, a_dimy - ct, k, a_dimx, a_dimy)] =
					(float) out_rev[ I(i, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			}
		}

		// Edges:

		for (i = ct; i < (a_dimx - ct); i++) {
			out_rev[ I(i, ct - 1, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(i, ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(i, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(i, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];
			out_rev[ I(i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (j = ct; j < (a_dimy - ct); j++) {
			out_rev[ I(ct - 1, j, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] =
				(float) out_rev[ I(a_dimx - 1 - ct, j, ct, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy) ];
		}

		for (k = ct; k < (a_dimz - ct); k++) {
			out_rev[ I(ct - 1, ct - 1, k, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, ct, k, a_dimx, a_dimy) ];
			out_rev[ I(ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] =
				(float) out_rev[ I(ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
			out_rev[ I(a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] =
				(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy) ];
		}

		// Corners:

		out_rev[ I(ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];


		out_rev[ I(a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy) ];

		out_rev[ I(ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

		out_rev[ I(a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] =
			(float) out_rev[ I(a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy) ];

	}

	// Return OK:
	return P3D_SUCCESS;
}

