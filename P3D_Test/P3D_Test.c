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