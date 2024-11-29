#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dBlob.h"
#include "p3dTime.h"

#include "Common/p3dConnectedComponentsLabeling.h"
#include "Common/p3dUtils.h"out_rev

int p3dWriteRaw8_local1(
	unsigned char* in_im,
	char* filename,
	const int dimx,
	const int dimy,
	const int dimz,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	FILE* fvol;

	// Start tracking computational time:
	if (wr_log != NULL) {
		p3dResetStartTime();
		wr_log("Pore3D - Writing 8-bit RAW file %s ...", filename);
	}

	/* Get a handler for the input file */
	if ((fvol = fopen(filename, "wb")) == NULL) {
		wr_log("Pore3D - IO error: cannot open output file %s. Program will exit.", filename);

		return 0;
	}

	/* Write raw data to file: */
	if (fwrite(in_im, sizeof (unsigned char), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
		wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

		return 0;
	}

	/* Close file handler: */
	fclose(fvol);

	// Print elapsed time (if required):
	if (wr_log != NULL) {
		wr_log("Pore3D - RAW file written successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec()) ;
	}


	// Return OK:
	return P3D_SUCCESS;
}

int p3dWriteRaw32_local1(
	unsigned int* in_im,
	char* filename,
	const int dimx,
	const int dimy,
	const int dimz,
	const int flagLittle,
	const int flagSigned,
	int (*wr_log)(const char*, ...),
	int (*wr_progress)(const int, ...)
	) 
{
	FILE* fvol;
	int* s_tmp_im = NULL;
	unsigned int* us_tmp_im = NULL;
	int ct;

	// Start tracking computational time:
	
	if (wr_log != NULL) {
		p3dResetStartTime();
		wr_log("Pore3D - Writing unsigned integer 32-bit RAW file %s ...", filename);
		if (flagSigned == P3D_TRUE)
			wr_log("\tSigned/Unsigned: Signed.");
		else
			wr_log("\tSigned/Unsigned: Unsigned.");
		if (flagLittle == P3D_TRUE)
			wr_log("\tLittle/Big Endian: Little.");
		else
			wr_log("\tLittle/Big Endian: Big.");
	}

	/* Get a handler for the input file */
	if ((fvol = fopen(filename, "wb")) == NULL) {
		wr_log("Pore3D - IO error: cannot open output file %s. Program will exit.", filename);
		in_im = NULL;

		return 0;
	}

	/* Read data signed/unsigned: */
	if (flagSigned == P3D_TRUE) {
		s_tmp_im = (int*) malloc(dimx * dimy * dimz * sizeof (int));

		/* Convert to signed: */
		for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
			if (flagLittle == P3D_FALSE) 
			{
				//s_tmp_im[ct] = _EndianSwapInt((int) (in_im[ct] - UINT_MAX / 2));
			} 
			else 
			{
				s_tmp_im[ct] = (int) (in_im[ct] - UINT_MAX / 2);
			}
		}

		/* Write raw data to file: */
		//fwrite(s_tmp_im, sizeof (short), dimx*dimy, fvol);
		if (fwrite(s_tmp_im, sizeof (int), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
			wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

			return 0;
		}

		/* Free: */
		free(s_tmp_im);
	} else {
		/* Swap endian if necessary: */
		if (flagLittle == P3D_FALSE) {
			us_tmp_im = (unsigned int*) malloc(dimx * dimy * dimz * sizeof (unsigned int));

			for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
				//us_tmp_im[ct] = _EndianSwapUnsignedInt(in_im[ct]);
			}

			//fwrite(us_tmp_im, sizeof (unsigned short), dimx*dimy, fvol);
			if (fwrite(us_tmp_im, sizeof (unsigned int), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
				wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

				return 0;
			}

			free(us_tmp_im);
		} else {
			fwrite(in_im, sizeof (unsigned int), dimx * dimy * dimz, fvol);
		}
	}

	/* Close file handler: */
	fclose(fvol);

	// Print elapsed time (if required):
	if (wr_log != NULL) {
		wr_log("Pore3D - RAW file written successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}
    

	return P3D_SUCCESS;
}


// Wrapper for the extern:



int p3dBlobLabeling_ushort(
        unsigned char* in_im,
        unsigned short* out_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz,
        const int conn,
        const int random_lbl, // Flag for random labels
        const int skip_borders,
        int (*wr_log)(const char*, ...)
        ) {
	
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing blob labeling...");
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
        if (random_lbl == P3D_TRUE)
            wr_log("\tRandom labels used. ");
        if (skip_borders == P3D_TRUE)
            wr_log("\tBorders skipped. ");
    }


    P3D_TRY( p3dConnectedComponentsLabeling_ushort(in_im, out_im, NULL, NULL, NULL, dimx, dimy, dimz,
            conn, random_lbl, skip_borders));

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Blob labeling performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    if (filename!=NULL)
        p3dWriteRaw8_local1(out_im, filename, dimx, dimy,dimz,0,0,NULL,NULL);

    // Return OK:
    return P3D_ERROR;    
}

int p3dBlobLabeling_uint(
        unsigned char* in_im,
        unsigned int* out_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz,
        const int conn,
        const int random_lbl, // Flag for random labels
        const int skip_borders,
        int (*wr_log)(const char*, ...)
        ) {

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing blob labeling...");
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
        if (random_lbl == P3D_TRUE)
            wr_log("\tRandom labels used. ");
        if (skip_borders == P3D_TRUE)
            wr_log("\tBorders skipped. ");
    }

    P3D_MEM_TRY(out_im = (unsigned int*) malloc(dimx * dimy * dimz * sizeof (unsigned int)));
    P3D_TRY( p3dConnectedComponentsLabeling_uint(in_im, out_im, NULL, NULL, NULL, dimx, dimy, dimz,conn, random_lbl, skip_borders) );

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Blob labeling performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }
    
    if (filename!=NULL)
        p3dWriteRaw32_local1(out_im, filename, dimx, dimy,dimz,0,0,NULL,NULL);
    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Return OK:
    return P3D_ERROR;
    
}

