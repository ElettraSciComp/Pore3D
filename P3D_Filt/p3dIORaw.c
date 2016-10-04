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
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"

short _EndianSwapSignedShort(short s) {
    unsigned char b1, b2;

    b1 = s & 255;
    b2 = (s >> 8) & 255;

    return (b1 << 8) +b2;
}

unsigned short _EndianSwapUnsignedShort(unsigned short s) {
    unsigned char b1, b2;

    b1 = s & 255;
    b2 = (s >> 8) & 255;

    return (b1 << 8) +b2;
}

unsigned int _EndianSwapUnsignedInt( unsigned int val )
{
    val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF ); 
    return (val << 16) | (val >> 16);
}

int _EndianSwapInt( int val )
{
    val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
    return (val << 16) | ((val >> 16) & 0xFFFF);
}


int p3dReadRaw8(
        char* filename,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    FILE* fvol;

    /*char* auth_code;*/

    //
    // Authenticate:
    //
    /* auth_code = authenticate("p3dReadRaw8");
     if (strcmp(auth_code, "1") != 0) goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Reading RAW file %s ...", filename);        
    }

    /* Get a handler for the input file */
    if ((fvol = fopen(filename, "rb")) == NULL) {
        printf("Cannot open input file %s.", filename);
        out_im = NULL;

        return P3D_IO_ERROR;
    }

    /* Read raw data from file: */
    if (fread(out_im, sizeof (unsigned char), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
        wr_log("Pore3D - IO error: error on reading file %s. Program will exit.", filename);

        return P3D_IO_ERROR;
    }

    /* Close file handler: */
    fclose(fvol);

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - RAW file read successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec()) ;
    }


    // Return OK:
    return P3D_SUCCESS;


    /*AUTH_ERROR:

        if (wr_log != NULL) {
            wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
        }

        return P3D_AUTH_ERROR;*/
}

int p3dWriteRaw8(
        unsigned char* in_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    FILE* fvol;

    /*char* auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dWriteRaw8");
    if (strcmp(auth_code, "1") != 0) goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Writing 8-bit RAW file %s ...", filename);
    }

    /* Get a handler for the input file */
    if ((fvol = fopen(filename, "wb")) == NULL) {
        wr_log("Pore3D - IO error: cannot open output file %s. Program will exit.", filename);

        return P3D_IO_ERROR;
    }

    /* Write raw data to file: */
    if (fwrite(in_im, sizeof (unsigned char), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
        wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

        return P3D_IO_ERROR;
    }

    /* Close file handler: */
    fclose(fvol);

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - RAW file written successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec()) ;
    }


    // Return OK:
    return P3D_SUCCESS;


    /*AUTH_ERROR:

        if (wr_log != NULL) {
            wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
        }

        return P3D_AUTH_ERROR;*/
}

int p3dReadRaw16(
        char* filename,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int flagLittle,
        const int flagSigned,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    FILE* fvol;
    short* s_tmp_im = NULL;
    unsigned short* us_tmp_im = NULL;
    int ct;

    /*char* auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dReadRaw16");
    if (strcmp(auth_code, "1") != 0) goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Reading RAW file %s ...", filename);
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
    if ((fvol = fopen(filename, "rb")) == NULL) {
        wr_log("Pore3D - IO error: cannot open input file %s. Program will exit.", filename);
        out_im = NULL;

        return P3D_IO_ERROR;
    }


    /* Read data signed/unsigned: */
    if (flagSigned == P3D_TRUE) {
        s_tmp_im = (short*) malloc(dimx * dimy * dimz * sizeof (short));

        /* Read raw data from file: */
        if (fread(s_tmp_im, sizeof (short), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
            wr_log("Pore3D - IO error: error on reading file %s. Program will exit.", filename);

            return P3D_IO_ERROR;
        }

        /* Convert to unsigned: */
        for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
            if (flagLittle == P3D_FALSE) {
                out_im[ct] = (unsigned short) (_EndianSwapSignedShort(s_tmp_im[ct]) + USHRT_MAX / 2);
            } else {
                out_im[ct] = (unsigned short) (s_tmp_im[ct] + USHRT_MAX / 2);
            }
        }

        /* Free: */
        free(s_tmp_im);
    } else {
        /* Swap endian if necessary: */
        if (flagLittle == P3D_FALSE) {
            us_tmp_im = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short));

            /* Read raw data from file: */
            if (fread(us_tmp_im, sizeof (unsigned short), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
                wr_log("Pore3D - IO error: error on reading file %s. Program will exit.", filename);

                return P3D_IO_ERROR;
            }

            for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
                out_im[ct] = (unsigned short) (_EndianSwapUnsignedShort(us_tmp_im[ct]));
            }

            /* Free: */
            free(us_tmp_im);
        } else {
            if (fread(out_im, sizeof (unsigned short), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
                wr_log("Pore3D - IO error: error on reading file %s. Program will exit.", filename);

                return P3D_IO_ERROR;
            }
        }
    }



    /* Close file handler: */
    fclose(fvol);

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - RAW file read successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }


    // Return OK:
    return P3D_SUCCESS;


    /*AUTH_ERROR:

        if (wr_log != NULL) {
            wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
        }

        return P3D_AUTH_ERROR;*/
}

int p3dWriteRaw16(
        unsigned short* in_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz,
        const int flagLittle,
        const int flagSigned,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    FILE* fvol;
    short* s_tmp_im = NULL;
    unsigned short* us_tmp_im = NULL;
    int ct;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Writing 16-bit RAW file %s ...", filename);
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

        return P3D_IO_ERROR;
    }

    /* Read data signed/unsigned: */
    if (flagSigned == P3D_TRUE) {
        s_tmp_im = (short*) malloc(dimx * dimy * dimz * sizeof (short));

        /* Convert to signed: */
        for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
            if (flagLittle == P3D_FALSE) {
                s_tmp_im[ct] = _EndianSwapSignedShort((short) (in_im[ct] - USHRT_MAX / 2));
            } else {
                s_tmp_im[ct] = (short) (in_im[ct] - USHRT_MAX / 2);
            }
        }

        /* Write raw data to file: */
        //fwrite(s_tmp_im, sizeof (short), dimx*dimy, fvol);
        if (fwrite(s_tmp_im, sizeof (short), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
            wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

            return P3D_IO_ERROR;
        }

        /* Free: */
        free(s_tmp_im);
    } else {
        /* Swap endian if necessary: */
        if (flagLittle == P3D_FALSE) {
            us_tmp_im = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short));

            for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
                us_tmp_im[ct] = _EndianSwapUnsignedShort(in_im[ct]);
            }

            //fwrite(us_tmp_im, sizeof (unsigned short), dimx*dimy, fvol);
            if (fwrite(us_tmp_im, sizeof (unsigned short), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
                wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

                return P3D_IO_ERROR;
            }

            free(us_tmp_im);
        } else {
            fwrite(in_im, sizeof (unsigned short), dimx * dimy * dimz, fvol);
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

int p3dWriteRaw32(
        unsigned int* in_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz,
        const int flagLittle,
        const int flagSigned,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    FILE* fvol;
    int* s_tmp_im = NULL;
    unsigned int* us_tmp_im = NULL;
    int ct;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Writing 32-bit RAW file %s ...", filename);
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

        return P3D_IO_ERROR;
    }

    /* Read data signed/unsigned: */
    if (flagSigned == P3D_TRUE) {
        s_tmp_im = (int*) malloc(dimx * dimy * dimz * sizeof (int));

        /* Convert to signed: */
        for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
            if (flagLittle == P3D_FALSE) {
                s_tmp_im[ct] = _EndianSwapInt((int) (in_im[ct] - UINT_MAX / 2));
            } else {
                s_tmp_im[ct] = (int) (in_im[ct] - UINT_MAX / 2);
            }
        }

        /* Write raw data to file: */
        //fwrite(s_tmp_im, sizeof (short), dimx*dimy, fvol);
        if (fwrite(s_tmp_im, sizeof (int), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
            wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

            return P3D_IO_ERROR;
        }

        /* Free: */
        free(s_tmp_im);
    } else {
        /* Swap endian if necessary: */
        if (flagLittle == P3D_FALSE) {
            us_tmp_im = (unsigned int*) malloc(dimx * dimy * dimz * sizeof (unsigned int));

            for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
                us_tmp_im[ct] = _EndianSwapUnsignedInt(in_im[ct]);
            }

            //fwrite(us_tmp_im, sizeof (unsigned short), dimx*dimy, fvol);
            if (fwrite(us_tmp_im, sizeof (unsigned int), dimx * dimy * dimz, fvol) < (dimx * dimy * dimz)) {
                wr_log("Pore3D - IO error: error on writing file %s. Program will exit.", filename);

                return P3D_IO_ERROR;
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
