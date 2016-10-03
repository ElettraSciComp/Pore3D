#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <stdio.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dBoundingBoxList.h"
#include "Common/p3dConnectedComponentsLabeling.h"
#include "Common/p3dSquaredEuclideanDT.h"
#include "Common/p3dUtils.h"
#include "Common/p3dThinning.h"

int __p3dTmpWriteRaw8(
        unsigned char* in_im,
        char* filename,
        const int dimx,
        const int dimy,
        const int dimz
        ) {
    FILE* fvol;

    /* Get a handler for the input file */
    if ((fvol = fopen(filename, "wb")) == NULL) {
        printf("Cannot open output file %s.", filename);
        in_im = NULL;

        return P3D_MEM_ERROR;
    }

    /* Write raw data to file: */
    fwrite(in_im, sizeof (unsigned char), dimx * dimy*dimz, fvol);

    /* Close file handler: */
    fclose(fvol);

    return P3D_SUCCESS;
}

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

int __p3dTmpWriteRaw16(
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

        return P3D_ERROR;
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

            return P3D_ERROR;
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

                return P3D_ERROR;
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


int _findNode(
        unsigned short* im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int i,
        const int j,
        const int k
        ) {
    int a, b, c;
    int ct = 0;

    // 6-connection:
    c = k - 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    // Do other 12 tests for 18-connection

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];


    // On k plane:
    c = k;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    // Do other 8 tests for 26-connectivity

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k - 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    c = k + 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        if (im[ I(a, b, c, dimx, dimy) ] != BACKGROUND)
            return im[ I(a, b, c, dimx, dimy) ];

    return BACKGROUND;
}

int _countNeighbors(
        unsigned short* im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int i,
        const int j,
        const int k,
        unsigned short lbl
        ) {
    int a, b, c;
    int ct = 0;

    // 6-connection:
    c = k - 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    // Do other 12 tests for 18-connection

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;


    // On k plane:
    c = k;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    // Do other 8 tests for 26-connectivity

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k - 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    c = k + 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx))
        ct += (im[ I(a, b, c, dimx, dimy) ] == lbl) ? 1 : 0;

    // Return number of voxels in the neighborhood:
    return ct;
}

int _findNeighbor(
        unsigned short* im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int i,
        const int j,
        const int k,
        coords_t* coords
        ) {
    int a, b, c;
    int ct = 0;

    // 6-connection:
    c = k - 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    // Do other 12 tests for 18-connection

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j + 1;
    a = i;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }


    // On k plane:
    c = k;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    // Do other 8 tests for 26-connectivity

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k - 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j - 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j + 1;
    a = i - 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    c = k + 1;
    b = j + 1;
    a = i + 1;
    if ((c >= 0) && (c < dimz) && (b >= 0) && (b < dimy) && (a >= 0) && (a < dimx)) {
        ct += (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) ? 1 : 0;
        if (im[ I(a, b, c, dimx, dimy) ] == USHRT_MAX) {
            coords->x = a;
            coords->y = b;
            coords->z = c;
        }
    }

    // Return number of voxels in the neighborhood:
    return ct;
}

int _p3dSkeletonAnalysis_EndPoints(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* ends_im, // OUT: Image with maximal balls filled on endpoints
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize // IN: voxel resolution
        ) {
    unsigned short* tmp_im = NULL;

    double max_width;

    int cc_array_numel;
    unsigned int* cc_array = NULL;

    bb_t* bbs = NULL;
    bb_t curr_bb;
    int ct, i, j, k, a, b, c;
    int rad;
    double delta;

    // Allocate memory:
    P3D_TRY(tmp_im = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));

    // Fill the balls:
#pragma omp parallel for private(i, j, rad, a, b, c, delta)
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                if (lbl_skl_im[ I(i, j, k, dimx, dimy) ] == END_LABEL) {
                    // Get the radius:
                    rad = (int) (sqrt((double) dt_im [ I(i, j, k, dimx, dimy) ]));

                    // Fill the ball:
                    for (c = k - rad; c <= k + rad; c++)
                        for (b = j - rad; b <= j + rad; b++)
                            for (a = i - rad; a <= i + rad; a++) {
                                // We are scanning the bounding box, so we need to be sure
                                // if current position (a,b,c) is inside the ball:
                                delta = sqrt((double) ((a - i)*(a - i) + (b - j)*(b - j)
                                        + (c - k)*(c - k)));

                                if ((a >= 0) && (b >= 0) && (c >= 0) &&
                                        (a < dimx) && (b < dimy) && (c < dimz) &&
                                        (delta < rad)) {
                                    ends_im [ I(a, b, c, dimx, dimy) ] = OBJECT;
                                }
                            }
                }
            }


    // Perform connected components labeling of END points:
    P3D_TRY(p3dConnectedComponentsLabeling(ends_im, tmp_im, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN6, P3D_TRUE));



    // If there are endpoints:
    if (cc_array_numel != 0) {
        // Allocate memory for the distribution of widths on endpoints:
        P3D_TRY(out_stats->End_Width = (double*) malloc(cc_array_numel * sizeof (double)));
        out_stats->End_Counter = cc_array_numel;

        // Compute max width on endpoints:
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            max_width = 0.0;

            // Scan the bounding box for maximum value:			
            for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
                for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                    for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                        if (ends_im[ I(i, j, k, dimx, dimy) ] == OBJECT) {
                            max_width = MAX(max_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                        }
                    }


            // Copy statistics to output structure:
            out_stats->End_Width[ct] = 2 * sqrt(max_width) * voxelsize;
        }
    } else {
        // Copy empty statistics to output structure:
        out_stats->End_Counter = 0;
        out_stats->End_Width = NULL;
    }



    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);

    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);

    return P3D_MEM_ERROR;
}

int _p3dSkeletonAnalysis_NodePoints(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* nodes_im, // OUT: Image with maximal balls filled on nodepoints
        unsigned char* pores_im, // OUT: Image with only the maximal balls for pores
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const double merging_factor,
        const double voxelsize // IN: voxel resolution
        ) {
    unsigned short* tmp_im = NULL;
    unsigned short* tmp_roi = NULL;

    double max_width;

    int cc_array_numel;
    int max_i, max_j, max_k;
    unsigned int* cc_array = NULL;

    bb_t* bbs = NULL;
    bb_t curr_bb;
    int i, j, k, a, b, c;
    int ct, rad, delta;
    coords_t coords;

    int offset = 1; // Extra-marging for bounding box


    // Allocate memory:	
    P3D_TRY(tmp_im = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));
    P3D_TRY(tmp_roi = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));

    // A ball is filled on each node voxel. The radius of this ball is the value of 
    // the distance transform on the node voxel. While, in principle, the medialness
    // property should guarantee that the node voxel lies on the center of the maximal
    // inscribed sphere this is actually not true. Therefore, assuming the value of 
    // the distance transform in that voxel as pore thickness might lead to underestimated
    // measurements. This will be fixed by the merging criterion and further processing.


    // Scan the image and look for each node:
#pragma omp parallel for private(i, j, rad, a, b, c, delta)
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                if (lbl_skl_im[ I(i, j, k, dimx, dimy) ] == NODE_LABEL) {
                    // Now that we are on a node, get the radius. If a very small
                    // merging factor has been specified, rad could be too small.
                    // A radius lower than 1 is corrected to 1.
                    rad = (int) (merging_factor * sqrt((double) dt_im [ I(i, j, k, dimx, dimy) ]));
                    if (rad < 1) rad = 1;

                    // Scan the "squared" bounding box:
                    for (c = k - rad; c <= k + rad; c++)
                        for (b = j - rad; b <= j + rad; b++)
                            for (a = i - rad; a <= i + rad; a++) {
                                // We are scanning the bounding box, so we need to be sure
                                // if current position (a,b,c) is inside the ball:
                                delta = (a - i)*(a - i) + (b - j)*(b - j) + (c - k)*(c - k);

                                if ((a >= 0) && (b >= 0) && (c >= 0) &&
                                        (a < dimx) && (b < dimy) && (c < dimz) &&
                                        (delta <= rad * rad)) {
                                    // Fill the ball:
                                    nodes_im [ I(a, b, c, dimx, dimy) ] = OBJECT;
                                }
                            }
                }
            }

    // Now nodes_im contains the cluster of balls and some of these balls overlap as 
    // more than one node might occur within a pore. In order to correctly assess the
    // number of pores, we assume a 1-1 correspondence of a cluster (set of overlapped
    // balls) with a pore. We now count the number of clusters:

    P3D_TRY(p3dConnectedComponentsLabeling(nodes_im, tmp_im, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN6, P3D_TRUE));

    // Now that we know the number of pores we can allocate memory:
    if (cc_array_numel != 0) {

        // Start stuff for connectivity density:
        out_stats->ConnectivityDensity = (double) (cc_array_numel);

        // Allocate memory for the pore-thickness distribution:
        P3D_TRY(out_stats->Node_Width = (double*) malloc(cc_array_numel * sizeof (double)));
        out_stats->Node_Counter = cc_array_numel;

        // Allocate memory for the coordination number distribution:
        P3D_TRY(out_stats->CoordinationNumber = (int*) calloc(cc_array_numel, sizeof (int)));

        // Compute pore thickness distribution for each pore (i.e. cluster of balls):
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            max_width = 0.0;

            // Reset the copy of the ROI of the bounding box:
            memset(tmp_roi, 0, dimx * dimy * dimz * sizeof (unsigned short));

            // Scan the bounding box of the cluster of balls:			
            for (k = (curr_bb.min_z - offset); k <= (curr_bb.max_z + offset); k++)
                for (j = (curr_bb.min_y - offset); j <= (curr_bb.max_y + offset); j++)
                    for (i = (curr_bb.min_x - offset); i <= (curr_bb.max_x + offset); i++) {
                        if ((i >= 0) && (j >= 0) && (k >= 0) &&
                                (i < dimx) && (j < dimy) && (k < dimz)) {

                            if (tmp_im[ I(i, j, k, dimx, dimy) ] == ((unsigned short) (ct + 3))) {
                                // The maximum value of the distance transform is assumed as pore
                                // thickness. This value is not necessarily the value of the distance
                                // transform on one of the skeleton nodes that have originated the 
                                // cluster. We should record also the position of this maximum.
                                if (((double) dt_im [ I(i, j, k, dimx, dimy) ]) > max_width) {
                                    max_width = (double) dt_im [ I(i, j, k, dimx, dimy) ];
                                    max_i = i;
                                    max_j = j;
                                    max_k = k;
                                }
                            }

                            // The following code is for further coordination number assessment:

                            // Create a temporary copy of the ROI (initialized with the same dimension of the 
                            // whole image for simplicity even if it's a waste of memory) of the bounding box.
                            // At this point, tmp_im is unsigned short with labeled nodes and lbl_skl_im is 
                            // unsigned short with classification of branches. So the ROI is created with the 
                            // labeled current node (direct copy from tmp_im) and NODE-TO-NODE and NODE-TO-END 
                            // assigned to USHRT_MAX.

                            if ((lbl_skl_im[ I(i, j, k, dimx, dimy) ] == NODETONODE_LABEL) ||
                                    (lbl_skl_im[ I(i, j, k, dimx, dimy) ] == NODETOEND_LABEL)) {
                                tmp_roi[ I(i, j, k, dimx, dimy) ] = USHRT_MAX;
                            }


                            if (tmp_im[ I(i, j, k, dimx, dimy) ] != BACKGROUND)
                                tmp_roi[ I(i, j, k, dimx, dimy) ] = tmp_im[ I(i, j, k, dimx, dimy) ];

                        }
                    }

            // Compute the coordination number re-scanning bounding box:
            for (k = (curr_bb.min_z - offset); k <= (curr_bb.max_z + offset); k++)
                for (j = (curr_bb.min_y - offset); j <= (curr_bb.max_y + offset); j++)
                    for (i = (curr_bb.min_x - offset); i <= (curr_bb.max_x + offset); i++) {
                        if ((i >= 0) && (j >= 0) && (k >= 0) &&
                                (i < dimx) && (j < dimy) && (k < dimz)) {
                            // If a branch is found (USHRT_MAX voxel):
                            if (tmp_roi[ I(i, j, k, dimx, dimy) ] == USHRT_MAX) {
                                // Increment coordination number if the current node label is found in
                                // the neighborhood:
                                if (_countNeighbors(tmp_roi, dimx, dimy, dimz, i, j, k,
                                        (unsigned short) (ct + 3)) >= 1) {
                                    out_stats->CoordinationNumber[ct]++;

                                    // Remove current branch in order to avoid counting more than once
                                    // non-perfect interconnections with the maximal ball (simple points 
                                    // occur):							
                                    tmp_roi[ I(i, j, k, dimx, dimy) ] = BACKGROUND;

                                    while (_findNeighbor(tmp_roi, dimx, dimy, dimz, i, j, k, &coords) != 0) {
                                        a = i;
                                        b = j;
                                        c = k;
                                        while (_findNeighbor(tmp_roi, dimx, dimy, dimz, a, b, c, &coords) != 0) {
                                            a = coords.x;
                                            b = coords.y;
                                            c = coords.z;
                                            tmp_roi[ I(a, b, c, dimx, dimy) ] = BACKGROUND;
                                        }
                                    }
                                }
                            }
                        }
                    }


            // Create the pore image:            
            rad = (int) sqrt(max_width);

            // Scan the "squared" bounding box:
            for (c = max_k - rad; c <= max_k + rad; c++)
                for (b = max_j - rad; b <= max_j + rad; b++)
                    for (a = max_i - rad; a <= max_i + rad; a++) {
                        // We are scanning the bounding box, so we need to be sure
                        // if current position (a,b,c) is inside the ball:
                        delta = (a - max_i)*(a - max_i) + (b - max_j)*(b - max_j)
                                + (c - max_k)*(c - max_k);

                        if ((a >= 0) && (b >= 0) && (c >= 0) &&
                                (a < dimx) && (b < dimy) && (c < dimz) &&
                                (delta <= rad * rad)) {
                            // Fill the ball:
                            pores_im [ I(a, b, c, dimx, dimy) ] = OBJECT;
                        }
                    }


            // Copy statistics to output structure:
            out_stats->Node_Width[ct] = 2 * sqrt(max_width) * voxelsize;
        }
    } else {
        // Copy empty statistics to output structure:
        out_stats->Node_Counter = 0;
        out_stats->Node_Width = NULL;
        out_stats->ConnectivityDensity = 0.0;
    }



    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_roi != NULL) free(tmp_roi);

    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_roi != NULL) free(tmp_roi);


    return P3D_MEM_ERROR;
}

int _p3dSkeletonAnalysis_NodeToNodeBranches(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* nodes_im, // IN: Input image of identified nodes with merging criteria
        struct SkeletonStats* out_stats, // IN/OUT: Skeleton statistics
        unsigned char* throats_im, // OUT: Input image of identified throats with merging criteria
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize // IN: voxel resolution
        ) {
    unsigned char* tmp_im = NULL;
    unsigned short* tmp_im2 = NULL;

    double min_width, mean_width, max_width, delta;

    unsigned int* cc_array = NULL;
    int cc_array_numel;

    bb_t* bbs = NULL;
    bb_t curr_bb;

    int ct, i, j, k, a, b, c, rad;
    int min_coord_x, min_coord_y, min_coord_z;


    //
    // Isolate NODE-TO-NODE branches:
    //

    // Allocate memory for temp image skeleton:
    P3D_TRY(tmp_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    P3D_TRY(tmp_im2 = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));



    // Create temporary matrix removing the filled nodes. Doing so, only the part of a 
    // branch outside nodes_im, i.e. the image with filled balls on skeleton nodes, is 
    // taken into acccount.
#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if (lbl_skl_im[ ct ] == NODETONODE_LABEL) {
            tmp_im [ ct ] = OBJECT;
        }
        if (nodes_im[ ct ] == OBJECT) {
            tmp_im [ ct ] = BACKGROUND;
        }
    }


    // Perform connected components labeling of NODE-TO-NODE branches:
    P3D_TRY(p3dConnectedComponentsLabeling(tmp_im, tmp_im2, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN26, P3D_TRUE));

    // Tmp_im2 is no needed anymore:
    if (tmp_im2 != NULL) free(tmp_im2);
    tmp_im2 = NULL;

    if (cc_array_numel != 0) {
        // Allocate memory for the distribution of widths on endpoints:
        P3D_TRY(out_stats->NodeToNode_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToNode_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToNode_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToNode_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

        // Copy NODE-TO-END branches statistics to output structure:
        out_stats->NodeToNode_Counter = cc_array_numel;

        // At this point connectivity density contains the number of node points, so
        // it suffices to subtract the number of NODE-TO-NODE branches:
        out_stats->ConnectivityDensity = out_stats->ConnectivityDensity - (double) (cc_array_numel);

        // Compute mean length on NODE-TO-NODE branches:
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            min_width = LONG_MAX;
            mean_width = 0.0;
            max_width = LONG_MIN;

            // Scan the bounding box for minimum, mean and maximum value:			
            for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
                for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                    for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                        if (tmp_im[ I(i, j, k, dimx, dimy) ] == OBJECT) {
                            max_width = MAX(max_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                            mean_width += (double) dt_im [ I(i, j, k, dimx, dimy) ];

                            if ((double) (dt_im [ I(i, j, k, dimx, dimy) ]) < min_width) {
                                min_width = (double) (dt_im [ I(i, j, k, dimx, dimy) ]);
                                min_coord_x = i;
                                min_coord_y = j;
                                min_coord_z = k;
                            }
                        }
                    }

            // Get the radius:
            rad = (int) (sqrt(min_width));

            // Fill the ball on throat:		
            for (c = min_coord_z - rad; c <= min_coord_z + rad; c++)
                for (b = min_coord_y - rad; b <= min_coord_y + rad; b++)
                    for (a = min_coord_x - rad; a <= min_coord_x + rad; a++) {
                        //skip_throat = FALSE;

                        // We are scanning the bounding box, so we need to be sure
                        // if current position (a,b,c) is inside the ball:
                        delta = (double) ((a - min_coord_x)*(a - min_coord_x) +
                                (b - min_coord_y)*(b - min_coord_y) + (c - min_coord_z)*(c - min_coord_z));

                        if ((a >= 0) && (b >= 0) && (c >= 0) &&
                                (a < dimx) && (b < dimy) && (c < dimz) &&
                                (delta < min_width)) {
                            throats_im [ I(a, b, c, dimx, dimy) ] = OBJECT;
                            /*if ( ( a == 0 ) || ( b == 0 ) || ( c == 0 ) ||
                            ( a == (dimx-1) ) || ( b == (dimy-1) ) || ( c == (dimz-1) ) )
                            {
                            skip_throat = TRUE;
                            }*/
                        }
                    }

            mean_width /= cc_array[ ct ];

            // Copy NODE-TO-NODE branches statistics to output structure:
            out_stats->NodeToNode_Length[ct] = cc_array[ ct ] * voxelsize;
            out_stats->NodeToNode_MinWidth[ct] = 2 * sqrt(min_width) * voxelsize;
            out_stats->NodeToNode_MeanWidth[ct] = 2 * sqrt(mean_width) * voxelsize;
            out_stats->NodeToNode_MaxWidth[ct] = 2 * sqrt(max_width) * voxelsize;
        }
    } else {
        // Copy NODE-TO-NODE branches statistics to output structure:
        out_stats->NodeToNode_Counter = 0;
        out_stats->NodeToNode_Length = NULL;
        out_stats->NodeToNode_MinWidth = NULL;
        out_stats->NodeToNode_MeanWidth = NULL;
        out_stats->NodeToNode_MaxWidth = NULL;
    }

    // TODO: Reshape arrays due to the skip_throat flag

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);

    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    return P3D_MEM_ERROR;
}

int _p3dSkeletonAnalysis_NodeToEndBranches(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* nodes_im, // IN: Input image of identified nodes with merging criteria
        unsigned char* ends_im, // IN: Input image of identified nodes with merging criteria
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize // IN: voxel resolution
        ) {
    unsigned char* tmp_im = NULL;
    unsigned short* tmp_im2 = NULL;

    double min_width, mean_width, max_width;

    unsigned int* cc_array = NULL;
    int cc_array_numel;

    bb_t* bbs = NULL;
    bb_t curr_bb;

    int ct, i, j, k;


    //
    // Isolate NODE-TO-NODE branches:
    //

    // Allocate memory for temp image skeleton:
    P3D_TRY(tmp_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    P3D_TRY(tmp_im2 = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));


    // Create temporary matrix removing the filled balls. Doing so, only the part of a 
    // branch outside nodes_im, i.e. the image with filled balls on skeleton nodes, and 
    // outside ends_im, i.e. the image with filled balls on skeleton ends, is taken into 
    // acccount.
#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if (lbl_skl_im[ ct ] == NODETOEND_LABEL) {
            tmp_im [ ct ] = OBJECT;
        }
        if (nodes_im[ ct ] == OBJECT) {
            tmp_im [ ct ] = BACKGROUND;
        }
        if (ends_im[ ct ] == OBJECT) {
            tmp_im [ ct ] = BACKGROUND;
        }
    }


    // Perform connected components labeling of NODE-TO-END branches:
    P3D_TRY(p3dConnectedComponentsLabeling(tmp_im, tmp_im2, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN26, P3D_TRUE));


    if (cc_array_numel != 0) {
        // Allocate memory for the distribution of widths on endpoints:
        P3D_TRY(out_stats->NodeToEnd_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToEnd_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToEnd_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->NodeToEnd_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

        // Copy NODE-TO-END branches statistics to output structure:
        out_stats->NodeToEnd_Counter = cc_array_numel;

        // Uncomment this line for taking into account NODE-TO-END branches
        // in the computation of the connectivity density. However, taking into 
        // account only NODE-TO-NODE branches the connectivity density can be 
        // correctly interpreted as the number of redundant connections per 
        // volume element.
        //out_stats->ConnectivityDensity = out_stats->ConnectivityDensity - cc_array_numel;


        // Compute mean length on NODE-TO-NODE branches:
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            min_width = LONG_MAX;
            mean_width = 0.0;
            max_width = LONG_MIN;

            // Scan the bounding box for minimum, mean and maximum value:			
            for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
                for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                    for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                        if (tmp_im[ I(i, j, k, dimx, dimy) ] == OBJECT) {
                            min_width = MIN(min_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                            max_width = MAX(max_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                            mean_width += (double) dt_im [ I(i, j, k, dimx, dimy) ];
                        }
                    }

            mean_width /= cc_array[ ct ];

            // Copy NODE-TO-NODE branches statistics to output structure:
            out_stats->NodeToEnd_Length[ct] = cc_array[ ct ] * voxelsize;
            out_stats->NodeToEnd_MinWidth[ct] = 2 * sqrt(min_width) * voxelsize;
            out_stats->NodeToEnd_MeanWidth[ct] = 2 * sqrt(mean_width) * voxelsize;
            out_stats->NodeToEnd_MaxWidth[ct] = 2 * sqrt(max_width) * voxelsize;
        }
    } else {
        // Copy NODE-TO-NODE branches statistics to output structure:
        out_stats->NodeToEnd_Counter = 0;
        out_stats->NodeToEnd_Length = NULL;
        out_stats->NodeToEnd_MinWidth = NULL;
        out_stats->NodeToEnd_MeanWidth = NULL;
        out_stats->NodeToEnd_MaxWidth = NULL;
    }

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);


    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    return P3D_MEM_ERROR;
}

int _p3dSkeletonAnalysis_EndToEndBranches(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* ends_im, // IN: Input image of identified nodes with merging criteria
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize // IN: voxel resolution
        ) {
    unsigned char* tmp_im = NULL;
    unsigned short* tmp_im2 = NULL;

    double min_width, mean_width, max_width;

    unsigned int* cc_array = NULL;
    int cc_array_numel;

    bb_t* bbs = NULL;
    bb_t curr_bb;

    int ct, i, j, k;




    //
    // Isolate END-TO-END branches:
    //

    // Allocate memory for labeled skeleton:
    P3D_TRY(tmp_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));
    P3D_TRY(tmp_im2 = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));


    // Set memory of tmp_im:
    memset(tmp_im, BACKGROUND, dimx * dimy * dimz * sizeof (unsigned char));


    // Create temporary matrix removing the filled balls. Doing so, only the part of a 
    // branch outside ends_im, i.e. the image with filled balls on skeleton ends, is 
    // taken into acccount.
#pragma omp parallel for
    for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if (lbl_skl_im[ ct ] == ENDTOEND_LABEL) {
            // Assign temporary image:
            tmp_im [ ct ] = OBJECT;
        }
        if (ends_im[ ct ] == OBJECT) {
            // Assign temporary image:
            tmp_im [ ct ] = BACKGROUND;
        }
    }


    // Perform connected components labeling of NODE-TO-NODE branches:
    P3D_TRY(p3dConnectedComponentsLabeling(tmp_im, tmp_im2, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN26, P3D_TRUE));


    if (cc_array_numel != 0) {
        // Allocate memory for the distribution of widths on endpoints:
        P3D_TRY(out_stats->EndToEnd_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->EndToEnd_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->EndToEnd_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_TRY(out_stats->EndToEnd_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

        // Copy NODE-TO-END branches statistics to output structure:
        out_stats->EndToEnd_Counter = cc_array_numel;

        // Compute mean length on NODE-TO-NODE branches:
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            min_width = LONG_MAX;
            mean_width = 0.0;
            max_width = LONG_MIN;

            // Scan the bounding box for minimum, mean and maximum value:			
            for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
                for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                    for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                        if (tmp_im[ I(i, j, k, dimx, dimy) ] == OBJECT) {
                            min_width = MIN(min_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                            max_width = MAX(max_width, (double) dt_im [ I(i, j, k, dimx, dimy) ]);
                            mean_width += (double) dt_im [ I(i, j, k, dimx, dimy) ];
                        }
                    }

            mean_width /= cc_array[ ct ];

            // Copy NODE-TO-NODE branches statistics to output structure:
            out_stats->EndToEnd_Length[ct] = cc_array[ ct ] * voxelsize;
            out_stats->EndToEnd_MinWidth[ct] = 2 * sqrt(min_width) * voxelsize;
            out_stats->EndToEnd_MeanWidth[ct] = 2 * sqrt(mean_width) * voxelsize;
            out_stats->EndToEnd_MaxWidth[ct] = 2 * sqrt(max_width) * voxelsize;
        }
    } else {
        // Copy NODE-TO-NODE branches statistics to output structure:
        out_stats->EndToEnd_Counter = 0;
        out_stats->EndToEnd_Length = NULL;
        out_stats->EndToEnd_MinWidth = NULL;
        out_stats->EndToEnd_MeanWidth = NULL;
        out_stats->EndToEnd_MaxWidth = NULL;
    }

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);


    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources after each call to _p3dConnectedComponentsLabeling:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    return P3D_MEM_ERROR;
}


int p3dSkeletonAnalysis(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned char* skl_im, // IN: Input (binary) skeleton of the segmented volume
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        unsigned char* nodes_im,
        unsigned char* pores_im,
        unsigned char* ends_im,
        unsigned char* throats_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const double merging_factor,
        const double voxelsize, // IN: voxel resolution
        int (*wr_log)(const char*, ...)
        ) {

    // Temporary matrices:
    unsigned char* max_skl_im;
    unsigned char* lbl_skl_im;
    unsigned short* dt_im;

    double mean = 0.0;
    double mean_sqr = 0.0;
    double glob_mean = 0.0;
    double glob_mean_sqr = 0.0;
    int ct, i, j;
    int flag_nodes_null = P3D_FALSE;
    int flag_pores_null = P3D_FALSE;
    int flag_ends_null = P3D_FALSE;
    int flag_throats_null = P3D_FALSE;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing skeleton analysis...");
        wr_log("\tVoxelsize: %0.6f mm.", voxelsize);
        wr_log("\tMerging factor: %0.2f.", merging_factor);
    }


    // Allocate memory:
    if (nodes_im == NULL) {
        flag_nodes_null = P3D_TRUE;
        P3D_TRY(nodes_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (pores_im == NULL) {
        flag_pores_null = P3D_TRUE;
        P3D_TRY(pores_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (ends_im == NULL) {
        flag_ends_null = P3D_TRUE;
        P3D_TRY(ends_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (throats_im == NULL) {
        flag_throats_null = P3D_TRUE;
        P3D_TRY(throats_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }

    P3D_TRY(dt_im = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));
    P3D_TRY(max_skl_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));
    P3D_TRY(lbl_skl_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));


    // Compute distance transform for further use:
    P3D_TRY(p3dSquaredEuclideanDT(vol_im, dt_im, dimx, dimy, dimz));


    // Ensure that only one skeleton image is presented (the maximum):
    P3D_TRY(p3dGetMaxVolumeRegion(skl_im, max_skl_im, dimx, dimy, dimz, CONN26));
    P3D_TRY(p3dThinning(max_skl_im, dimx, dimy, dimz));


    // Perform skeleton labeling:
    P3D_TRY(p3dSkeletonLabeling(max_skl_im, lbl_skl_im, dimx, dimy, dimz, NULL));


    // ANALYSIS STEP. Computes number and lengths of each identified component.
    // This implementation uses connected component labeling for segments but 
    // also for nodes to mitigate thinning algorithm defects:


    P3D_TRY(_p3dSkeletonAnalysis_EndPoints(vol_im, dt_im, lbl_skl_im, ends_im, out_stats,
            dimx, dimy, dimz, voxelsize));
    P3D_TRY(_p3dSkeletonAnalysis_NodePoints(vol_im, dt_im, lbl_skl_im, nodes_im, pores_im, out_stats,
            dimx, dimy, dimz, merging_factor, voxelsize));
    P3D_TRY(_p3dSkeletonAnalysis_NodeToNodeBranches(vol_im, dt_im, lbl_skl_im, nodes_im,
            out_stats, throats_im, dimx, dimy, dimz, voxelsize));
    P3D_TRY(_p3dSkeletonAnalysis_NodeToEndBranches(vol_im, dt_im, lbl_skl_im, nodes_im, ends_im,
            out_stats, dimx, dimy, dimz, voxelsize));
    P3D_TRY(_p3dSkeletonAnalysis_EndToEndBranches(vol_im, dt_im, lbl_skl_im, ends_im,
            out_stats, dimx, dimy, dimz, voxelsize));
    //P3D_TRY(_p3dSkeletonAnalysis_Tortousity(vol_im, dt_im, lbl_skl_im, nodes_im, ends_im,
    //        out_stats, dimx, dimy, dimz, tortuosity_depth));

    // Set connectivity index. At this point, out_stats->ConnectivityDensity contains the 
    // Euler number, i.e. nodes - branches (only NODE-TO-NODEs are considered). We need to 
    // compute the final connectivity index as (1 - Euler#)/Volume:
    out_stats->ConnectivityDensity = (1.0 - out_stats->ConnectivityDensity) /
            (dimx * voxelsize * dimy * voxelsize * dimz * voxelsize);


    // Print out number of connected components and mean values of parameters:
    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tNumber of PORES: %d. ", out_stats->Node_Counter);
        wr_log("\tNumber of ENDS: %d.", out_stats->End_Counter);
        wr_log("\tNumber of NODE-TO-NODE branches (throats): %d.", out_stats->NodeToNode_Counter);
        wr_log("\tNumber of NODE-TO-END branches: %d.", out_stats->NodeToEnd_Counter);
        wr_log("\tNumber of END-TO-END branches: %d.", out_stats->EndToEnd_Counter);
        wr_log("\t----");


        if (out_stats->Node_Counter > 0) {
            // Compute mean values of volume:
            mean = 0.0;
            mean_sqr = 0.0;
            for (ct = 0; ct < (out_stats->Node_Counter); ct++) {
                mean_sqr += (double) (out_stats->Node_Width[ct] * out_stats->Node_Width[ct]);
                mean += (double) (out_stats->Node_Width[ct]);
            }
            mean = mean / ((double) (out_stats->Node_Counter));
            mean_sqr = mean_sqr / ((double) (out_stats->Node_Counter));
            mean_sqr = mean_sqr - mean*mean;

            wr_log("\tPore size: %0.3f +/- %0.3f [mm].", mean, sqrt(mean_sqr));
        }

        if (out_stats->NodeToNode_Counter > 0) {
            // Compute mean values of volume:
            mean = 0.0;
            mean_sqr = 0.0;
            for (ct = 0; ct < out_stats->NodeToNode_Counter; ct++) {
                mean_sqr += (double) (out_stats->NodeToNode_MinWidth[ct] * out_stats->NodeToNode_MinWidth[ct]);
                mean = mean + (double) (out_stats->NodeToNode_MinWidth[ct]);
                
                //wr_log("\tThroat size: %0.3f +/- %0.3f [mm].", (double) (out_stats->NodeToNode_MinWidth[ct]));
            }
            mean = mean / ((double) (out_stats->NodeToNode_Counter));
            mean_sqr = mean_sqr / ((double) (out_stats->NodeToNode_Counter));
            mean_sqr = mean_sqr - mean*mean;

            wr_log("\tThroat size: %0.3f +/- %0.3f [mm].", mean, sqrt(mean_sqr));
        }

        if (out_stats->Node_Counter > 0) {

            wr_log("\tConnectivity density: %0.3f [mm^-3].", out_stats->ConnectivityDensity);
        }

        if (out_stats->Node_Counter > 0) {
            // Compute mean values of volume:
            mean = 0.0;
            mean_sqr = 0.0;
            for (ct = 0; ct < (out_stats->Node_Counter); ct++) {
                mean = mean + (double) (out_stats->CoordinationNumber[ct]);
                mean_sqr += (double) (out_stats->CoordinationNumber[ct] * out_stats->CoordinationNumber[ct]*1.0);
                //wr_log ("\t\tCoordination number: %d ", out_stats->CoordinationNumber[ct] );
            }
            mean = mean / ((double) (out_stats->Node_Counter));
            mean_sqr = mean_sqr / ((double) (out_stats->Node_Counter));
            mean_sqr = mean_sqr - mean*mean;

            wr_log("\tCoordination number: %0.3f +/- %0.3f [-].", mean, sqrt(mean_sqr));
        }       
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Skeleton analysis performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }


    // Release resources:
    if (dt_im != NULL) free(dt_im);
    if (max_skl_im != NULL) free(max_skl_im);
    if (lbl_skl_im != NULL) free(lbl_skl_im);
    if (flag_nodes_null == P3D_TRUE) {
        if (nodes_im != NULL) free(nodes_im);
        nodes_im = NULL;
    }
    if (flag_ends_null == P3D_TRUE) {
        if (ends_im != NULL) free(ends_im);
        ends_im = NULL;
    }
    if (flag_throats_null == P3D_TRUE) {
        if (throats_im != NULL) free(throats_im);
        throats_im = NULL;
    }

    // Return OK:
    return P3D_SUCCESS;


MEM_ERROR:

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Free memory if previous malloc were successfully:
    if (dt_im != NULL) free(dt_im);
    if (max_skl_im != NULL) free(max_skl_im);
    if (lbl_skl_im != NULL) free(lbl_skl_im);
    if (flag_nodes_null == P3D_TRUE) {
        if (nodes_im != NULL) free(nodes_im);
        nodes_im = NULL;
    }
    if (flag_ends_null == P3D_TRUE) {
        if (ends_im != NULL) free(ends_im);
        ends_im = NULL;
    }
    if (flag_throats_null == P3D_TRUE) {
        if (throats_im != NULL) free(throats_im);
        throats_im = NULL;
    }


    // Return error code and exit:
    return P3D_MEM_ERROR;
}







