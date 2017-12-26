#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdio.h>

#include "p3dSkel.h"
#include "p3dTime.h"
//#include "p3dGraph.h"

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

        return P3D_ERROR;
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

/*float _dijsktra(float* cost, int dim, int source, int target) {
    float* dist = NULL;
    int* prev = NULL;
    int* selected = NULL;
    char* path = NULL;

    float min, d, out;
    int i, m, start, ct;


    P3D_TRY(dist = (float*) calloc(dim, sizeof (float)));
    P3D_TRY(prev = (int*) calloc(dim, sizeof (int)));
    P3D_TRY(selected = (int*) calloc(dim, sizeof (int)));
    P3D_TRY(path = (char*) calloc(dim, sizeof (char)));

    for (i = 0; i < dim; i++) {
        dist[i] = (float) _DIJKSTRA_IN;
        prev[i] = (float) _DIJKSTRA_IN;
    }
    start = source;
    selected[start] = 1;
    dist[start] = 0.0;

    //while( selected[target] == 0 )
    ct = 0;
    // To avoid infinite loop in case of "imperfect network"
    while ((start != target) && (ct < (2 * dim))) {
        min = (float) _DIJKSTRA_IN;
        m = 0;
        for (i = 0; i < dim; i++) {
            if (selected[i] == 0) {
                d = dist[start] + cost[ I2(start, i, dim) ];
                if (d < dist[i]) {
                    dist[i] = d;
                    prev[i] = start;
                }
                if (min > dist[i]) {
                    min = dist[i];
                    m = i;
                }
            }
        }
        start = m;
        selected[start] = 1;
        ct++;
    }
    if (ct >= (dim * dim))
        out = -1.0;
    else
        out = dist[target];

    if (dist != NULL) free(dist);
    if (prev != NULL) free(prev);
    if (selected != NULL) free(selected);
    if (path != NULL) free(path);

    return out;

MEM_ERROR:

    if (dist != NULL) free(dist);
    if (prev != NULL) free(prev);
    if (selected != NULL) free(selected);
    if (path != NULL) free(path);

    return 0.0;
}*/

/*double _computeTortousity(float* tort_matrix, float* tort_array, int dim,
        int source, int dest) {
    float dist, real_dist;
    float a1, a2, b1, b2, c1, c2;
    //int i, j;

    // Get Euclidean distance:
    a1 = tort_array [ I2(source, 0, dim) ];
    a2 = tort_array [ I2(dest, 0, dim) ];
    b1 = tort_array [ I2(source, 1, dim) ];
    b2 = tort_array [ I2(dest, 1, dim) ];
    c1 = tort_array [ I2(source, 2, dim) ];
    c2 = tort_array [ I2(dest, 2, dim) ];
    dist = (float) sqrt((a1 - a2)*(a1 - a2) + (b1 - b2)*(b1 - b2) + (c1 - c2)*(c1 - c2));

    // Get real distance
    real_dist = p3dShortestPath(tort_matrix, dim, source, dest);
    //real_dist = _dijsktra(tort_matrix, dim, source, dest);

    /*printf("Branch [%0.1f,%0.1f,%0.1f] - [%0.1f,%0.1f,%0.1f] = %0.1f / %0.1f\n", 
    tort_array [ I2(source, 0, dim) ],
    tort_array [ I2(source, 1, dim) ],
    tort_array [ I2(source, 2, dim) ],
    tort_array [ I2(dest, 0, dim) ],
    tort_array [ I2(dest, 1, dim) ],
    tort_array [ I2(dest, 2, dim) ],
    real_dist, 
    dist
    );*/

    /*if ( fabs( real_dist + 1.0 ) < 1E-4)
        return -1.0;
    else
        return real_dist / dist;
}*/

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
    P3D_MEM_TRY(tmp_im = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));

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
        P3D_MEM_TRY(out_stats->End_Width = (double*) malloc(cc_array_numel * sizeof (double)));
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

    return P3D_ERROR;
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
    P3D_MEM_TRY(tmp_im = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));
    P3D_MEM_TRY(tmp_roi = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));

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
                    for (c = k - rad; c <= k + rad; c++) {
                        for (b = j - rad; b <= j + rad; b++) {
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
        P3D_MEM_TRY(out_stats->Node_Width = (double*) malloc(cc_array_numel * sizeof (double)));
        out_stats->Node_Counter = cc_array_numel;

        // Allocate memory for the coordination number distribution:
        P3D_MEM_TRY(out_stats->CoordinationNumber = (int*) calloc(cc_array_numel, sizeof (int)));

        // Compute pore thickness distribution for each pore (i.e. cluster of balls):
        for (ct = 0; ct < cc_array_numel; ct++) {
            curr_bb = bbs[ct];

            max_width = 0.0;

            // Reset the copy of the ROI of the bounding box:
            memset(tmp_roi, 0, dimx * dimy * dimz * sizeof (unsigned short));

            // Scan the bounding box of the cluster of balls:			
            for (k = (curr_bb.min_z - offset); k <= (curr_bb.max_z + offset); k++) {
                for (j = (curr_bb.min_y - offset); j <= (curr_bb.max_y + offset); j++) {
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
				}
			}

            // Compute the coordination number re-scanning bounding box:
            for (k = (curr_bb.min_z - offset); k <= (curr_bb.max_z + offset); k++) {
                for (j = (curr_bb.min_y - offset); j <= (curr_bb.max_y + offset); j++) {
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


    return P3D_ERROR;
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
    P3D_MEM_TRY(tmp_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    P3D_MEM_TRY(tmp_im2 = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));



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
        P3D_MEM_TRY(out_stats->NodeToNode_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToNode_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToNode_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToNode_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

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

    return P3D_ERROR;
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
    P3D_MEM_TRY(tmp_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    P3D_MEM_TRY(tmp_im2 = (unsigned short*) calloc(dimx * dimy*dimz, sizeof (unsigned short)));


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
        P3D_MEM_TRY(out_stats->NodeToEnd_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToEnd_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToEnd_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->NodeToEnd_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

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

    return P3D_ERROR;
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
    P3D_MEM_TRY(tmp_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));
    P3D_MEM_TRY(tmp_im2 = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));


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
        P3D_MEM_TRY(out_stats->EndToEnd_Length = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->EndToEnd_MinWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->EndToEnd_MeanWidth = (double*) malloc(cc_array_numel * sizeof (double)));
        P3D_MEM_TRY(out_stats->EndToEnd_MaxWidth = (double*) malloc(cc_array_numel * sizeof (double)));

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

    return P3D_ERROR;
}

int _p3dSkeletonAnalysis_Tortousity(
        unsigned char* vol_im, // IN: Input segmented (binary) volume
        unsigned short* dt_im, // IN: Input squared euclidean distance transform of vol_im
        unsigned char* lbl_skl_im, // IN: Input labeled skeleton of the segmented volume
        unsigned char* pores_im, // OUT: Image with only the maximal balls for pores
        unsigned char* ends_im, // IN: Input image of identified nodes with merging criteria
        struct SkeletonStats* out_stats, // OUT: Skeleton statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const int tortuosity_depth
        ) {

    unsigned short* tmp_im = NULL;
    unsigned short* tmp_im2 = NULL;
    unsigned char* poresplusend_im = NULL;

    float* tort_matrix = NULL;
    float* tort_array = NULL;
    int* min_ct_array = NULL;
    int* max_ct_array = NULL;
    int* coord_array = NULL;

    double max_width, dist1, dist2;

    int cc_array_numel;
    double cc_length; //, min, max;
    int max_i, max_j, max_k;
    unsigned int* cc_array = NULL;

    bb_t* bbs = NULL;
    bb_t curr_bb, curr_bb2;
    int i, j, k, a, b, c, r, s, t, ct/*, ct_num_nodes*/;
    double cen_i, cen_j, cen_k;
    int cen_ct;
    coords_t coords;

    int other_lbl;
    int offset = 1; // Extra-marging for bounding box

    int aa = 0;
    //int flag;
    int num_nodes = tortuosity_depth;


    // Allocate memory:	
    P3D_MEM_TRY(poresplusend_im = (unsigned char*) calloc(dimx * dimy * dimz, sizeof (unsigned char)));
    P3D_MEM_TRY(tmp_im = (unsigned short*) calloc(dimx * dimy * dimz, sizeof (unsigned short)));
    P3D_MEM_TRY(tmp_im2 = (unsigned short*) calloc(dimx * dimy * dimz, sizeof (unsigned short)));

    // Create an additional image with pores and also ends:
    for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if ((pores_im[ ct ] != BACKGROUND) || (ends_im[ ct ] != BACKGROUND)) {
            // Assign temporary image:
            poresplusend_im[ ct ] = UCHAR_MAX;
        }
    }
    /*for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if ((lbl_skl_im[ ct ] != NODE_LABEL) || (lbl_skl_im[ ct ] != END_LABEL)) {
            // Assign temporary image:
            poresplusend_im[ ct ] = UCHAR_MAX;
        }
    }*/

    // Get a labeled image for the pores:    
    P3D_TRY(p3dConnectedComponentsLabeling(poresplusend_im, tmp_im, &cc_array_numel, &cc_array, &bbs,
            dimx, dimy, dimz, CONN6, P3D_FALSE));
    //P3D_TRY(p3dConnectedComponentsLabeling(pores_im, tmp_im, &cc_array_numel, &cc_array, &bbs,
    //	dimx, dimy, dimz, CONN6, P3D_FALSE));

    //__p3dTmpWriteRaw16 ( tmp_im, "R:\\TEMP\\Pipp_tmp.raw", dimx, dimy, dimz, P3D_TRUE, P3D_TRUE, NULL, NULL );

    // Create a second temporary copy of the image. At this point, tmp_im is unsigned 
    // short with labeled pores and lbl_skl_im is unsigned short with classification 
    // of branches. Now, tmp_im2 is created with the labeled pores (direct copy from 
    // tmp_im) and NODE-TO-NODE and NODE-TO-END assigned to USHRT_MAX:
    for (ct = 0; ct < (dimx * dimy * dimz); ct++) {
        if ((lbl_skl_im[ ct ] == ENDTOEND_LABEL) ||
                (lbl_skl_im[ ct ] == NODETOEND_LABEL) ||
                (lbl_skl_im[ ct ] == NODETONODE_LABEL))
            //(lbl_skl_im[ ct ] == NODE_LABEL) )		            
            tmp_im2[ ct ] = USHRT_MAX;
        if (tmp_im[ ct ] != BACKGROUND)
            tmp_im2[ ct ] = tmp_im[ ct ];

    }


    // If there are pores:
    if (cc_array_numel != 0) {

        // Allocate memory for saving the center of mass of each pore:
        P3D_MEM_TRY(tort_matrix = (float*) calloc(cc_array_numel * cc_array_numel, sizeof (float)));
        P3D_MEM_TRY(tort_array = (float*) calloc(cc_array_numel * 4, sizeof (float)));
        P3D_MEM_TRY(coord_array = (int*) calloc(cc_array_numel, sizeof (int)));

        // For each pore:
        for (ct = 0; ct < cc_array_numel; ct++) {
            //for (ct = (cc_array_numel - 1); ct >= 0; ct--) {

            // Initialize the variable that will contain the value for thickness:
            max_width = 0.0;

            // Get the current bounding box:
            curr_bb = bbs[ct];

            cen_i = 0;
            cen_j = 0;
            cen_k = 0;
            cen_ct = 0;

            // Scan the bounding box of current pore:			
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
                                if (((float) (dt_im [ I(i, j, k, dimx, dimy) ])) > max_width) {
                                    max_width = (float) dt_im [ I(i, j, k, dimx, dimy) ];
                                    max_i = i;
                                    max_j = j;
                                    max_k = k;
                                }

                                // Compute also the center of mass of this pore:
                                cen_i += (float) i;
                                cen_j += (float) j;
                                cen_k += (float) k;
                                cen_ct++;
                            }
                        }
                    }

            // Baricenter:				
            cen_i = cen_i / ((float) (cen_ct));
            cen_j = cen_j / ((float) (cen_ct));
            cen_k = cen_k / ((float) (cen_ct));
            //printf("Baricenter1 = (%0.3f, %0.3f, %0.3f)\n", cen_i, cen_j, cen_k);

            tort_array[ I2(ct, 0, cc_array_numel) ] = (float) cen_i;
            tort_array[ I2(ct, 1, cc_array_numel) ] = (float) cen_j;
            tort_array[ I2(ct, 2, cc_array_numel) ] = (float) cen_k;
            tort_array[ I2(ct, 3, cc_array_numel) ] = (float) max_width;

            // Compute the length for further assessment of tortuosity by re-scanning bounding box:
            for (k = (curr_bb.min_z - offset); k <= (curr_bb.max_z + offset); k++)
                for (j = (curr_bb.min_y - offset); j <= (curr_bb.max_y + offset); j++)
                    for (i = (curr_bb.min_x - offset); i <= (curr_bb.max_x + offset); i++) {
                        if ((i >= 0) && (j >= 0) && (k >= 0) &&
                                (i < dimx) && (j < dimy) && (k < dimz)) {
                            // If a branch is found (USHRT_MAX voxel):
                            if (tmp_im2[ I(i, j, k, dimx, dimy) ] == USHRT_MAX) {
                                // Check if the branch voxel is connected to the pore ball:
                                if (_countNeighbors(tmp_im2, dimx, dimy, dimz, i, j, k,
                                        (unsigned short) (ct + 3)) >= 1) {
                                    // This counter is not coordination number because
                                    // the branches are removed after their identification
                                    // in order to speed up the computation:
                                    coord_array[ct]++;

                                    // At this point we know the first skeleton voxel "in touch"
                                    // with pore ball and we still know the baricenter, so let's
                                    // compute the euclidean distance from this two points since
                                    // this value could be greater than the size of the maximal ball:
                                    //printf("Dist1 = (%d, %d, %d)\n", i, j, k);
                                    dist1 = sqrt((i * 1.0 - cen_i)*(i * 1.0 - cen_i) +
                                            (j * 1.0 - cen_j)*(j * 1.0 - cen_j) +
                                            (k * 1.0 - cen_k)*(k * 1.0 - cen_k));

                                    // Now that we are close to the border with current node, go 
                                    // and run along the branch, removing it in order to avoid 
                                    // counting the considered segment more than once:			
                                    tmp_im2[ I(i, j, k, dimx, dimy) ] = BACKGROUND;

                                    cc_length = 0.0;

                                    // Go until the end:
                                    a = i;
                                    b = j;
                                    c = k;
                                    while (_findNeighbor(tmp_im2, dimx, dimy, dimz, a, b, c, &coords) != 0) {

                                        // Count the length with Euclidean criteria:
                                        cc_length += sqrt((double) ((a - coords.x)*(a - coords.x) +
                                                (b - coords.y)*(b - coords.y) +
                                                (c - coords.z)*(c - coords.z)));
                                        //printf("cc_length = %0.3f\n", cc_length );


                                        a = coords.x;
                                        b = coords.y;
                                        c = coords.z;
                                        tmp_im2[ I(a, b, c, dimx, dimy) ] = BACKGROUND;
                                    }

                                    // At this point (a,b,c) are the coordinates of the last point
                                    // "in touch" with the other pore ball. Again we have to compute 
                                    // the euclidean distance from this point and the baricenter of
                                    // the destination pore since this value could be greater than 
                                    // the size of the maximal ball. So let's first compute the baricenter
                                    // of the destination pore:

                                    // Get the destination pore:
                                    other_lbl = _findNode(tmp_im2, dimx, dimy, dimz, a, b, c) - 3;

                                    // Get its bounding box:
                                    curr_bb2 = bbs[other_lbl];

                                    cen_i = 0;
                                    cen_j = 0;
                                    cen_k = 0;
                                    cen_ct = 0;

                                    // Scan the bounding box of destination pore:
                                    for (t = (curr_bb2.min_z - offset); t <= (curr_bb2.max_z + offset); t++)
                                        for (s = (curr_bb2.min_y - offset); s <= (curr_bb2.max_y + offset); s++)
                                            for (r = (curr_bb2.min_x - offset); r <= (curr_bb2.max_x + offset); r++) {
                                                if ((r >= 0) && (s >= 0) && (t >= 0) &&
                                                        (r < dimx) && (s < dimy) && (t < dimz)) {
                                                    if (tmp_im[ I(r, s, t, dimx, dimy) ] == ((unsigned short) (other_lbl + 3))) {

                                                        // Compute also the center of mass of this pore:
                                                        cen_i += (float) r;
                                                        cen_j += (float) s;
                                                        cen_k += (float) t;
                                                        cen_ct++;
                                                    }
                                                }
                                            }

                                    cen_i = cen_i / ((float) (cen_ct));
                                    cen_j = cen_j / ((float) (cen_ct));
                                    cen_k = cen_k / ((float) (cen_ct));

                                    // Compute the Euclidean distance:
                                    //printf("Baricenter2 = (%0.3f, %0.3f, %0.3f)\n", cen_i, cen_j, cen_k);
                                    //printf("Dist2 = (%d, %d, %d)\n", a, b, c);
                                    dist2 = sqrt((a * 1.0 - cen_i)*(a * 1.0 - cen_i) +
                                            (b * 1.0 - cen_j)*(b * 1.0 - cen_j) +
                                            (c * 1.0 - cen_k)*(c * 1.0 - cen_k));



                                    cc_length += (dist1 + dist2 + 2*sqrt(3.0));
                                    //wr_log("\tfrom: %d, to: %d, length: %0.3f + 0.3f + 0.3f\n", ct, other_lbl, cc_length, dist1, dist2);
                                    //printf("from: %d, to: %d, length: %0.3f + %0.3f + %0.3f\n", ct, other_lbl, cc_length, dist1, dist2);

                                    if ( fabs(cc_length - 0.0) < 1E-4)
                                        cc_length = 1.0;

                                    // Now add cc_length to output matrix:		
                                    if ((other_lbl != -3) && (ct != other_lbl)) {
                                        if ( fabs(tort_matrix[ I2(ct, other_lbl, cc_array_numel) ] - 0.0) < 1E-4) {
                                            tort_matrix[ I2(ct, other_lbl, cc_array_numel) ] = (float) cc_length;
                                            //printf("%d - from: %d, to: %d, length: %0.3f\n", aa++, ct, other_lbl, cc_length);
                                        } else {
                                            if (cc_length < tort_matrix[ I2(ct, other_lbl, cc_array_numel) ]) {
                                                tort_matrix[ I2(ct, other_lbl, cc_array_numel) ] = (float) cc_length;
                                                //printf("%d - from: %d, to: %d, length: %0.3f\n", aa++, ct, other_lbl, cc_length);
                                            }
                                            //aa++;
                                        }
                                    }
                                }
                            }
                        }
                    }

        }

        // Correct tort_matrix and prepare it for Dijkstra:
        for (i = 0; i < cc_array_numel; i++) {
            for (j = 0; j < cc_array_numel; j++) {
                if (i == j) {
                    tort_matrix[ I2(i, j, cc_array_numel) ] = 0.0;
                } /*else {
                    if (abs(tort_matrix[ I2(i, j, cc_array_numel) ] - 0.0) < 1E-4) {
                        /*tort_matrix[ I2(i, j, cc_array_numel) ] += (float) sqrt(tort_array[ I2(i, 3, cc_array_numel) ] / 2.0);
                        tort_matrix[ I2(i, j, cc_array_numel) ] += (float) sqrt(tort_array[ I2(j, 3, cc_array_numel) ] / 2.0);
                        //tort_matrix[ I2( i, j, cc_array_numel) ]*= (float) 1.0;	
                    } else {
                        tort_matrix[ I2(i, j, cc_array_numel) ] = (float) (_DIJKSTRA_IN);
                    }
                }*/
            }
        }

        // Make symmetric:
        for (i = 0; i < cc_array_numel; i++) {
            for (j = 0; j < cc_array_numel; j++) {
                if ((j > i) && !((fabs(tort_matrix[ I2(i, j, cc_array_numel) ] - 0.0) < 1E-4))) {
                    tort_matrix[ I2(j, i, cc_array_numel) ] = tort_matrix[ I2(i, j, cc_array_numel) ];
                }
                if ((j < i) && !((fabs(tort_matrix[ I2(i, j, cc_array_numel) ] - 0.0) < 1E-4))) {
                    tort_matrix[ I2(j, i, cc_array_numel) ] = tort_matrix[ I2(i, j, cc_array_numel) ];
                }
            }
        }

        //printf("Array of dim %d:\n", cc_array_numel);
        /*for ( i = 0; i < cc_array_numel; i++ )
        {
        printf("%0.1f %0.1f %0.1f - %0.1f\n", tort_array [ I2(i, 0, cc_array_numel) ], 
                tort_array [ I2(i, 1, cc_array_numel) ], tort_array [ I2(i, 2, cc_array_numel) ], 
                tort_array [ I2(i, 3, cc_array_numel) ] );
        }
			
			
        fvol = fopen("R:\\TEMP\\matrix.dat", "wb");
        //printf("Matrix of %d x %d:\n", cc_array_numel, cc_array_numel);
        for ( i = 0; i < cc_array_numel; i++ )
        {
                for ( j = 0; j < cc_array_numel; j++ )
                {
                        if ( tort_matrix[ I2(i, j, cc_array_numel) ] == _DIJKSTRA_IN )
                        {
                                //printf("-1.0 ");
                                fprintf(fvol,"%0.1f ",0.0);
                                //fwrite(-1.0f, sizeof (float), 1, fvol);
                        }
                        else
                        {
                                //printf("%0.1f ");
                                fprintf(fvol,"%0.1f ",tort_matrix[ I2(i, j, cc_array_numel) ]);
                                //fwrite( (float) (tort_matrix[ I2(i, j, cc_array_numel) ]), sizeof (float), 1, fvol);
                        }
                }
                //printf("\n");	
                fprintf(fvol,"%s",";\n");
        }
        fclose(fvol);*/

        //out_stats->Tort_Z = _computeTortousity( tort_matrix, tort_array, cc_array_numel, 0, 5 );

        //
        // Compute actual tortuosity:
        //					

        // Check if tortousity depth is greater than the square of number of nodes:
        /*num_nodes = (int) sqrt( (float) (MIN(num_nodes*num_nodes, cc_array_numel)));
        out_stats->Tort_Counter = num_nodes*num_nodes;

        // Allocate memory:
        P3D_MEM_TRY(min_ct_array = (int*) calloc(num_nodes, sizeof (int)));
        P3D_MEM_TRY(max_ct_array = (int*) calloc(num_nodes, sizeof (int)));

        P3D_MEM_TRY(out_stats->Tort_X = (double*) calloc(out_stats->Tort_Counter, sizeof (double)));
        P3D_MEM_TRY(out_stats->Tort_Y = (double*) calloc(out_stats->Tort_Counter, sizeof (double)));
        P3D_MEM_TRY(out_stats->Tort_Z = (double*) calloc(out_stats->Tort_Counter, sizeof (double)));*/


        //
        // Turtousity along X:
        //

        // Get num_nodes minimum and maximum values along x (num_nodes is at
        // least 1 and no more than a threshold):
        /*for (ct_num_nodes = 0; ct_num_nodes < num_nodes; ct_num_nodes++) {
            min = (float) (INT_MAX);
            max = -1.0;

            // Scan the array:
            for (ct = 0; ct < cc_array_numel; ct++) {
                // Check if probable new min is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (min_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous min save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 0, cc_array_numel) ] < min) {
                        min = tort_array[ I2(ct, 0, cc_array_numel) ];
                        min_ct_array[ct_num_nodes] = ct;
                    }
                }
                // Check if probable new max is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (max_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous max save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 0, cc_array_numel) ] > max) {

                        max = tort_array[ I2(ct, 0, cc_array_numel) ];
                        max_ct_array[ct_num_nodes] =  ct;
                    }
                }
            }
        }

        // Compute tortousity in X:		
        ct = 0;
        for (i = 0; i < num_nodes; i++) {
            for (j = 0; j < num_nodes; j++) {
                out_stats->Tort_X[ct++] = _computeTortousity(tort_matrix, tort_array,
                        cc_array_numel, min_ct_array[i], max_ct_array[j]);
            }
        }



        //
        // Tortousity along Y:
        //

        // Get num_nodes minimum and maximum values along x (num_nodes is at
        // least 1 and no more than a threshold):
        // Get num_nodes minimum and maximum values along x (num_nodes is at
        // least 1 and no more than a threshold):
        for (ct_num_nodes = 0; ct_num_nodes < num_nodes; ct_num_nodes++) {
            min = (float) (INT_MAX);
            max = -1.0;

            // Scan the array:
            for (ct = 0; ct < cc_array_numel; ct++) {
                // Check if probable new min is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (min_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous min save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 1, cc_array_numel) ] < min) {
                        min = tort_array[ I2(ct, 1, cc_array_numel) ];
                        min_ct_array[ct_num_nodes] = ct;
                    }
                }
                // Check if probable new max is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (max_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous max save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 1, cc_array_numel) ] > max) {
                        max = tort_array[ I2(ct, 1, cc_array_numel) ];
                        max_ct_array[ct_num_nodes] = ct;
                    }
                }
            }
        }

        // Compute tortousity in Y:		
        ct = 0;
        for (i = 0; i < num_nodes; i++)
            for (j = 0; j < num_nodes; j++) {
                out_stats->Tort_Y[ct++] = _computeTortousity(tort_matrix, tort_array,
                        cc_array_numel, min_ct_array[i], max_ct_array[j]);
            }


        //
        // Tortousity along Z:
        //

        // Get num_nodes minimum and maximum values along Z (num_nodes is at
        // least 1 and no more than a threshold):
        for (ct_num_nodes = 0; ct_num_nodes < num_nodes; ct_num_nodes++) {
            min = (float) (INT_MAX);
            max = -1.0;

            // Scan the array:
            for (ct = 0; ct < cc_array_numel; ct++) {
                // Check if probable new min is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (min_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous min save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 2, cc_array_numel) ] < min) {
                        min = tort_array[ I2(ct, 2, cc_array_numel) ];
                        min_ct_array[ct_num_nodes] = ct;
                    }
                }
                // Check if probable new max is actually one the previous mins:
                flag = P3D_FALSE;
                for (k = 0; k <= (ct_num_nodes - 1); k++)
                    if (max_ct_array[k] == ct)
                        flag = P3D_TRUE;
                // If it's not a previous max save it:
                if (flag == P3D_FALSE) {
                    if (tort_array[ I2(ct, 2, cc_array_numel) ] > max) {
                        max = tort_array[ I2(ct, 2, cc_array_numel) ];
                        max_ct_array[ct_num_nodes] = ct;
                    }
                }
            }
        }


        // Compute tortousity in Z:		
        ct = 0;
        for (i = 0; i < num_nodes; i++)
            for (j = 0; j < num_nodes; j++)
                out_stats->Tort_Z[ct++] = _computeTortousity(tort_matrix, tort_array,
                    cc_array_numel, min_ct_array[i], max_ct_array[j]);
		*/

    } else {
        // Copy empty statistics to output structure:
        out_stats->Node_Counter = 0;
        out_stats->Node_Width = NULL;
        out_stats->ConnectivityDensity = 0.0;

       /* out_stats->Tort_Counter = 0;
        out_stats->Tort_X = NULL;
        out_stats->Tort_Y = NULL;
        out_stats->Tort_Z = NULL;*/
    }

    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);
    if (poresplusend_im != NULL) free(poresplusend_im);
    if (tort_matrix != NULL) free(tort_matrix);
    if (tort_array != NULL) free(tort_array);
    if (min_ct_array != NULL) free(min_ct_array);
    if (max_ct_array != NULL) free(max_ct_array);
    if (coord_array != NULL) free(coord_array);

    return P3D_SUCCESS;

MEM_ERROR:

    // Release resources:
    if (cc_array != NULL) free(cc_array);
    if (bbs != NULL) free(bbs);
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);
    if (poresplusend_im != NULL) free(poresplusend_im);
    if (tort_matrix != NULL) free(tort_matrix);
    if (tort_array != NULL) free(tort_array);
    if (min_ct_array != NULL) free(min_ct_array);
    if (max_ct_array != NULL) free(max_ct_array);
    if (coord_array != NULL) free(coord_array);

    return P3D_ERROR;
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
        const int tortuosity_depth,
        const double voxelsize, // IN: voxel resolution
        int (*wr_log)(const char*, ...)
        ) 
{

    // Temporary matrices:
    unsigned char* max_skl_im;
    unsigned char* lbl_skl_im;
    unsigned short* dt_im;

    double mean = 0.0;
    double mean_sqr = 0.0;
    double glob_mean = 0.0;
    double glob_mean_sqr = 0.0;
    int /*i,j,*/ct;
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
        //wr_log("\tTortuosity depth: %d.", tortuosity_depth);
    }


    // Allocate memory:
    if (nodes_im == NULL) {
        flag_nodes_null = P3D_TRUE;
        P3D_MEM_TRY(nodes_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (pores_im == NULL) {
        flag_pores_null = P3D_TRUE;
        P3D_MEM_TRY(pores_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (ends_im == NULL) {
        flag_ends_null = P3D_TRUE;
        P3D_MEM_TRY(ends_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }
    if (throats_im == NULL) {
        flag_throats_null = P3D_TRUE;
        P3D_MEM_TRY(throats_im = (unsigned char*) calloc(dimx * dimy*dimz, sizeof (unsigned char)));
    }

    P3D_MEM_TRY(dt_im = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));
    P3D_MEM_TRY(max_skl_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));
    P3D_MEM_TRY(lbl_skl_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char)));


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

        /*if (out_stats->NodeToNode_Counter > 0) {

            // Compute mean values of tortuosity:
            glob_mean = 0.0;
            glob_mean_sqr = 0.0;
            mean = 0.0;
            mean_sqr = 0.0;
            i = 0;
            j = 0;
            for (ct = 0; ct < (out_stats->Tort_Counter); ct++) {
                if (out_stats->Tort_X[ct] > -1.0) {
                    mean = mean + (double) (out_stats->Tort_X[ct]);
                    glob_mean += (double) (out_stats->Tort_X[ct]);
                    mean_sqr += (double) (out_stats->Tort_X[ct] * out_stats->Tort_X[ct]);
                    glob_mean_sqr += (double) (out_stats->Tort_X[ct] * out_stats->Tort_X[ct]);
                    i++;
                    j++;
                } else {
                    //wr_log("\tWARNING: Unable to compute a tortuosity path.");
                }
            }
            mean = mean / i;
            mean_sqr = mean_sqr / i;
            mean_sqr = sqrt(mean_sqr - mean * mean);
            wr_log("\tTortuosity in X: %0.3f +/- %0.3f [-] (%d considered paths).", mean, mean_sqr, i);

            mean = 0.0;
            mean_sqr = 0.0;
            i = 0;
            for (ct = 0; ct < (out_stats->Tort_Counter); ct++) {
                if (out_stats->Tort_Y[ct] > -1.0) {
                    mean = mean + (double) (out_stats->Tort_Y[ct]);
                    glob_mean += (double) (out_stats->Tort_Y[ct]);
                    mean_sqr += (double) (out_stats->Tort_Y[ct] * out_stats->Tort_Y[ct]);
                    glob_mean_sqr += (double) (out_stats->Tort_Y[ct] * out_stats->Tort_Y[ct]);
                    i++;
                    j++;
                } else {
                    //wr_log("\tWARNING: Unable to compute a tortuosity path.");
                }
            }
            mean = mean / i;
            mean_sqr = mean_sqr / i;
            mean_sqr = sqrt(mean_sqr - mean * mean);
            wr_log("\tTortuosity in Y: %0.3f +/- %0.3f [-] (%d considered paths).", mean, mean_sqr, i);


            mean = 0.0;
            mean_sqr = 0.0;
            i = 0;
            for (ct = 0; ct < (out_stats->Tort_Counter); ct++) {
                if (out_stats->Tort_Z[ct] > -1.0) {
                    mean = mean + (double) (out_stats->Tort_Z[ct]);
                    glob_mean += (double) (out_stats->Tort_Z[ct]);
                    mean_sqr += (double) (out_stats->Tort_Z[ct] * out_stats->Tort_Z[ct]);
                    glob_mean_sqr += (double) (out_stats->Tort_Z[ct] * out_stats->Tort_Z[ct]);
                    i++;
                    j++;
                } else {
                    //wr_log("\tWARNING: Unable to compute a tortuosity path.");
                }
            }
            mean = mean / i;
            mean_sqr = mean_sqr / i;
            mean_sqr = sqrt(mean_sqr - mean * mean);
            glob_mean = glob_mean / j;
            glob_mean_sqr = glob_mean_sqr / j;
            glob_mean_sqr = sqrt(glob_mean_sqr - glob_mean * glob_mean);
            wr_log("\tTortuosity in Z: %0.3f +/- %0.3f [-] (%d considered paths).", mean, mean_sqr, i);
            wr_log("\tGlobal tortuosity: %0.3f +/- %0.3f [-] (%d considered paths).", glob_mean, glob_mean_sqr, j);
        }*/
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
    return P3D_ERROR;

}


