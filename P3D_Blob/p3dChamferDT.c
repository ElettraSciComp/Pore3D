/* 
 * p3dChamferDistanceTransform  
 * 
 * Computes the Chamfer distance transform of the input volume. For each 
 * non-zero voxel in input image, the distance transform assigns a number 
 * that is the distance between that voxel and the nearest zero pixel 
 * (background). The Chamfer distance transform is an integer approximation 
 * for the Euclidean distance metric.
 *
 * Remarks
 * -------
 * Current implementation uses <1,1,1> weights. See [1] for details.
 * Computational cost is O(N^3) for a NxNxN input volume. Memory 
 * requirement is two extra volumes of a NxNxN, i.e. O(N^3). Current 
 * implementation uses unsigned short elements (see compiler specifications).
 *
 * References
 * ----------
 * [1] Gunilla Borgefors. "On Digital Distance Transforms in Three Dimensions", 
 * Computer Vision and Image Understanding, Vol. 64, No. 3, pp. 368-376, 1996.
 *

/****/
/* Separare i tre casi "cityblock", "chessboard" e "euclidean". C'e da usare
 la connettivit� 6 per cityblock riscrivendo le funzioni in modo da tenere 
 conto solo della connettivit� 6.
 */
/***/


/* 
 * Copyright 2008, SYRMEP Group - Sincrotrone Trieste S.C.p.A.
 *
 * Author:	Brun Francesco
 * Version:	1.0
 * Date:		19 june 2008
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dBlob.h"
#include "p3dTime.h"
#include "p3dAuth.h"


#define INF USHRT_MAX		// infinity

unsigned short minOfPrevious13(
        unsigned char* in_rev,
        unsigned short* out_rev,
        unsigned short* weights,
        const int x,
        const int y,
        const int z,
        const int dimx,
        const int dimy,
        const int dimz
        ) {
    int a, b, c;

    // Temporary values:
    unsigned short tmp, min;

    // Initialize current minimum:
    min = INF;


    // Computes the sums of the values of the 13 already visited
    // neighbors (in direct order) and the corresponding local 
    // distances:	

    // From 0 to 8:
    c = z - 1;
    for (b = (y - 1); b <= (y + 1); b++)
        for (a = (x - 1); a <= (x + 1); a++) {
            tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                    weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

            if (tmp < min) {
                min = tmp;
            }
        }
    // From 9 to 11:
    c = z;
    b = y - 1;
    for (a = (x - 1); a <= (x + 1); a++) {
        tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

        if (tmp < min) {
            min = tmp;
        }
    }
    // The 12th:
    c = z;
    b = y;
    a = x - 1;


    tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
            weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

    if (tmp < min) {
        min = tmp;
    }

    // Return the minimum of the 13 sums:
    return min;
}

unsigned short minOfPrevious14(
        unsigned short* out_rev,
        unsigned short* weights,
        const int x,
        const int y,
        const int z,
        const int dimx,
        const int dimy,
        const int dimz
        ) {
    int a, b, c;

    // Temporary values:
    unsigned short tmp, min;


    // Initialize current minimum:
    min = out_rev[ I(x, y, z, dimx, dimy) ];

    // Computes the sums of the values of the 14 already visited
    // neighbors (in reverse order plus the center voxel) and the
    // corresponding local distances:	

    // From 0 to 8:
    c = z + 1;
    for (b = (y + 1); b >= (y - 1); b--)
        for (a = (x + 1); a >= (x - 1); a--) {
            tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                    weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

            if (tmp < min) {
                min = tmp;
            }
        }
    // From 9 to 11:
    c = z;
    b = y + 1;
    for (a = (x + 1); a >= (x - 1); a--) {
        tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];


        if (tmp < min) {
            min = tmp;
        }
    }
    // The 12th voxel:
    c = z;
    b = y;
    a = x + 1;


    tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
            weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

    if (tmp < min) {
        min = tmp;
    }


    // The center voxel (without a sum):
    if (out_rev[ I(x, y, z, dimx, dimy) ] < min) {
        min = tmp;
    }

    // Return minimum:
    return min;

}

unsigned short minOfPrevious13_bt(
        unsigned char* in_rev,
        unsigned short* out_rev,
        unsigned short* weights,
        const int x,
        const int y,
        const int z,
        const int dimx,
        const int dimy,
        const int dimz
        ) {
    int a, b, c;

    // Temporary values:
    unsigned short tmp, min;

    // Initialize current minimum:
    min = INF;


    // Computes the sums of the values of the 13 already visited
    // neighbors (in direct order) and the corresponding local 
    // distances:	

    // From 0 to 8:
    c = z - 1;
    for (b = (y - 1); b <= (y + 1); b++)
        for (a = (x - 1); a <= (x + 1); a++) {
            // If inside volume margins:
            if ((a >= 0) && (a < dimx) &&
                    (b >= 0) && (b < dimy) &&
                    (c >= 0) && (c < dimz)) {
                tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                        weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

                if (tmp < min) {
                    min = tmp;
                }
            }
        }
    // From 9 to 11:
    c = z;
    b = y - 1;
    for (a = (x - 1); a <= (x + 1); a++) {
        // If inside volume margins:
        if ((a >= 0) && (a < dimx) &&
                (b >= 0) && (b < dimy) &&
                (c >= 0) && (c < dimz)) {
            tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                    weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

            if (tmp < min) {
                min = tmp;
            }
        }
    }
    // The 12th:
    c = z;
    b = y;
    a = x - 1;

    // If inside volume margins:
    if ((a >= 0) && (a < dimx) &&
            (b >= 0) && (b < dimy) &&
            (c >= 0) && (c < dimz)) {
        tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

        if (tmp < min) {
            min = tmp;
        }
    }

    // Return the minimum of the 13 sums:
    return min;
}

unsigned short minOfPrevious14_bt(
        unsigned short* out_rev,
        unsigned short* weights,
        const int x,
        const int y,
        const int z,
        const int dimx,
        const int dimy,
        const int dimz
        ) {
    int a, b, c;

    // Temporary values:
    unsigned short tmp, min;


    // Initialize current minimum:
    min = out_rev[ I(x, y, z, dimx, dimy) ];

    // Computes the sums of the values of the 14 already visited
    // neighbors (in reverse order plus the center voxel) and the
    // corresponding local distances:	

    // From 0 to 8:
    c = z + 1;
    for (b = (y + 1); b >= (y - 1); b--)
        for (a = (x + 1); a >= (x - 1); a--) {
            // If inside volume margins:
            if ((a >= 0) && (a < dimx) &&
                    (b >= 0) && (b < dimy) &&
                    (c >= 0) && (c < dimz)) {
                tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                        weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

                if (tmp < min) {
                    min = tmp;
                }
            }
        }
    // From 9 to 11:
    c = z;
    b = y + 1;
    for (a = (x + 1); a >= (x - 1); a--) {
        // If inside volume margins:
        if ((a >= 0) && (a < dimx) &&
                (b >= 0) && (b < dimy) &&
                (c >= 0) && (c < dimz)) {
            tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                    weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];


            if (tmp < min) {
                min = tmp;
            }
        }
    }
    // The 12th voxel:
    c = z;
    b = y;
    a = x + 1;

    // If inside volume margins:
    if ((a >= 0) && (a < dimx) &&
            (b >= 0) && (b < dimy) &&
            (c >= 0) && (c < dimz)) {
        tmp = out_rev[ I(a, b, c, dimx, dimy) ] +
                weights[ I((a - x) + 1, (b - y) + 1, (c - z) + 1, 3, 3) ];

        if (tmp < min) {
            min = tmp;
        }
    }

    // The center voxel (without a sum):
    if (out_rev[ I(x, y, z, dimx, dimy) ] < min) {
        min = tmp;
    }

    // Return minimum:
    return min;

}

/**
 * p3dChamferDistanceTransform
 *
 * Computes the Chamfer distance transform of the input volume.
 *
 */
int p3dChamferDT(
        unsigned char* in_rev,
        unsigned short* out_rev,
        const int dimx,
        const int dimy,
        const int dimz,
        const int w1,
        const int w2,
        const int w3,
        int (*wr_log)(const char*, ...)
        ) {

    // Counters:
    int i;
    int x, y, z;


    // Temporary values:
    unsigned short weights[27];

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dChamferDT");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Computing Chamfer distance transform...");
        wr_log("\tWeights: [%d,%d,%d].", w1, w2, w3);
    }

    // Apply algorithm:

    i = 0;
    // z = -1 :
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;
    // z = 0 : 
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) 0;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) w2;
    // z = 1 :
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w1;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;
    weights[i++] = (unsigned short) w2;
    weights[i++] = (unsigned short) w3;


    // Start forward scanning of padded volume:
    for (z = 0; z < dimz; z++)
        for (y = 0; y < dimy; y++)
            for (x = 0; x < dimx; x++) {
                // If input REV is not zero:
                if (in_rev[ I(x, y, z, dimx, dimy) ]) {
                    // For performance reasons a distinction is made between
                    // internal voxels and border voxels:
                    if ((x >= 1) && (x < (dimx - 1)) &&
                            (y >= 1) && (y < (dimy - 1)) &&
                            (z >= 1) && (z < (dimz - 1))) {
                        // Computes the sums of the values of the 13 already visited
                        // neighbors (in direct order) and the corresponding local 
                        // distances. The new voxel value is the minimum of the 13 sums:
                        out_rev[ I(x, y, z, dimx, dimy) ] =
                                minOfPrevious13(in_rev, out_rev, weights, x, y, z, dimx, dimy, dimz);

                    } else {
                        // The new voxel value is the minimum of the sums of 
                        // neighborhood already visited voxels within volume:
                        out_rev[ I(x, y, z, dimx, dimy) ] =
                                minOfPrevious13_bt(in_rev, out_rev, weights, x, y, z, dimx, dimy, dimz);
                    }
                } else {
                    // Set output voxel to zero:
                    out_rev[ I(x, y, z, dimx, dimy) ] = 0;
                }
            }


    // Start backward scanning of 1-voxel zero-padded volume:
    for (z = (dimz - 1); z >= 0; z--)
        for (y = (dimy - 1); y >= 0; y--)
            for (x = (dimx - 1); x >= 0; x--) {
                // If input REV (output REV of previous step) is not zero:
                if (out_rev[ I(x, y, z, dimx, dimy) ]) {
                    // For performance reasons a distinction is made between
                    // internal voxels and border voxels:
                    if ((x >= 1) && (x < (dimx - 1)) &&
                            (y >= 1) && (y < (dimy - 1)) &&
                            (z >= 1) && (z < (dimz - 1))) {
                        // Computes the sums of the values of the 14 already visited
                        // neighbors (in reverse order plus the center voxel) and the
                        // corresponding local distances. The new voxel value is the 
                        // minimum of the 14 sums:
                        out_rev[ I(x, y, z, dimx, dimy) ] =
                                minOfPrevious14(out_rev, weights, x, y, z, dimx, dimy, dimz);
                    } else {
                        // The new voxel value is the minimum of the sum of the 14 
                        // already visited neighbors voxels within volumes:
                        out_rev[ I(x, y, z, dimx, dimy) ] =
                                minOfPrevious14_bt(out_rev, weights, x, y, z, dimx, dimy, dimz);
                    }
                }
            }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Chamfer distance transform computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

