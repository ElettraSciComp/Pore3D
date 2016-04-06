#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

/// To be implemented for square, rectangular, circle.... At first: circle

int p3dCreateBinaryCircle(
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int cenx,
        const int ceny,
        const int radius,
        int (*wr_log)(const char*, ...)
        ) {
    unsigned int i, j;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dCreateBinaryCircle");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Log a message:
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Creating binary circle...");
        wr_log("\tCenter: [%d, %d].", cenx, ceny);
        wr_log("\tRadius: %d.", radius);
    }


    // Remove connected component:
    for (j = 0; j < dimy; j++)
        for (i = 0; i < dimx; i++) {
            // Labels start from 3:
            if (sqrt((double) ((i - cenx)*(i - cenx)+(j - ceny)*(j - ceny))) <= radius) {
                out_im[ I2(i, j, dimx) ] = OBJECT;
            } else {
                out_im[ I2(i, j, dimx) ] = BACKGROUND;
            }
        }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Binary circle created successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dCreateBinaryCylinder(
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int cenx,
        const int ceny,
        const int radius,
        int (*wr_log)(const char*, ...)
        ) {
    int i, j, k;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dCreateBinaryCylinder");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Log a message:
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Creating binary cylinder...");
        wr_log("\tCenter: [%d, %d].", cenx, ceny);
        wr_log("\tRadius: %d.", radius);
    }


    // Remove connected component:
#pragma omp parallel for private(i, j)
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                if (sqrt((double) ((i - cenx)*(i - cenx)+(j - ceny)*(j - ceny))) <= radius) {
                    out_im[ I(i, j, k, dimx, dimy) ] = OBJECT;
                } else {
                    out_im[ I(i, j, k, dimx, dimy) ] = BACKGROUND;
                }
            }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Binary cylinder created successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dCreateBinarySphere(
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int cenx,
        const int ceny,
        const int cenz,
        const int radius,
        int (*wr_log)(const char*, ...)
        ) {
    int i, j, k;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dCreateBinarySphere");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Log a message:
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Creating binary sphere...");
        wr_log("\tCenter: [%d, %d, %d].", cenx, ceny, cenz);
        wr_log("\tRadius: %d.", radius);
    }


    // Remove connected component:
#pragma omp parallel for private(i, j)
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                if (sqrt((double) ((i - cenx)*(i - cenx)+(j - ceny)*(j - ceny)+(k - cenz)*(k - cenz))) <= radius) {
                    out_im[ I(i, j, k, dimx, dimy) ] = OBJECT;
                } else {
                    out_im[ I(i, j, k, dimx, dimy) ] = BACKGROUND;
                }
            }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Binary sphere created successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}







