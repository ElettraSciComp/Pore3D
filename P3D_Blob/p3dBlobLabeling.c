#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dBlob.h"
#include "p3dTime.h"
#include "p3dAuth.h"

#include "Common/p3dConnectedComponentsLabeling.h"
#include "Common/p3dUtils.h"

// Wrapper for the extern:

int p3dBlobLabeling_ushort(
        unsigned char* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int conn,
        const int random_lbl, // Flag for random labels
        const int skip_borders,
        int (*wr_log)(const char*, ...)
        ) {

    //char auth_code;

    //
    // Authenticate:
    //
    //auth_code = authenticate("p3dBlobLabeling_ushort");
    //if (auth_code == '0') goto AUTH_ERROR;

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

    // Return OK:
    return P3D_MEM_ERROR;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dBlobLabeling_uint(
        unsigned char* in_im,
        unsigned int* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int conn,
        const int random_lbl, // Flag for random labels
        const int skip_borders,
        int (*wr_log)(const char*, ...)
        ) {

    //char auth_code;

    //
    // Authenticate:
    //
    //auth_code = authenticate("p3dBlobLabeling_uint");
    //if (auth_code == '0') goto AUTH_ERROR;
    //if (strcmp(auth_code, "1") != 0) goto AUTH_ERROR;

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


    P3D_TRY( p3dConnectedComponentsLabeling_uint(in_im, out_im, NULL, NULL, NULL, dimx, dimy, dimz,
            conn, random_lbl, skip_borders) );

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

    // Return OK:
    return P3D_MEM_ERROR;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

