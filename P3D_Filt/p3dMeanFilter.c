#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

int p3dMeanFilter2D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned char* tmp_im;

    int a_dimx, a_dimy;
    int i, j;
    int x, y;
    int pr;

    // Variables for computing gaussian kernel:
    int a_rad;
    double sum;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying mean filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;

    // Initialize input:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * sizeof (unsigned char)));
    P3D_TRY(p3dReplicatePadding2D_8(in_im, tmp_im, dimx, dimy, a_rad, NULL, NULL));

    pr = 0;

    // Volume scanning:
#pragma omp parallel for private(i, x, y, sum) reduction( + : pr)
    for (j = a_rad; j < (a_dimy - a_rad); j++) {
        for (i = a_rad; i < (a_dimx - a_rad); i++) {
            sum = 0;
            // Fill temporary array:
            for (y = (j - a_rad); y <= (j + a_rad); y++)
                for (x = (i - a_rad); x <= (i + a_rad); x++)
                    sum += (double) tmp_im[ I2(x, y, a_dimx) ];

            // Set out voxel with the mean of the sorted temporary array:
            out_im[ I2(i - a_rad, j - a_rad, dimx) ] = (unsigned char) (sum / (size * size) + 0.5);

            // Increase progress counter:
            pr++;
        }

        // Update any progress bar:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy)*100 + 0.5));
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Mean filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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
    return P3D_MEM_ERROR;
}

int p3dMeanFilter2D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned short* tmp_im;

    int a_dimx, a_dimy;
    int i, j;
    int x, y;
    int pr;

    // Variables for computing gaussian kernel:
    int a_rad;
    double sum;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying mean filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;

    // Initialize input:
    P3D_TRY(tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * sizeof (unsigned short)));
    P3D_TRY(p3dReplicatePadding2D_16(in_im, tmp_im, dimx, dimy, a_rad, NULL, NULL));

    pr = 0;

    // Volume scanning:
#pragma omp parallel for private(i, x, y, sum) reduction( + : pr)
    for (j = a_rad; j < (a_dimy - a_rad); j++) {
        for (i = a_rad; i < (a_dimx - a_rad); i++) {
            sum = 0;
            // Fill temporary array:
            for (y = (j - a_rad); y <= (j + a_rad); y++)
                for (x = (i - a_rad); x <= (i + a_rad); x++)
                    sum += (double) tmp_im[ I2(x, y, a_dimx) ];

            // Set out voxel with the mean of the sorted temporary array:
            out_im[ I2(i - a_rad, j - a_rad, dimx) ] = (unsigned short) (sum / (size * size) + 0.5);

            // Increase progress counter:
            pr++;
        }

        // Update any progress bar:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy)*100 + 0.5));
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Mean filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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
    return P3D_MEM_ERROR;
}

int p3dMeanFilter3D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned char* tmp_im;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int pr;

    // Variables for computing gaussian kernel:
    int a_rad;
    double sum;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying mean filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Initialize input:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
    P3D_TRY(p3dReplicatePadding3D_8(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));

    pr = 0;

    // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, sum) reduction( + : pr)
    for (k = a_rad; k < (a_dimz - a_rad); k++) {
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0;
                // Fill temporary array:
                for (z = (k - a_rad); z <= (k + a_rad); z++)
                    for (y = (j - a_rad); y <= (j + a_rad); y++)
                        for (x = (i - a_rad); x <= (i + a_rad); x++)
                            sum += (double) tmp_im[ I(x, y, z, a_dimx, a_dimy) ];

                // Set out voxel with the mean of the sorted temporary array:
                out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = (unsigned char)
                        (sum / (size * size * size));

                // Increment progress counter:
                pr++;
            }

        // Update any progress counter:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy * dimz)*100 + 0.5));
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Mean filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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
    return P3D_MEM_ERROR;
}

int p3dMeanFilter3D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    unsigned short* tmp_im;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int pr;

    // Variables for computing gaussian kernel:
    int a_rad;
    double sum;


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying mean filter...");
        wr_log("\tKernel size: %d.", size);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Initialize input:
    P3D_TRY(tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned short)));
    P3D_TRY(p3dReplicatePadding3D_16(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));

    pr = 0;

    // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, sum) reduction( + : pr)
    for (k = a_rad; k < (a_dimz - a_rad); k++) {
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0;
                // Fill temporary array:
                for (z = (k - a_rad); z <= (k + a_rad); z++)
                    for (y = (j - a_rad); y <= (j + a_rad); y++)
                        for (x = (i - a_rad); x <= (i + a_rad); x++)
                            sum += (double) tmp_im[ I(x, y, z, a_dimx, a_dimy) ];

                // Set out voxel with the mean of the sorted temporary array:
                out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = (unsigned short)
                        (sum / (size * size * size));

                // Increment progress counter:
                pr++;
            }

        // Update any progress counter:
        if (wr_progress != NULL) wr_progress((int) ((double) (pr) / (dimx * dimy * dimz)*100 + 0.5));
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Mean filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
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
    return P3D_MEM_ERROR;
}

