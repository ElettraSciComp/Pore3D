#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

int p3dGaussianFilter3D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        const double sigma,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    float* tmp_im = NULL;
    float* kernel = NULL;


    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int ct, a_rad, a_size;
    double d;
    double f, w, sum, sum_w;
    
    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dGaussianFilter3D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying gaussian filter...");
        wr_log("\tKernel size: %d.", size);
        wr_log("\tSigma: %0.3f.", sigma);
    }

    // Set kernel size and variance:
    if (size < 1.0)
        a_rad = 1;
    else
        a_rad = ceil(size / 2);
    a_size = 2 * a_rad + 1;

    d = (2.0 * sigma * sigma);


    // Compute dimensions of padded volume:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;


    // Initialize input:
    tmp_im = (float*) malloc(a_dimx * a_dimy * a_dimz * sizeof (float));
    kernel = (float*) malloc(a_size * sizeof (float));

    _p3dReplicatePadding3D_uchar2float(in_im, tmp_im, dimx, dimy, dimz, a_rad);

    ct = 0;
    sum_w = 0.0;

    // Create gausssian kernel:
    for (x = (-a_rad); x <= a_rad; x++) {
       
       f = (x)*(x) / d;
       kernel[ct] = exp(-f);
       sum_w = sum_w + kernel[ct++];
    }
    ct = 0;
    // Normalize kernel
    for (x = (-a_rad); x <= a_rad; x++) 
        kernel[ct] = kernel[ct++]/sum_w;

    // X-direction scanning
#pragma omp parallel for private(i, j, x, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (x = (i - a_rad); x <= (i + a_rad); x++) {
                    // Gausssian kernel:                       
                    sum = sum + kernel[ct++] * tmp_im[ I(x, j, k, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array:                
                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = sum;
            }

    // Y-direction scanning
#pragma omp parallel for private(i, j, y, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (y = (j - a_rad); y <= (j + a_rad); y++) {
                    // Gausssian kernel:
                    sum = sum + kernel[ct++] * tmp_im[ I(i, y, k, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array:             
                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = sum;
            }

    // Z-direction scanning
#pragma omp parallel for private(i, j, z, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (z = (k - a_rad); z <= (k + a_rad); z++) {
                    // Gausssian kernel:
                    sum = sum + kernel[ct++] * tmp_im[ I(i, j, z, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array: 
                if (sum < 0)
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = 0;
                else if (sum > UCHAR_MAX)
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = UCHAR_MAX;
                else
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = (unsigned char) sum;

            }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Gaussian filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (kernel != NULL) free(kernel);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (kernel != NULL) free(kernel);

    // Return error:
    return P3D_MEM_ERROR;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}

int p3dGaussianFilter3D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        const double sigma,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Padded input and related dims:
    float* tmp_im = NULL;
    float* kernel = NULL;


    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int ct, a_rad, a_size;
    double d;
    double f, w, sum, sum_w;
    
    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dGaussianFilter3D_16");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying gaussian filter...");
        wr_log("\tKernel size: %d.", size);
        wr_log("\tSigma: %0.3f.", sigma);
    }

    // Set kernel size and variance:
    if (size < 1.0)
        a_rad = 1;
    else
        a_rad = ceil(size / 2);
    a_size = 2 * a_rad + 1;

    d = (2.0 * sigma * sigma);


    // Compute dimensions of padded volume:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;


    // Initialize input:
    tmp_im = (float*) malloc(a_dimx * a_dimy * a_dimz * sizeof (float));
    kernel = (float*) malloc(a_size * sizeof (float));

    _p3dReplicatePadding3D_ushort2float(in_im, tmp_im, dimx, dimy, dimz, a_rad);

    ct = 0;
    sum_w = 0.0;

    // Create gausssian kernel:
    for (x = (-a_rad); x <= a_rad; x++) {
       
       f = (x)*(x) / d;
       kernel[ct] = exp(-f);
       sum_w = sum_w + kernel[ct++];
    }
    ct = 0;
    // Normalize kernel
    for (x = (-a_rad); x <= a_rad; x++) 
        kernel[ct] = kernel[ct++]/sum_w;

    // X-direction scanning
#pragma omp parallel for private(i, j, x, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (x = (i - a_rad); x <= (i + a_rad); x++) {
                    // Gausssian kernel:                       
                    sum = sum + kernel[ct++] * tmp_im[ I(x, j, k, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array:                
                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = sum;
            }

    // Y-direction scanning
#pragma omp parallel for private(i, j, y, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (y = (j - a_rad); y <= (j + a_rad); y++) {
                    // Gausssian kernel:
                    sum = sum + kernel[ct++] * tmp_im[ I(i, y, k, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array:             
                tmp_im[ I(i, j, k, a_dimx, a_dimy) ] = sum;
            }

    // Z-direction scanning
#pragma omp parallel for private(i, j, z, sum, ct)
    for (k = a_rad; k < (a_dimz - a_rad); k++)
        for (j = a_rad; j < (a_dimy - a_rad); j++)
            for (i = a_rad; i < (a_dimx - a_rad); i++) {
                sum = 0.0;
                ct = 0;

                // Process kernel:
                for (z = (k - a_rad); z <= (k + a_rad); z++) {
                    // Gausssian kernel:
                    sum = sum + kernel[ct++] * tmp_im[ I(i, j, z, a_dimx, a_dimy) ];
                }

                // Set out voxel with the mean of the sorted temporary array:               
                if (sum < 0)
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = 0;
                else if (sum > USHRT_MAX)
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = USHRT_MAX;
                else
                    out_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] = (unsigned short) sum;
            }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Gaussian filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (kernel != NULL) free(kernel);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (kernel != NULL) free(kernel);

    // Return error:
    return P3D_MEM_ERROR;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/

}
