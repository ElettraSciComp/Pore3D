#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"

int p3dBilateralFilter3D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        const double sigma_d,
        const double sigma_r,
        const int iter,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Temporary padded input-volume:
    unsigned char* tmp_im;
    unsigned char* tmp_im2;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int ct;

    // Variables for computing gaussian kernel:
    int a_rad;
    double tmp;

    // Variables for filter management:
    double w, sum_f, sum_fi;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dBilateralFilter3D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying bilateral filter...");
        wr_log("\tKernel size: %d.", size);
        wr_log("\tDomain sigma: %0.3f.", sigma_d);
        wr_log("\tRange sigma: %0.3f.", sigma_r);
        wr_log("\tNumber of iterations: %d.", iter);
    }


    // Init variables:
    a_rad = size / 2; // integer division

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Try to allocate memory:
    P3D_TRY(tmp_im = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
    P3D_TRY(tmp_im2 = (unsigned char*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned char)));
    P3D_TRY(p3dReplicatePadding3D_8(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));

    for (ct = 0; ct < iter; ct++) {
        // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, sum_f, sum_fi, w)
        for (k = a_rad; k < (a_dimz - a_rad); k++)
            for (j = a_rad; j < (a_dimy - a_rad); j++)
                for (i = a_rad; i < (a_dimx - a_rad); i++) {
                    // Init variables:
                    sum_f = 0.0;
                    sum_fi = 0.0;

                    // Convolve (i,j,k) voxel:
                    for (z = (k - a_rad); z <= (k + a_rad); z++)
                        for (y = (j - a_rad); y <= (j + a_rad); y++)
                            for (x = (i - a_rad); x <= (i + a_rad); x++) {
                                // Gaussian intensity weights:
                                w = exp(-(
                                        (tmp_im[ I(x, y, z, a_dimx, a_dimy) ] - tmp_im [ I(i, j, k, a_dimx, a_dimy) ])*
                                        (tmp_im[ I(x, y, z, a_dimx, a_dimy) ] - tmp_im [ I(i, j, k, a_dimx, a_dimy) ]) /
                                        (2.0 * sigma_r * sigma_r))
                                        - ((x - i)*(x - i) + (y - j)*(y - j) + (z - k)*(z - k)) / (2.0 * sigma_d * sigma_d)
                                        );

                                // Bilateral filter response:
                                sum_f  += w;
                                sum_fi += w * tmp_im[ I(x, y, z, a_dimx, a_dimy) ];      
                            }

                    // Set out voxel:                    
                    if ((sum_fi / sum_f) < 0)
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = 0;
                    else if (tmp > UCHAR_MAX)
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = UCHAR_MAX;
                    else
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = (unsigned char) (sum_fi / sum_f);
                }

        // Prepare for next iteration:
        memcpy(tmp_im, tmp_im2, a_dimx * a_dimy * a_dimz * sizeof (unsigned char));
    }

    // Crop image:
    P3D_TRY(p3dCrop3D_8(tmp_im2, out_im, a_dimx, a_dimy, a_dimz, a_rad, NULL, NULL));


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Bilateral filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:	
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    // Return error:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dBilateralFilter3D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int size,
        const double sigma_d,
        const double sigma_r,
        const int iter,
        int (*wr_log)(const char*, ...),
        int (*wr_progress)(const int, ...)
        ) {
    // Temporary padded input-volume:
    unsigned short* tmp_im;
    unsigned short* tmp_im2;

    int a_dimx, a_dimy, a_dimz;
    int i, j, k;
    int x, y, z;
    int ct;

    // Variables for computing gaussian kernel:
    int a_rad;
    double tmp;

    // Variables for filter management:
    double w, sum_f, sum_fi;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dBilateralFilter3D_16");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Applying bilateral filter...");
        wr_log("\tKernel size: %d.", size);
        wr_log("\tDomain sigma: %0.3f.", sigma_d);
        wr_log("\tRange sigma: %0.3f.", sigma_r);
        wr_log("\tNumber of iterations: %d.", iter);
    }


    // Init variables:
    a_rad = size / 2; // integer division   

    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad * 2;
    a_dimy = dimy + a_rad * 2;
    a_dimz = dimz + a_rad * 2;

    // Try to allocate memory:
    P3D_TRY(tmp_im = (unsigned short*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned short)));
    P3D_TRY(tmp_im2 = (unsigned short*) malloc(a_dimx * a_dimy * a_dimz * sizeof (unsigned short)));
    P3D_TRY(p3dReplicatePadding3D_16(in_im, tmp_im, dimx, dimy, dimz, a_rad, NULL, NULL));

    for (ct = 0; ct < iter; ct++) {
        // Volume scanning:
#pragma omp parallel for private(i, j, x, y, z, sum_f, sum_fi, w, tmp)
        for (k = a_rad; k < (a_dimz - a_rad); k++)
            for (j = a_rad; j < (a_dimy - a_rad); j++)
                for (i = a_rad; i < (a_dimx - a_rad); i++) {
                    // Init variables:
                    sum_f = 0.0;
                    sum_fi = 0.0;

                    // Convolve (i,j,k) voxel:
                    for (z = (k - a_rad); z <= (k + a_rad); z++)
                        for (y = (j - a_rad); y <= (j + a_rad); y++)
                            for (x = (i - a_rad); x <= (i + a_rad); x++) {
                                // Gaussian intensity weights:
                                w = exp(-(
                                        (tmp_im[ I(x, y, z, a_dimx, a_dimy) ] - tmp_im [ I(i, j, k, a_dimx, a_dimy) ])*
                                        (tmp_im[ I(x, y, z, a_dimx, a_dimy) ] - tmp_im [ I(i, j, k, a_dimx, a_dimy) ]) /
                                        (2.0 * sigma_r * sigma_r))
                                        - ((x - i)*(x - i) + (y - j)*(y - j) + (z - k)*(z - k)) / (2.0 * sigma_d * sigma_d)
                                        );

                                // Bilateral filter response:
                                sum_f  += w;
                                sum_fi += w * tmp_im[ I(x, y, z, a_dimx, a_dimy) ];      
                            }

                    // Set out voxel:                    
                    if ((sum_fi / sum_f) < 0)
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = 0;
                    else if (tmp > USHRT_MAX)
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = USHRT_MAX;
                    else
                        tmp_im2[ I(i, j, k, a_dimx, a_dimy) ] = (unsigned short) (sum_fi / sum_f);
                }

        // Prepare for next iteration:
        memcpy(tmp_im, tmp_im2, a_dimx * a_dimy * a_dimz * sizeof (unsigned short));
    }

    // Crop image:
    P3D_TRY(p3dCrop3D_16(tmp_im2, out_im, a_dimx, a_dimy, a_dimz, a_rad, NULL, NULL));


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Bilateral filter applied successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    // Return success:
    return P3D_SUCCESS;


MEM_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }

    // Release resources:	
    if (tmp_im != NULL) free(tmp_im);
    if (tmp_im2 != NULL) free(tmp_im2);

    // Return error:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

