#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>

#include "p3dFilt.h"
#include "p3dTime.h"
#include "p3dAuth.h"
#include "Common/p3dRingRemoverCommon.h"

#define EPS 0.0001

#define P3DBOINHAIBELRINGREMOVER_ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
#define P3DBOINHAIBELRINGREMOVER_MEDIAN(a,n) _p3dBoinHaibelRingRemover_kth_smallest(a,n,((n)/2))

double _p3dBoinHaibelRingRemover_kth_smallest(double* a, int n, int k) {
    register int i, j, l, m;
    register double x;

    l = 0;
    m = n - 1;
    while (l < m) {
        x = a[k];
        i = l;
        j = m;
        do {
            while (a[i] < x) i++;
            while (x < a[j]) j--;
            if (i <= j) {
                P3DBOINHAIBELRINGREMOVER_ELEM_SWAP(a[i], a[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (j < k) l = i;
        if (k < i) m = j;
    }
    return a[k];
}


// This procedure removes ring artifacts from CT images.

int p3dBoinHaibelRingRemover2D_8(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int centerX,
        const int centerY,
        const int winsize, // IN: width of the moving average
        const int iterations, // IN: filter can be re-iterated
        const double precision,
        int (*wr_log)(const char*, ...)
        ) {
    int i, j; // generic counters
    int k, tmp_k;
    int sum;

    double tmp_val;

    // Polar images:
    unsigned char* p_in_im; // no need for malloc
    unsigned char* c_out_im; // no need for malloc

    int p_dim, ct, it_ct;

    // Artifacts vectors:
    double* glob_art;
    double* v;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dBoinHaibelRingRemover2D_8");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // STEP2: Transform in polar coordinates. Remember that memory allocation
    // for output images is performed within the procedure.

    p3dCartesian2polar_8(in_im, &p_in_im, dimx, dimy, centerX, centerY, precision, &p_dim);

    // STEP3: Artifact template selection.

    // Allocate dynamic memory:

    // Now that we know the dimensions of polar image we can allocate
    // the global artifacts correction vector:
    glob_art = (double*) calloc(p_dim, sizeof (double));

    // Allocate memory for the column:
    // v = (double*) malloc(row_ct * sizeof (double));

    for (it_ct = 0; it_ct < iterations; it_ct++) {
        // Compute median for each column and store the value in the
        // artifact vector:
#pragma omp parallel for private(k, sum)
        for (j = 0; j < p_dim; j++) {
            // Fill the column array:
            sum = 0;

            for (k = 0; k < p_dim; k++) {
                sum = sum + p_in_im[ I2(j, k, p_dim) ];
            }

            // Put the sum in artifact vector:
            glob_art[j] = sum;
        }

        // Free memory for the column:
        // Cycle for each element in line to compute moving average:
#pragma omp parallel for private(v, ct, k, tmp_k)
        for (j = 0; j < p_dim; j++) {
            /*sum = 0;

            // Moving average:
            for (k = ( j - (win_size / 2)); k < ( j + (win_size / 2)); k++)
            {
                    tmp_k = k;

                    // Replicate padding:
                    if (tmp_k < 0) tmp_k = 0;
                    if (tmp_k > (p_dim - 1)) tmp_k = p_dim - 1;

                    sum = sum + glob_art[ tmp_k ];
            }

            // Moving average:
            glob_art[j] = (sum / win_size) / (glob_art[j] + EPS);*/

            // Allocate memory for the column:
            v = (double*) malloc(winsize * sizeof (double));
            ct = 0;

            for (k = (j - (winsize / 2)); k < (j + (winsize / 2)); k++) {
                tmp_k = k;

                // Replicate padding:
                if (tmp_k < 0) tmp_k = 0;
                if (tmp_k > (p_dim - 1)) tmp_k = p_dim - 1;

                v[ct++] = glob_art[ tmp_k ];
            }

            glob_art[j] = (P3DBOINHAIBELRINGREMOVER_MEDIAN(v, winsize) / (glob_art[j] + EPS));

            // Free memory for the column:
            free(v);
        }


        // The artifact template vector is multiplied for each
        // row of the polar image P:
#pragma omp parallel for private(i, tmp_val)
        for (j = 0; j < p_dim; j++)
            for (i = 0; i < p_dim; i++) {
                // Take care of intensity bounds:
                tmp_val = p_in_im[ I2(i, j, p_dim) ] * glob_art[i];

                if (tmp_val < 0.0) tmp_val = 0.0;
                if (tmp_val > UCHAR_MAX) tmp_val = UCHAR_MAX;

                p_in_im[ I2(i, j, p_dim) ] = (unsigned char) tmp_val;
            }


    }

    //p3dWriteRaw8 ( p_in_im, "C:\\p_out_im.raw", p_dim, p_dim, 1, printf );

    // Return in cartesian coordinates:
    p3dPolar2cartesian_8(p_in_im, &c_out_im, p_dim, centerX, centerY, dimx, dimy);

    // Copy C_OUT_IM to the output of the procedure:
    memcpy(out_im, c_out_im, dimx * dimy * sizeof (unsigned char));

    // Free memory:
    free(glob_art);

    free(p_in_im);
    free(c_out_im);

    return P3D_SUCCESS;
    
   /* AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

int p3dBoinHaibelRingRemover2D_16(
        unsigned short* in_im,
        unsigned short* out_im,
        const int dimx,
        const int dimy,
        const int centerX,
        const int centerY,
        const int winsize, // IN: width of the moving average
        const int iterations, // IN: filter can be re-iterated
        const double precision,
        int (*wr_log)(const char*, ...)
        ) {
    int i, j; // generic counters
    int k, tmp_k;
    int sum;

    double tmp_val;

    // Polar images:
    unsigned short* p_in_im; // no need for malloc
    unsigned short* c_out_im; // no need for malloc

    int p_dim, ct, it_ct;

    // Artifacts vectors:
    double* glob_art;
    double* v;

   /* char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dBoinHaibelRingRemover2D_16");
    if (auth_code == '0') goto AUTH_ERROR;*/


    // STEP2: Transform in polar coordinates. Remember that memory allocation
    // for output images is performed within the procedure.

    p3dCartesian2polar_16(in_im, &p_in_im, dimx, dimy, centerX, centerY, precision, &p_dim);

    // STEP3: Artifact template selection.

    // Allocate dynamic memory:

    // Now that we know the dimensions of polar image we can allocate
    // the global artifacts correction vector:
    glob_art = (double*) calloc(p_dim, sizeof (double));

    // Allocate memory for the column:
    // v = (double*) malloc(row_ct * sizeof (double));

    // Compute median for each column and store the value in the
    // artifact vector:

        for (it_ct = 0; it_ct < iterations; it_ct++) {
#pragma omp parallel for private(k, sum)
    for (j = 0; j < p_dim; j++) {
        // Fill the column array:
        sum = 0;

        for (k = 0; k < p_dim; k++) {
            sum = sum + p_in_im[ I2(j, k, p_dim) ];
        }

        // Put the sum in artifact vector:
        glob_art[j] = sum;
    }

    // Free memory for the column:
    // Cycle for each element in line to compute moving average:
#pragma omp parallel for private(v, ct, k, tmp_k)
    for (j = 0; j < p_dim; j++) {
        /*sum = 0;

        // Moving average:
        for (k = ( j - (win_size / 2)); k < ( j + (win_size / 2)); k++)
        {
                tmp_k = k;

                // Replicate padding:
                if (tmp_k < 0) tmp_k = 0;
                if (tmp_k > (p_dim - 1)) tmp_k = p_dim - 1;

                sum = sum + glob_art[ tmp_k ];
        }

        // Moving average:
        glob_art[j] = (sum / win_size) / (glob_art[j] + EPS);*/

        // Allocate memory for the column:
        v = (double*) malloc(winsize * sizeof (double));
        ct = 0;

        for (k = (j - (winsize / 2)); k < (j + (winsize / 2)); k++) {
            tmp_k = k;

            // Replicate padding:
            if (tmp_k < 0) tmp_k = 0;
            if (tmp_k > (p_dim - 1)) tmp_k = p_dim - 1;

            v[ct++] = glob_art[ tmp_k ];
        }

        glob_art[j] = (P3DBOINHAIBELRINGREMOVER_MEDIAN(v, winsize) / (glob_art[j] + EPS));

        // Free memory for the column:
        free(v);
    }


    // The artifact template vector is multiplied for each
    // row of the polar image P:
#pragma omp parallel for private(i, tmp_val)
    for (j = 0; j < p_dim; j++)
        for (i = 0; i < p_dim; i++) {
            // Take care of intensity bounds:
            tmp_val = p_in_im[ I2(i, j, p_dim) ] * glob_art[i];

            if (tmp_val < 0.0) tmp_val = 0.0;
            if (tmp_val > USHRT_MAX) tmp_val = USHRT_MAX;

            p_in_im[ I2(i, j, p_dim) ] = (unsigned short) tmp_val;
        }

        }

    //p3dWriteRaw8 ( p_in_im, "C:\\p_out_im.raw", p_dim, p_dim, 1, printf );

    // Return in cartesian coordinates:
    p3dPolar2cartesian_16(p_in_im, &c_out_im, p_dim, centerX, centerY, dimx, dimy);

    // Copy C_OUT_IM to the output of the procedure:
    memcpy(out_im, c_out_im, dimx * dimy * sizeof (unsigned short));

    // Free memory:
    free(glob_art);

    free(p_in_im);
    free(c_out_im);

    return P3D_SUCCESS;
    
/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}