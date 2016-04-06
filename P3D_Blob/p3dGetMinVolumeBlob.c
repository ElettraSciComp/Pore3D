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

int p3dGetMinVolumeBlob3D(
        unsigned char* in_rev,
        unsigned char* out_rev,
        const int dimx,
        const int dimy,
        const int dimz,
        int conn,
        int (*wr_log)(const char*, ...)
        ) {
    unsigned short* lbl_rev;

    unsigned int* volumes;
    unsigned short num_el, lbl;

    int i, j, k;

    int lbl_min;
    int vol_min;



    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dGetMinVolumeBlob3D");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Extracting blob having minimum volume...");
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
    }
    // Initialize output cloning input:
    memcpy(out_rev, in_rev, dimx * dimy * dimz * sizeof (unsigned char));

    // Allocate memory for labels:
    P3D_TRY(lbl_rev = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));

    // Perform connected component labeling:
    P3D_TRY(p3dConnectedComponentsLabeling_ushort(in_rev, lbl_rev, &num_el, &volumes, NULL, dimx,
            dimy, dimz, conn, P3D_FALSE, P3D_FALSE));


    lbl_min = 0;
    vol_min = UINT_MAX;

    // For each connected component labeled:
    for (lbl = 0; lbl < num_el; lbl++) {
        if (volumes[lbl] < vol_min) {
            lbl_min = lbl;
            vol_min = volumes[lbl];
        }
    }

    // Remove connected component:
    for (k = 0; k < dimz; k++)
        for (j = 0; j < dimy; j++)
            for (i = 0; i < dimx; i++) {
                // Labels start from 3:
                if (lbl_rev[ I(i, j, k, dimx, dimy) ] != (lbl_min + 3)) {
                    out_rev[ I(i, j, k, dimx, dimy) ] = 0;
                }
            }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tPore3D - The volume of extracted blob is %d voxels.", vol_min);
        wr_log("Pore3D - Blob having minimum area extracted successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (lbl_rev != NULL) free(lbl_rev);
    if (volumes != NULL) free(volumes);

    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Free memory:
    if (lbl_rev != NULL) free(lbl_rev);
    if (volumes != NULL) free(volumes);

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












