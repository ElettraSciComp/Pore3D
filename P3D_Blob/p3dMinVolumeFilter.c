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

/*int p3dMinVolumeFilter2D(
        unsigned char* in_rev,
        unsigned char* out_rev,
        const unsigned int dimx,
        const unsigned int dimy,
        const unsigned int min_area,
        int conn,
        int (*wr_log)(const char*, ...)
        ) {
    unsigned short* lbl_rev;

    unsigned int* volumes;
    unsigned short num_el, lbl;

    unsigned int i, j, ct;

    char* auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dMinVolumeFilter3D");
    if (strcmp(auth_code, "1") != 0) goto AUTH_ERROR;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Removing blobs having area below %d pixels...", min_area);
        if (conn == CONN4)
            wr_log("\t4-connectivity used. ");
        else // default:
            wr_log("\t8-connectivity used. ");
    }

    // Initialize output cloning input:
    memcpy(out_rev, in_rev, dimx * dimy * sizeof (unsigned char));

    // Allocate memory for labels REV:
    lbl_rev = (unsigned short*) malloc(dimx * dimy * sizeof (unsigned short));


    // Perform connected component labeling:
    p3dConnectedComponentsLabeling2D(in_rev, lbl_rev, &volumes, &num_el, dimx,
            dimy, conn, NULL);

    // Init counter:
    ct = 0;

    // For each connected component labeled:
    for (lbl = 0; lbl < num_el; lbl++) {
        // If volume of connected component is below minimum volume:
        if (volumes[lbl] < min_area) {
            // Remove connected component:
            for (j = 0; j < dimy; j++)
                for (i = 0; i < dimx; i++) {
                    // Labels start from 3:
                    if (lbl_rev[ I(i, j, dimx) ] == (lbl + 3)) {
                        out_rev[ I(i, j, dimx) ] = 0;
                    }
                }

            // Keep trace of connected components removed:
            ct++;
        }
    }



    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tPore3D - %d of %f blobs removed.", ct, num_el);
        wr_log("Pore3D - Removal of blobs having area below %d pixels performed successfully in %s.", min_area, p3dGetElapsedTime_s());
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

AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;
}*/

int p3dMinVolumeFilter3D(
        unsigned char* in_im,
        unsigned char* out_im,
        const int dimx,
        const int dimy,
        const int dimz,
        const int min_volume,
        int conn,
        int (*wr_log)(const char*, ...)
        ) {
    unsigned short* lbl_im;

    unsigned int* volumes;
    bb_t* bbs;

    bb_t curr_bb;
    unsigned short num_el, curr_lbl;

    int i, j, k, ccl_removed_ct;

    int ct, err_code;



    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dMinVolumeFilter3D");
    if (auth_code == '0') goto AUTH_ERROR;*/

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Removing blobs having volume below %d voxels...", min_volume);
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
    }

    // Initialize output cloning input:
    memcpy(out_im, in_im, dimx * dimy * dimz * sizeof (unsigned char));


    // Allocate memory for labels:
    P3D_TRY(lbl_im = (unsigned short*) malloc(dimx * dimy * dimz * sizeof (unsigned short)));

    // Perform connected component labeling:
    P3D_TRY(p3dConnectedComponentsLabeling_ushort(in_im, lbl_im, &num_el, &volumes, &bbs, dimx,
            dimy, dimz, conn, P3D_FALSE, P3D_FALSE));


    // Init counter:
    ccl_removed_ct = 0;

    // For each connected component labeled:
#pragma omp parallel for private(i, j, k, curr_lbl, curr_bb)
    for (ct = 0; ct < num_el; ct++) {
        // If volume of connected component is below minimum volume:
        if (volumes[ct] < min_volume) {
            // Get bounding box:
            curr_bb.min_x = bbs[ct].min_x;
            curr_bb.max_x = bbs[ct].max_x;
            curr_bb.min_y = bbs[ct].min_y;
            curr_bb.max_y = bbs[ct].max_y;
            curr_bb.min_z = bbs[ct].min_z;
            curr_bb.max_z = bbs[ct].max_z;

            // Set current label:
            curr_lbl = ct + 3;

            // Scan bounding box:
            for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
                for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                    for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                        if (lbl_im[ I(i, j, k, dimx, dimy) ] == curr_lbl) {
                            out_im[ I(i, j, k, dimx, dimy) ] = 0;
                        }
                    }

            // Keep trace of connected components removed:
            ccl_removed_ct++;
        }
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("\tPore3D - %d of %d blobs removed.", ccl_removed_ct, num_el);
        wr_log("Pore3D - Removal of blobs having volume below %d voxels performed successfully in %dm%0.3fs.", min_volume, p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Free memory:
    if (lbl_im != NULL) free(lbl_im);
    if (volumes != NULL) free(volumes);
    if (bbs != NULL) free(bbs);

    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Free memory:
    if (lbl_im != NULL) free(lbl_im);
    if (volumes != NULL) free(volumes);
    if (bbs != NULL) free(bbs);


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














