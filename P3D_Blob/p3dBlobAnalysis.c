#include <omp.h>


#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include <stdio.h>

#include "p3dBlob.h"
#include "p3dTime.h"

#include "Common/p3dBoundingBoxList.h"
#include "Common/p3dConnectedComponentsLabeling.h"
#include "Common/p3dUtils.h"

void _drawline(unsigned char* in_im,
        unsigned int* lbl_im,
        unsigned int curr_lbl,
        unsigned char* out_im,
        const double theta, const double phi,
        const int cen_i, const int cen_j, const int cen_k,
        const int dimx, const int dimy, const int dimz,
        unsigned char lbl) {
    // Indexes for volume scanning:
    double x, y, z;
    double prec_x, prec_y, prec_z;

    double end_x1, end_y1, end_z1;
    double end_x2, end_y2, end_z2;


    unsigned char prec_status, curr_status;

    // Test line direction vectors:
    double r_0, r_1, r_2;
    double normr_0, normr_1, normr_2;

    /*fprintf(stderr, "rot_theta[ct_rot]: %0.5f\n", rot_theta[ct_rot]);
    fprintf(stderr, "rot_phi[ct_rot]: %0.5f\n", rot_phi[ct_rot]);*/

    //   Define direction vector (versor):
    r_0 = cos(phi) * sin(theta);
    r_1 = sin(phi) * sin(theta);
    r_2 = cos(theta);

    // Normalize direction vector to +1 for the maximum component,
    // doing so we can save iteration on next cycle:
    if ((fabs(r_0) >= fabs(r_1)) && (fabs(r_0) >= fabs(r_2))) {
        normr_0 = 1.0;
        normr_1 = r_1 / r_0;
        normr_2 = r_2 / r_0;
    } else if ((fabs(r_1) >= fabs(r_0)) && (fabs(r_1) >= fabs(r_2))) {
        normr_1 = 1.0;
        normr_0 = r_0 / r_1;
        normr_2 = r_2 / r_1;
    } else {
        normr_2 = 1.0;
        normr_0 = r_0 / r_2;
        normr_1 = r_1 / r_2;
    }

    // Make sure signs are correct:		
    normr_0 = (r_0 > 0) ? fabs(normr_0) : (-fabs(normr_0));
    normr_1 = (r_1 > 0) ? fabs(normr_1) : (-fabs(normr_1));
    normr_2 = (r_2 > 0) ? fabs(normr_2) : (-fabs(normr_2));


    x = cen_i;
    y = cen_j;
    z = cen_k;

    // Reset status:
    prec_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];
    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

    // Init variables:
    prec_x = x;
    prec_y = y;
    prec_z = z;

    // Explore the incremental and decremental sides while edges of 
    // VOI are reached:

    // Explore the incremental and decremental sides while edges of 
    // mask are reached of image edges are reached:

    // Explore the incremental versus while edges of VOI are reached:
    while ((x >= 0) && (x < dimx) &&
            (y >= 0) && (y < dimy) &&
            (z >= 0) && (z < dimz) &&
            (lbl_im[ I((int) x, (int) y, (int) z, dimx, dimy) ] == curr_lbl)
            ) {
        out_im[ I((int) x, (int) y, (int) z, dimx, dimy)] = lbl;
        curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

        // If we reach the background save coords:
        if (curr_status != prec_status) {
            prec_status = curr_status;
        }

        prec_x = x;
        prec_y = y;
        prec_z = z;

        x = x + normr_0;
        y = y + normr_1;
        z = z + normr_2;
    }

    // Get end point of the test line:
    end_x1 = prec_x;
    end_y1 = prec_y;
    end_z1 = prec_z;

    // Reset "center" of the line:
    x = cen_i;
    y = cen_j;
    z = cen_k;

    // Reset status:
    prec_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];
    curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

    // Init variables:
    prec_x = x;
    prec_y = y;
    prec_z = z;

    // Explore the decremental versus while edges of VOI are reached:
    while ((x >= 0) && (x < dimx) &&
            (y >= 0) && (y < dimy) &&
            (z >= 0) && (z < dimz) &&
            (lbl_im[ I((int) x, (int) y, (int) z, dimx, dimy) ] == curr_lbl)
            ) {

        out_im[ I((int) x, (int) y, (int) z, dimx, dimy)] = lbl;
        curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

        // If we reach the object save coords:
        if (curr_status != prec_status) {
            prec_status = curr_status;
        }

        prec_x = x;
        prec_y = y;
        prec_z = z;

        x = x - normr_0;
        y = y - normr_1;
        z = z - normr_2;
    }

    // Get the other end point of the test line:
    end_x2 = prec_x;
    end_y2 = prec_y;
    end_z2 = prec_z;

}

int p3dBlobAnalysis(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        BlobStats* out_stats, // OUT: Statistics
        unsigned char* blob_im, // OUT: Balls image
        unsigned char* star_im, // OUT: Balls image
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize, // IN: spatial resolution
        const int conn,
        const int max_rot,
        const int skip_borders,
        int (*wr_log)(const char*, ...)
        ) {	

    unsigned int* dt_im;
    unsigned int* lbl_im;
    int i, j, k;
    int a, b, c;
    int dist_x, dist_y, dist_z;
    int delta;
	unsigned int ct;
    int rad;
    int max_i, max_j, max_k;
    double cen_i, cen_j, cen_k, cen_ct;
    bb_t curr_bb;

    double* v_length = NULL;
    double* rot_theta = NULL;
    double* rot_phi = NULL;

    int free_flag = P3D_FALSE;
    unsigned int num_el;
    unsigned int curr_lbl;
    unsigned int max_sph = 0;
    int ct_rot;
    unsigned int* volumes;
    bb_t* bbs;

    double mean, mean_sq;
    double bb_vol;
    double l_min, l_max;
    int min_idx, max_idx;

    // Indexes for volume scanning:
    double x, y, z;
    double prec_x, prec_y, prec_z;

    double end_x1, end_y1, end_z1;
    double end_x2, end_y2, end_z2;


    // Test line direction vectors:
    double r_0, r_1, r_2;
    double normr_0, normr_1, normr_2;

    // Variables for length management:
    double length;


    unsigned char prec_status, curr_status;
    unsigned int iseed;
	time_t t;
	
    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing blob analysis...");
        wr_log("\tAdopted voxelsize: %0.6f mm.", voxelsize);
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
        if (skip_borders == P3D_TRUE)
            wr_log("\tBorders skipped. ");
        wr_log("\tNumber of random orientations: %d.", max_rot);
    }

    // Get distance transform and allocate memory for related image:
    P3D_MEM_TRY(dt_im = (unsigned int*) malloc(dimx * dimy * dimz * sizeof (unsigned int)));
    p3dSquaredEuclideanDT(in_im, dt_im, dimx, dimy, dimz, NULL);

	if (wr_log != NULL) {
		wr_log("\t----");
        wr_log("\tDistance transform successfully computed.");
	}

    // Get connected components labeled image and allocate memory for related image:
    P3D_MEM_TRY(lbl_im = (unsigned int*) malloc(dimx * dimy * dimz * sizeof (unsigned int)));
	P3D_TRY(p3dConnectedComponentsLabeling_uint(in_im, lbl_im, &num_el, &volumes, &bbs,
            dimx, dimy, dimz, conn, P3D_FALSE, skip_borders));	

    if (blob_im == NULL) {
        free_flag = P3D_TRUE;
        P3D_MEM_TRY(blob_im = (unsigned char*) calloc(dimx * dimy * dimz, sizeof (unsigned char)));
    }

	if (wr_log != NULL) {
        wr_log("\tBlob labeling successfully performed.");
	}
	

    // Allocate memory for output stats:
    out_stats->blobCount = num_el;

    if (num_el > 0) {
        out_stats->volume = (double*) malloc(num_el * sizeof (double));
        out_stats->aspect_ratio = (double*) malloc(num_el * sizeof (double));
        out_stats->extent = (double*) malloc(num_el * sizeof (double));
        out_stats->eq_sph = (double*) malloc(num_el * sizeof (double));
        out_stats->max_sph = (double*) malloc(num_el * sizeof (double));
        out_stats->l_min = (double*) malloc(num_el * sizeof (double));
        out_stats->l_max = (double*) malloc(num_el * sizeof (double));
        out_stats->sphericity = (double*) malloc(num_el * sizeof (double));

        rot_theta = (double*) malloc(max_rot * sizeof (double));
        rot_phi = (double*) malloc(max_rot * sizeof (double));
        v_length = (double*) malloc(max_rot * sizeof (double));

    } else {
        out_stats->volume = NULL;
        out_stats->aspect_ratio = NULL;
        out_stats->extent = NULL;
        out_stats->eq_sph = NULL;
        out_stats->max_sph = NULL;
        out_stats->l_min = NULL;
        out_stats->l_max = NULL;
        out_stats->sphericity = NULL;
    }



    // Scanning volume determining values for further use in stats computing:
    //#pragma omp parallel for private ( curr_bb, curr_lbl, max_sph, i, j, k, a, b, c, dist_x, dist_y, dist_z, dist_min, dist_max, max_i, max_j, max_k, rad, delta )
    for (ct = (unsigned int) 0; ct < num_el; ct = (unsigned int) (ct+1)) {
        // Get bounding box:
        curr_bb.min_x = bbs[ct].min_x;
        curr_bb.max_x = bbs[ct].max_x;
        curr_bb.min_y = bbs[ct].min_y;
        curr_bb.max_y = bbs[ct].max_y;
        curr_bb.min_z = bbs[ct].min_z;
        curr_bb.max_z = bbs[ct].max_z;

        // Compute distances in three directions:
        dist_x = curr_bb.max_x - curr_bb.min_x + 1;
        dist_y = curr_bb.max_y - curr_bb.min_y + 1;
        dist_z = curr_bb.max_z - curr_bb.min_z + 1;


        // Set current label:
        curr_lbl = (unsigned int) (ct + 3);

        // Init counters:
        max_sph = 0;

        ///
        /// Apply the "star" algorithm in order to get the minimum 
        /// and maximum axis length:
        ///

        // Scan bounding box to first get the baricenter:
        cen_i = 0.0;
        cen_j = 0.0;
        cen_k = 0.0;
        cen_ct = 0.0;
        for (k = curr_bb.min_z; k <= curr_bb.max_z; k++)
            for (j = curr_bb.min_y; j <= curr_bb.max_y; j++)
                for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                    if (lbl_im[ I(i, j, k, dimx, dimy) ] == curr_lbl) {
                        // Compute the center of mass:
                        cen_i = cen_i + (double) i;
                        cen_j = cen_j + (double) j;
                        cen_k = cen_k + (double) k;
                        cen_ct = cen_ct + 1.0;
                    }
                }

        cen_i = cen_i / cen_ct;
        cen_j = cen_j / cen_ct;
        cen_k = cen_k / cen_ct;

        // ---------------------------------------------------------------

		// Now apply the "star" algorithm:
		iseed = (unsigned int) time(&t);
		srand(iseed);

        // For each of the max_rot direction computed:
        for (ct_rot = 0; ct_rot < max_rot; ct_rot++) 
		{
            rot_theta[ct_rot] = (rand() / ((double) RAND_MAX + 1)) * M_PI;
            rot_phi[ct_rot] = (rand() / ((double) RAND_MAX + 1)) * M_PI;

            /*fprintf(stderr, "rot_theta[ct_rot]: %0.5f\n", rot_theta[ct_rot]);
            fprintf(stderr, "rot_phi[ct_rot]: %0.5f\n", rot_phi[ct_rot]);*/

            //   Define direction vector (versor):
            r_0 = cos(rot_phi[ct_rot]) * sin(rot_theta[ct_rot]);
            r_1 = sin(rot_phi[ct_rot]) * sin(rot_theta[ct_rot]);
            r_2 = cos(rot_theta[ct_rot]);

            // Normalize direction vector to +1 for the maximum component,
            // doing so we can save iteration on next cycle:
            if ((fabs(r_0) >= fabs(r_1)) && (fabs(r_0) >= fabs(r_2))) {
                normr_0 = 1.0;
                normr_1 = r_1 / r_0;
                normr_2 = r_2 / r_0;
            } else if ((fabs(r_1) >= fabs(r_0)) && (fabs(r_1) >= fabs(r_2))) {
                normr_1 = 1.0;
                normr_0 = r_0 / r_1;
                normr_2 = r_2 / r_1;
            } else {
                normr_2 = 1.0;
                normr_0 = r_0 / r_2;
                normr_1 = r_1 / r_2;
            }

            // Make sure signs are correct:		
            normr_0 = (r_0 > 0) ? fabs(normr_0) : (-fabs(normr_0));
            normr_1 = (r_1 > 0) ? fabs(normr_1) : (-fabs(normr_1));
            normr_2 = (r_2 > 0) ? fabs(normr_2) : (-fabs(normr_2));

            x = cen_i;
            y = cen_j;
            z = cen_k;
			

            // Reset status:
            prec_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];
            curr_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];

            // Init variables:
            prec_x = x;
            prec_y = y;
            prec_z = z;

            // Explore the incremental and decremental sides while edges of 
            // VOI are reached:

            // Explore the incremental and decremental sides while edges of 
            // mask are reached of image edges are reached:

            // Explore the incremental versus while edges of VOI are reached:
            while ( (((int) (x + 0.5)) >= 0) && (((int) (x + 0.5)) < dimx) &&
                    (((int) (y + 0.5)) >= 0) && (((int) (y + 0.5)) < dimy) &&
                    (((int) (z + 0.5)) >= 0) && (((int) (z + 0.5)) < dimz) &&
                    (lbl_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ] == curr_lbl)
                    ) {
                curr_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];

                // If we reach the background save coords:
                if (curr_status != prec_status) {
                    prec_status = curr_status;
                }

                prec_x = x;
                prec_y = y;
                prec_z = z;

                x = x + normr_0;
                y = y + normr_1;
                z = z + normr_2;
            }

            // Get end point of the test line:
            end_x1 = prec_x;
            end_y1 = prec_y;
            end_z1 = prec_z;

            // Reset "center" of the line:
            x = cen_i;
            y = cen_j;
            z = cen_k;

            // Reset status:
            prec_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];
            curr_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];

            // Init variables:
            prec_x = x;
            prec_y = y;
            prec_z = z;

            // Explore the decremental versus while edges of VOI are reached:
            while ( (((int) (x + 0.5)) >= 0) && (((int) (x + 0.5)) < dimx) &&
                    (((int) (y + 0.5)) >= 0) && (((int) (y + 0.5)) < dimy) &&
                    (((int) (z + 0.5)) >= 0) && (((int) (z + 0.5)) < dimz) &&
                    (lbl_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ] == curr_lbl)
                    ) {
                curr_status = in_im[ I((int) (x + 0.5), (int) (y + 0.5), (int) (z + 0.5), dimx, dimy) ];

                // If we reach the object save coords:
                if (curr_status != prec_status) {
                    prec_status = curr_status;
                }

                prec_x = x;
                prec_y = y;
                prec_z = z;

                x = x - normr_0;
                y = y - normr_1;
                z = z - normr_2;
            }

            // Get the other end point of the test line:
            end_x2 = prec_x;
            end_y2 = prec_y;
            end_z2 = prec_z;

            //}

            // Get distance:
            length = (end_x1 - end_x2)*(end_x1 - end_x2) +
                    (end_y1 - end_y2)*(end_y1 - end_y2) +
                    (end_z1 - end_z2)*(end_z1 - end_z2);
            length = sqrt(length) + 1.0;

            // Add current length to array for this grid:
            v_length[ct_rot] = length*voxelsize;
        }

        // ------------------------------------------------------------

        // Extract minimum and maximum length from the array 
        // of all the considered lengths:
        l_max = 0.0;
        max_idx = 0;
        l_min = (double) (UINT_MAX);
        min_idx = 0;
        for (i = 0; i < max_rot; i++) {
            if (v_length[i] > l_max) {
                l_max = v_length[i];
                max_idx = i;
            }
            if (v_length[i] < l_min) {
                l_min = v_length[i];
                min_idx = i;
            }
        }

        out_stats->l_max[ct] = l_max;
        out_stats->l_min[ct] = l_min;

        // Draw the axis lines (if required):
        if (star_im != NULL) {
            // Draw minor axis with label 2:
            _drawline(in_im, lbl_im, curr_lbl, star_im, rot_theta[min_idx], rot_phi[min_idx],
                    (int) cen_i, (int) cen_j, (int) cen_k, dimx, dimy, dimz, 2);
            // Draw major axis with label 3:
            _drawline(in_im, lbl_im, curr_lbl, star_im, rot_theta[max_idx], rot_phi[max_idx],
                    (int) cen_i, (int) cen_j, (int) cen_k, dimx, dimy, dimz, 3);
            // Draw center of mass with label 1:
            star_im[ I( (int) cen_i, (int) cen_j, (int) cen_k, dimx, dimy) ] = 1;
        }

        ///
        /// Maximal sphere part:
        ///

        // Scan bounding box:
        for (k = curr_bb.min_z; k <= curr_bb.max_z; k++) {
            for (j = curr_bb.min_y; j <= curr_bb.max_y; j++) {
                for (i = curr_bb.min_x; i <= curr_bb.max_x; i++) {
                    if (lbl_im[ I(i, j, k, dimx, dimy) ] == curr_lbl) {
                        // Compute maximum inscribed sphere:
                        if (dt_im[ I(i, j, k, dimx, dimy) ] > max_sph) {
                            max_sph = dt_im[ I(i, j, k, dimx, dimy) ];

                            max_i = i;
                            max_j = j;
                            max_k = k;
                        }
                    }
                }
            }
        }

        // Get the radius:
        rad = (int) (sqrt((double) max_sph));

        // Fill the ball (if required as output):
        if (free_flag == P3D_FALSE) {
            for (c = (max_k - rad); c <= (max_k + rad); c++)
                for (b = (max_j - rad); b <= (max_j + rad); b++)
                    for (a = (max_i - rad); a <= (max_i + rad); a++) {
                        // We are scanning the bounding box, so we need to be sure
                        // if current position (a,b,c) is inside the ball:
                        delta = ((a - max_i)*(a - max_i) + (b - max_j)*(b - max_j)
                                + (c - max_k)*(c - max_k));

                        if ((a >= 0) && (b >= 0) && (c >= 0) &&
                                (a < dimx) && (b < dimy) && (c < dimz) &&
                                (delta <= (rad * rad))) {
                            blob_im [ I(a, b, c, dimx, dimy) ] = OBJECT;
                        }
                    }
        }

        //printf("ct: %d - i: %d - j: %d - k: %d - rad: %d\n",ct,max_i,max_j,max_k,rad);

        // Set ouput parameters:


        ///
        /// Volume:		
        ///
        out_stats->volume[ct] = volumes[ct]*(voxelsize * voxelsize * voxelsize);

        // Aspect ratio (avoiding division by zero):
        if (out_stats->l_max[ct] != 0) {
            out_stats->aspect_ratio[ct] = out_stats->l_min[ct] / ((double) out_stats->l_max[ct]);
        } else {
            out_stats->aspect_ratio[ct] = 0;
        }

        ///
        /// Diameter of the maximum inscribed sphere:		
        ///
        if (max_sph != 0) {
            out_stats->max_sph[ct] = 2.0 * sqrt((double) max_sph) * voxelsize;
        } else {
            out_stats->max_sph[ct] = 0.0;
        }

        ///
        /// Equivalent diameter, i.e. the diameter of a sphere with
        /// the same volume as the blob, computed exploiting the
        /// inverse formula of the volume of a sphere:
        ///
        out_stats->eq_sph[ct] = 2 * pow((3.0 * out_stats->volume[ct] ) / (4.0 * M_PI), 1 / 3.0);

        ///
        /// Sphericity:
        ///
        out_stats->sphericity[ct] = out_stats->max_sph[ct] / out_stats->eq_sph[ct];

        ///
        /// Extent parameter, i.e. the ratio between the volume of the
        /// blob and the volume of the minimum bounding box (the 
        /// smallest parallelepiped oriented according to image axis 
        /// containing the blob).
        ///
        bb_vol = ((double) dist_x) * ((double) dist_y) * ((double) dist_z);
        if (bb_vol > 0)
            out_stats->extent[ct] = volumes[ct] / bb_vol;
        else
            out_stats->extent[ct] = 0.0;
            

    }

    // Print out number of connected components and mean values of parameters:
    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tNumber of identified blobs: %d. ", out_stats->blobCount);

        if (num_el > 0) {
            // Compute mean +/- std values of volume:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->volume[ct]);
                mean_sq += (double) (out_stats->volume[ct] * out_stats->volume[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tVolumes: %0.3f +/- %0.3f [mm^3].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of l_min:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->l_min[ct]);
                mean_sq += (double) (out_stats->l_min[ct] * out_stats->l_min[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tMinor axis length: %0.3f +/- %0.3f [mm].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of l_max:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->l_max[ct]);
                mean_sq += (double) (out_stats->l_max[ct] * out_stats->l_max[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tMajor axis length: %0.3f +/- %0.3f [mm].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of aspect ratio:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->aspect_ratio[ct]);
                mean_sq += (double) (out_stats->aspect_ratio[ct] * out_stats->aspect_ratio[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tAspect ratio: %0.3f +/- %0.3f [-].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of d_max:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->max_sph[ct]);
                mean_sq += (double) (out_stats->max_sph[ct] * out_stats->max_sph[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tDiameter of maximal inscribed sphere: %0.3f +/- %0.3f [mm].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of d_eq:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->eq_sph[ct]);
                mean_sq += (double) (out_stats->eq_sph[ct] * out_stats->eq_sph[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tDiameter of equivalent sphere: %0.3f +/- %0.3f [mm].", mean, mean_sq);

            // Compute mean +/- std values of sphericity:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->sphericity[ct]);
                mean_sq += (double) (out_stats->sphericity[ct] * out_stats->sphericity[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tSphericity: %0.3f +/- %0.3f [-].", mean, sqrt(mean_sq));

            // Compute mean +/- std values of extent:
            mean = 0.0;
            mean_sq = 0.0;
            for (ct = (unsigned int) 0; ct < out_stats->blobCount; ct = (unsigned int) (ct + 1)) {
                mean += (double) (out_stats->extent[ct]);
                mean_sq += (double) (out_stats->extent[ct] * out_stats->extent[ct]);
            }
            mean = mean / ((double) (out_stats->blobCount));
            mean_sq = mean_sq / ((double) (out_stats->blobCount)) - mean*mean;

            wr_log("\tExtent: %0.3f +/- %0.3f [-].", mean, sqrt(mean_sq));
        }
    }


    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Blob analysis performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Release resources:   
    if (free_flag == P3D_TRUE) free(blob_im);
    if (dt_im != NULL) free(dt_im);
    if (lbl_im != NULL) free(lbl_im);
    if (volumes != NULL) free(volumes);
    if (bbs != NULL) free(bbs);

	if (rot_theta != NULL) free(rot_theta);
    if (rot_phi != NULL) free(rot_phi);
    if (v_length != NULL) free(v_length);

    // Return OK:
    return P3D_SUCCESS;

MEM_ERROR:

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory or maximum number of blobs reached. Program will exit.");
    }

    // Release resources:    
    if (free_flag == P3D_TRUE) free(blob_im);
    if (dt_im != NULL) free(dt_im);
    if (lbl_im != NULL) free(lbl_im);
    if (volumes != NULL) free(volumes);
    if (bbs != NULL) free(bbs);

    if (rot_theta != NULL) free(rot_theta);
    if (rot_phi != NULL) free(rot_phi);
    if (v_length != NULL) free(v_length);

    // Return OK:
    return P3D_ERROR;

}
