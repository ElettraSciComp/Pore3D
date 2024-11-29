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

#define GRIDPOINTS 1024
#define MAXROT 128
#define BORDER 1

/*int p3dComputeFragmentationIndex( 
    unsigned char* in_im,				// IN: Input segmented (binary) volume
        double* fri,						// OUT: Fragmentation Index
        const unsigned int dimx,
        const unsigned int dimy, 
        const unsigned int dimz,    
        const double voxelsize				// IN: voxel resolution	
        )
{
        unsigned int a_dimx, a_dimy, a_dimz;
    unsigned int a_rad;

        unsigned int i,j,k;
        unsigned int x,y,z;
    double vol1, surf1;		// To avoid oveflow problem
        double vol2, surf2;		// To avoid oveflow problem

        unsigned short* dt_im;	// Squared Euclidean DT of input REV
        unsigned char*  tmp_im; // Temporary dilated version of input image
        unsigned char*  dil_im; // Dilated version of input image


        // Init variables:
    a_rad = 1;								// 1 voxel padding and also dilation
  
        // Compute dimensions of padded REV:
        a_dimx = dimx + a_rad*2;
        a_dimy = dimy + a_rad*2;
        a_dimz = dimz + a_rad*2;

        // Try to allocate memory for distance transform:
        dt_im = (unsigned short*) malloc(dimx*dimy*dimz*sizeof(unsigned short));
        if ( dt_im == NULL ) goto ON_MEM_ERROR;

        // Try to allocate memory for temporary dilated transform:
        tmp_im = (unsigned char*) calloc(a_dimx*a_dimy*a_dimz,sizeof(unsigned char));
        if ( tmp_im == NULL ) goto ON_MEM_ERROR;

        // Try to allocate memory for dilated transform:
        dil_im = (unsigned char*) malloc(dimx*dimy*dimz*sizeof(unsigned char));
        if ( dil_im == NULL ) goto ON_MEM_ERROR;


        // Determining distance transform:
        p3dSquaredEuclideanDT(in_im, dt_im, dimx, dimy, dimz, NULL);

        // Scanning volume determining values for further use in stats computing:
    vol1 = 0.0;
    surf1 = 0.0;   
    for( k = 0; k < dimz; k++ )    
                for( j = 0; j < dimy; j++ )
                        for( i = 0; i < dimx; i++ )            
            {  
                if ( in_im[ I(i,j,k,dimx,dimy) ] == OBJECT )
                    vol1 = vol1 + 1.0;
              
                if ( dt_im[ I(i,j,k,dimx,dimy) ] == BORDER )
                    surf1 = surf1 + 1.0;
            }


        // Compute dilation:
        for( k = a_rad; k < (a_dimz - a_rad); k++ )    
                for( j = a_rad; j < (a_dimy - a_rad); j++ )
                        for( i = a_rad; i < (a_dimx - a_rad); i++ )
                        {
                                if (  in_im[ I(i - a_rad, j - a_rad, k - a_rad, dimx, dimy) ] == OBJECT )
                                {
                                        // Dilate voxel filling the neighborhood in tmp_im:					
                                        for (z = (k - a_rad); z <= (k + a_rad); z++)                    
                                                for (y = (j - a_rad); y <= (j + a_rad); y++)   
                                                        for (x = (i - a_rad); x <= (i + a_rad); x++)
                                                        {
                                                                tmp_im[ I(x,y,z,a_dimx,a_dimy) ] = OBJECT;
                                                        }
                                }
                        }


        // Crop tmp_im:
        p3dCrop(tmp_im, dil_im, a_dimx, a_dimy, a_dimz, a_rad);	
	
        // Determining distance transform of dilated image:
        p3dSquaredEuclideanDT(dil_im, dt_im, dimx, dimy, dimz, NULL);

        // Scanning volume determining values for further use:
    vol2 = 0.0;
    surf2 = 0.0;
    for( k = 0; k < dimz; k++ )    
                for( j = 0; j < dimy; j++ )
                        for( i = 0; i < dimx; i++ )            
            {  
                if ( tmp_im[ I(i,j,k,dimx,dimy) ] == OBJECT )
                    vol2 = vol2 + 1.0;

                if ( dt_im[ I(i,j,k,dimx,dimy) ] == BORDER )
                    surf2 = surf2 + 1.0;
            }


        // Compute fragmentation index:
 *fri = ( (surf1 - surf2) / (vol1 - vol2) ) / voxelsize;

	
        // Release resources:
        if ( dt_im != NULL) free(dt_im);
        if ( tmp_im != NULL) free(tmp_im);
        if ( dil_im != NULL) free(dil_im);

        // Return OK:
        return P3D_SUCCESS;

ON_MEM_ERROR:
		
        // Free resources:	
        if ( dt_im != NULL) free(dt_im);
        if ( tmp_im != NULL) free(tmp_im);
        if ( dil_im != NULL) free(dil_im);

        return P3D_ERROR;
}*/


int p3dMorphometricAnalysis(
        unsigned char* in_im, // IN: Input segmented (binary) volume
        unsigned char* msk_im,
        MorphometricStats* out_stats, // OUT: Trabecular statistics
        const int dimx,
        const int dimy,
        const int dimz,
        const double voxelsize, // IN: voxel resolution
        int (*wr_log)(const char*, ...)
        ) {

    unsigned int s, t;
	int i;

    // Indexes for volume scanning:
    double x, y, z;
    double x_prec, y_prec, z_prec;

    int end_x, end_y, end_z;
    int curr_x, curr_y, curr_z;


    // Test line direction vectors:
    double r[3], norm_r[3];
    int ct;

    // Variables for length management:
    double length;
    double bvf; // Bone Volume Fraction (i.e. BV/TV)

    long intersect_ct;
    unsigned char prec_status, curr_status;

    double tot_length;
    double tot_intersect_ct; // WARNING: Should be a long but a double
    // is used to avoid overflow problem.

    double pl;
    double rot_theta, rot_phi;

    double sv; // BS/BV [mm^2 mm^-3]
    double tb; // Tb.Th [mm]
    double tpd; // Tb.N  [mm^-1]
    double tps; // Tb.Sp [mm] 


    unsigned int iseed;
	time_t tt;


    // Variables for rotation cycle:
    int ct_rot;

    // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Performing standard morphometric analysis...");
        wr_log("\tAdopted voxelsize: %0.6f mm.", voxelsize);
    }

    // Init random seeds:
    iseed = (unsigned int) time(&tt);
    srand(iseed);

    // Get BV/TV:
    if (msk_im == NULL) {
        s = 0;
        for (i = 0; i < (dimx * dimy * dimz); i++)
            if (in_im[i] == OBJECT) s++;

        bvf = s / ((double) (dimx * dimy * dimz));
    } else {
        s = 0;
        t = 0;
        for (i = 0; i < (dimx * dimy * dimz); i++) {
            if (msk_im[i] == OBJECT) {
                t++;
                if (in_im[i] == OBJECT) s++;
            }
        }

        bvf = s / ((double) t);
    }

    // Reset counters:
    tot_length = 0.0;
    tot_intersect_ct = 0.0;

    // Explore VOI on the random points:
    for (i = 0; i < GRIDPOINTS; i++) {
        curr_x = (int) (rand() / (((double) RAND_MAX + 1) / dimx));
        curr_y = (int) (rand() / (((double) RAND_MAX + 1) / dimy));
        curr_z = (int) (rand() / (((double) RAND_MAX + 1) / dimz));

        // Check if inside MSK (if mask):
        if (msk_im != NULL) {
            if (msk_im[ I(curr_x, curr_y, curr_z, dimx, dimy)] == BACKGROUND) {
                continue;
            }

        }

        // Skip point if we aren't into the marrow:
        if (in_im[ I(curr_x, curr_y, curr_z, dimx, dimy) ] == BACKGROUND) {
            // For each of the max_rot direction computed:
            for (ct_rot = 0; ct_rot < MAXROT; ct_rot++) {
                x = curr_x;
                y = curr_y;
                z = curr_z;

                rot_theta = ((double) rand() / (double) RAND_MAX) * M_PI;
                rot_phi = ((double) rand() / (double) RAND_MAX) * M_PI;

                //   Define direction vector (versor):
                r[0] = cos(rot_phi) * sin(rot_theta);
                r[1] = sin(rot_phi) * sin(rot_theta);
                r[2] = cos(rot_theta);

                // Normalize direction vector to +1 for the maximum component,
                // doing so we can save iteration on next cycle:
                if ((fabs(r[0]) >= fabs(r[1])) && (fabs(r[0]) >= fabs(r[2]))) {
                    norm_r[0] = 1.0;
                    norm_r[1] = r[1] / r[0];
                    norm_r[2] = r[2] / r[0];
                } else if ((fabs(r[1]) >= fabs(r[0])) && (fabs(r[1]) >= fabs(r[2]))) {
                    norm_r[1] = 1.0;
                    norm_r[0] = r[0] / r[1];
                    norm_r[2] = r[2] / r[1];
                } else {
                    norm_r[2] = 1.0;
                    norm_r[0] = r[0] / r[2];
                    norm_r[1] = r[1] / r[2];
                }

                // Make sure signs are correct:
                for (ct = 0; ct < 3; ct++) {
                    if (r[ct] > 0)
                        norm_r[ct] = fabs(norm_r[ct]);
                    else
                        norm_r[ct] = (-fabs(norm_r[ct]));
                }

                x_prec = x;
                y_prec = y;
                z_prec = z;

                // Reset status:
                intersect_ct = 0;
                prec_status = 0; // false
                curr_status = 0; // false

                // Explore the incremental versus while edges of VOI are 
                // reached:
                if (msk_im == NULL) {
                    while ((x >= 0) && (x < dimx) &&
                            (y >= 0) && (y < dimy) &&
                            (z >= 0) && (z < dimz)) {
                        curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                        // If we reach the object save coords:
                        if (curr_status != prec_status) {
                            intersect_ct++;
                            prec_status = curr_status;
                        }

                        x_prec = x;
                        y_prec = y;
                        z_prec = z;

                        x = x + norm_r[0];
                        y = y + norm_r[1];
                        z = z + norm_r[2];
                    }
                } else {
                    while ((x >= 0) && (x < dimx) &&
                            (y >= 0) && (y < dimy) &&
                            (z >= 0) && (z < dimz) &&
                            (msk_im[ I((int) x, (int) y, (int) z, dimx, dimy)] == OBJECT)
                            ) {

                        curr_status = in_im[ I((int) x, (int) y, (int) z, dimx, dimy)];

                        // If we reach the object save coords:
                        if (curr_status != prec_status) {
                            intersect_ct++;
                            prec_status = curr_status;
                        }

                        x_prec = x;
                        y_prec = y;
                        z_prec = z;

                        x = x + norm_r[0];
                        y = y + norm_r[1];
                        z = z + norm_r[2];
                    }

                }

                // We complete the cycle so save end coords:
                end_x = (int) x_prec;
                end_y = (int) y_prec;
                end_z = (int) z_prec;


                // Get distance:
                length = (end_x - curr_x)*(end_x - curr_x) +
                        (end_y - curr_y)*(end_y - curr_y) +
                        (end_z - curr_z)*(end_z - curr_z);
                length = sqrt(length);

                // Add current length to array for this "porus":
                tot_length += length;
                tot_intersect_ct += (double) (intersect_ct);

            } // end loop on each rotation
        } // end if we are on a "porus"
    }


    pl = tot_intersect_ct / tot_length;

    // Prepare output results:    
    sv = 2 * pl / bvf; // BS/BV [mm^2 mm^-3] = [mm^-1]
    tb = bvf / pl; // Tb.Th [mm]
    tpd = pl; // Tb.N  [mm^-1]
    tps = (1.0 - bvf) / pl; // Tb.Sp [mm] 

    // Filling the returned struct:	
    out_stats->BvTv = bvf;
    out_stats->BsBv = sv / voxelsize;
    out_stats->TbTh = tb * voxelsize;
    out_stats->TbSp = tps * voxelsize;
    out_stats->TbN = tpd / voxelsize;

    //p3dComputeFragmentationIndex ( in_im, &out_stats->FrI, dimx, dimy, dimz, voxelsize );   

    if (wr_log != NULL) {
        wr_log("\t----");
        wr_log("\tBone Volume / Total Volume (BV/TV): %0.3f [-].", out_stats->BvTv);
        wr_log("\tBone Surface / Bone Volume (BS/BV): %0.3f [mm^-1].", out_stats->BsBv);
        wr_log("\tTrabecular Number (Tb.N): %0.3f [mm^-1].", out_stats->TbN);
        wr_log("\tTrabecular Thickness (Tb.Th): %0.3f [mm].", out_stats->TbTh);
        wr_log("\tTrabecular Separation (Tb.Sp): %0.3f [mm].", out_stats->TbSp);
    }

    // Print elapsed time (if required):
    if (wr_log != NULL) {
        wr_log("Pore3D - Morphometric analysis computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }

    // Return OK:
    return P3D_SUCCESS;
}





