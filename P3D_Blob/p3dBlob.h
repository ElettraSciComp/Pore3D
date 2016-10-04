/***************************************************************************/
/* (C) 2016 Elettra - Sincrotrone Trieste S.C.p.A.. All rights reserved.   */
/*                                                                         */
/*                                                                         */
/* This file is part of Pore3D, a software library for quantitative        */
/* analysis of 3D (volume) images.                                         */
/*                                                                         */
/* Pore3D is free software: you can redistribute it and/or modify it       */
/* under the terms of the GNU General Public License as published by the   */
/* Free Software Foundation, either version 3 of the License, or (at your  */
/* option) any later version.                                              */
/*                                                                         */
/* Pore3D is distributed in the hope that it will be useful, but WITHOUT   */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   */
/* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    */
/* for more details.                                                       */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with Pore3D. If not, see <http://www.gnu.org/licenses/>.          */
/*                                                                         */
/***************************************************************************/

//
// Author: Francesco Brun
// Last modified: Sept, 28th 2016
//

#ifdef __cplusplus
extern "C" {
#endif


    /*
            Constants:
     */
#ifndef P3D_BLOB_DEFINED
#define P3D_BLOB_DEFINED

#define P3D_FALSE				-1 
#define P3D_TRUE				1 

#define P3D_AUTH_ERROR                          -1
#define P3D_MEM_ERROR			NULL	/* Leave it NULL for simplify tests */
#define P3D_SUCCESS				2		/* Any number */

#define BACKGROUND				0
#define OBJECT					UCHAR_MAX

#define CONN4   611
#define CONN8   612

    // Constants for 3D connectivity:
#define CONN6   711
#define CONN18  712
#define CONN26  713

#endif

    /*
            Macros:
     */

#ifndef P3D_BLOB_MACROS
#define P3D_BLOB_MACROS

#define I(i,j,k,N,M)    ( (j)*(N) + (i) + (k)*(N)*(M) ) 
#define MIN(x,y)        (((x) < (y))?(x):(y))
#define MAX(x,y)        (((x) > (y))?(x):(y))

    /* A sort of TRY-CATCH constructor: */
#define P3D_TRY( function ) if ( (function) == P3D_MEM_ERROR) { goto MEM_ERROR; }
#endif

#ifndef P3D_BLOB_STRUCTS_DEFINED
#define P3D_BLOB_STRUCTS_DEFINED

    // Struct for results of connected components analysis:

    typedef struct {
        double* aspect_ratio;
        unsigned int blobCount;
        double* max_sph;
        double* eq_sph;
        double* l_min;
        double* l_max;
        double* sphericity;
        double* extent;
        double* volume;
    } BlobStats;

    struct MorphometricStats {
        double BvTv;
        double BsBv;
        double TbN;
        double TbTh;
        double TbSp;
    };

    struct AnisotropyStats {
        double I; // Isotropy Index;
        double E; // Elongation Index;
    };

    struct BasicStats {
        double Vv;
        double Cv;
        double Mv;
        double Sv;
    };

    struct TextureStats {
        double FD;
    };

#endif

    int p3dBlobAnalysis(
            unsigned char* in_im,
            BlobStats* out_stats,
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
            );

    int p3dBasicAnalysis(
            unsigned char* in_im,
            struct BasicStats* out_stats,
            const int dimx,
            const int dimy,
            const int dimz,
            const double voxelsize,
            int (*wr_log)(const char*, ...)
            );

    int p3dTextureAnalysis(
            unsigned char* in_im,
            struct TextureStats* out_stats,
            const int dimx,
            const int dimy,
            const int dimz,
            int (*wr_log)(const char*, ...)
            );

    int p3dAnisotropyAnalysis(
            unsigned char* in_im,
            unsigned char* msk_im,
            struct AnisotropyStats* out_stats,
            const int dimx,
            const int dimy,
            const int dimz,
            const double voxelsize,
            const int verbose,
            int (*wr_log)(const char*, ...)
            );

    int p3dMorphometricAnalysis(
            unsigned char* in_im, // IN: Input segmented (binary) volume
            unsigned char* msk_im,
            struct MorphometricStats* out_stats, // OUT: Trabecular statistics
            const int dimx,
            const int dimy,
            const int dimz,
            const double voxelsize, // IN: voxel resolution
            int (*wr_log)(const char*, ...)
            );

    int p3dREVEstimation(
            unsigned char* in_rev, // IN: binary volume
            double** porosity, // OUT: array of porosity for the related cube side
            unsigned int** cubeEdges, // OUT: array of length of the side of the REV cube
            unsigned int* numel, // OUT: length of the arrays
            const int dimx,
            const int dimy,
            const int dimz,
            const int stepsize, // IN: stepsize for porosity
            const int centerx, // IN: X coordinate for the center of REV (-1 for 
            //     automatic determination of the center)
            const int centery, // IN: Y coordinate for the center of REV
            const int centerz, // IN: Z coordinate for the center of REV
            int (*wr_log)(const char*, ...)
            );

    int p3dChamferDT(
            unsigned char* in_im,
            unsigned short* out_im,
            const int dimx,
            const int dimy,
            const int dimz,
            const int w1,
            const int w2,
            const int w3,
            int (*wr_log)(const char*, ...)
            );

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
            );
    
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
            );


    int p3dGetMaxVolumeBlob3D(
            unsigned char* in_rev,
            unsigned char* out_rev,
            const int dimx,
            const int dimy,
            const int dimz,
            int conn,
            int (*wr_log)(const char*, ...)
            );

    int p3dGetMinVolumeBlob3D(
            unsigned char* in_rev,
            unsigned char* out_rev,
            const int dimx,
            const int dimy,
            const int dimz,
            int conn,
            int (*wr_log)(const char*, ...)
            );

    int p3dMinVolumeFilter3D(
            unsigned char* in_im,
            unsigned char* out_im,
            const int dimx,
            const int dimy,
            const int dimz,
            const int min_volume,
            int conn,
            int (*wr_log)(const char*, ...)
            );

	int p3dSquaredEuclideanDT(
        unsigned char* in_rev,
        unsigned short* out_rev,
        const int dimx,
        const int dimy,
        const int dimz,
        int (*wr_log)(const char*, ...)
        );

#ifdef __cplusplus
}
#endif