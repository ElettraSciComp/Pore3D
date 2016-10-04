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
#ifndef P3D_FILT_DEFINED
#define P3D_FILT_DEFINED

#define P3D_FALSE				-1 
#define P3D_TRUE				1 

#define P3D_AUTH_ERROR                          -1
#define P3D_MEM_ERROR                           NULL	/* Leave it NULL for simplify tests */
#define P3D_IO_ERROR                            1
#define P3D_SUCCESS				2	/* Any number */

#define BACKGROUND				0
#define OBJECT					UCHAR_MAX	

    // Constants for 3D connectivity:
#define CONN6   711
#define CONN18  712
#define CONN26  713

#endif

    /*
            Macros:
     */

#ifndef P3D_MACROS
#define P3D_MACROS

#define I(i,j,k,N,M)    ( (j)*(N) + (i) + (k)*(N)*(M) )
#define I2(i,j,N)       ( (j)*(N) + (i) )

    
#define MIN(x,y)        (((x) < (y))?(x):(y))
#define MAX(x,y)        (((x) > (y))?(x):(y))

    /* A sort of TRY-CATCH constructor: */
#define P3D_TRY( function ) if ( (function) == P3D_MEM_ERROR) { goto MEM_ERROR; }
#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
#endif

    // Input - output:
    int p3dReadRaw8(char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dReadRaw16(char*, unsigned short*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dWriteRaw8(unsigned char*, char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dWriteRaw16(unsigned short*, char*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
	int p3dWriteRaw32(unsigned int*, char*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));


    // Utils:
    int p3dCrop2D_8(unsigned char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dCrop2D_16(unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dCrop3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dCrop3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dZeroPadding2D_8(unsigned char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dZeroPadding2D_16(unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dZeroPadding3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dZeroPadding3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dReplicatePadding2D_8(unsigned char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dReplicatePadding2D_16(unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dReplicatePadding3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dReplicatePadding3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dFrom16To8_batch (char*, char*, const int, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...),  int (*wr_progress)(const int, ...), int );
    
    int p3dTIFFToRaw_batch (char*, char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    
    int p3dRaw16ToTIFF_batch (char*, char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    int p3dRaw8ToTIFF_batch (char*, char*, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    
    // Basic Filters:
    int p3dAnisotropicDiffusionFilter3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, const double, const double, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dAnisotropicDiffusionFilter3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, const double, const double, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dBilateralFilter3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, const double, const double, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dBilateralFilter3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, const double, const double, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dGaussianFilter3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dGaussianFilter3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dMeanFilter2D_8(unsigned char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMeanFilter2D_16(unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMeanFilter3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMeanFilter3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dMedianFilter2D_8(unsigned char*, unsigned char*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMedianFilter2D_16(unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMedianFilter3D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dMedianFilter3D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dBoinHaibelRingRemover2D_8(unsigned char*, unsigned char*, const int, const int, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dBoinHaibelRingRemover2D_16(unsigned short*, unsigned short*, const int, const int, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));

    int p3dMunchEtAlRingRemover2D_8(unsigned char*, unsigned char*, const int, const int, const double, const double, const int, const double, const int, const double, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int);
    int p3dMunchEtAlRingRemover2D_16(unsigned short*, unsigned short*, const int, const int, const double, const double, const int, const double, const int, const double, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int);
    
    int p3dMunchEtAlRingRemover2D_8_batch(char*, char*, const int, const int, const double, const double, const int, const double, const int, const double, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    int p3dMunchEtAlRingRemover2D_16_batch(char*, char*, const int, const int, const double, const double, const int, const double, const int,  const double, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    
    int p3dSijbersPostnovRingRemover2D_8(unsigned char*, unsigned char*, const int, const int, const double, const double, const int, const double, const int, const double, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dSijbersPostnovRingRemover2D_16(unsigned short*, unsigned short*, const int, const int, const double, const double, const int, const double, const int, const double, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    
    int p3dSijbersPostnovRingRemover2D_8_batch(char*, char*, const int, const int, const double, const double, const int, const double, const int, const double, char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
    int p3dSijbersPostnovRingRemover2D_16_batch(char*, char*, const int, const int, const double, const double, const int, const double, const int,  const double, const int, char*, const int, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...), int );
   


    // Segmentation:
    //int p3dSingleRegionGrowing3D_8( unsigned char*, unsigned char*, const int, const int, const int, const int, const int, const int, const double, const int, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dKittlerThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dKittlerThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dOtsuThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dOtsuThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dPunThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dPunThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dRidlerThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dRidlerThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dKapurThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dKapurThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dJohannsenThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dJohannsenThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));

    int p3dHuangYagerThresholding_8(unsigned char*, unsigned char*, const int, const int, const int, unsigned char*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));
    int p3dHuangYagerThresholding_16(unsigned short*, unsigned char*, const int, const int, const int, unsigned short*, int (*wr_log)(const char*, ...), int (*wr_progress)(const int, ...));


    // Binary:
    int p3dClearBorderFilter3D(unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dGetRegionByCoords3D(unsigned char*, unsigned char*, const int, const int, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    
    int p3dCreateBinaryCircle(unsigned char*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dCreateBinaryCylinder(unsigned char*, const int, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dCreateBinarySphere(unsigned char*, const int, const int, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    
    
    // Utils:    
    int p3dFrom16To8(unsigned short*, unsigned char*,const int, const int, const int, unsigned short, unsigned short, int (*wr_log)(const char*, ...),int (*wr_progress)(const int, ...));

#ifdef __cplusplus
}
#endif