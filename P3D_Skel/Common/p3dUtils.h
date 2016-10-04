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

#include <limits.h>

#include "p3dCoordsList.h"

#ifdef __cplusplus
    extern "C" {
#endif

/*
	Constants:
*/
#ifndef P3D_GVF_DEFINED
	#define P3D_GVF_DEFINED
	
	#define P3D_FALSE				-1 
	#define P3D_TRUE				1 

	#define P3D_MEM_ERROR			NULL	/* Leave it NULL for simplify tests */
	#define P3D_ERROR			    0	/* Leave it NULL for simplify tests */
	#define P3D_SUCCESS				2		/* Any number */

	#define BACKGROUND				0
	#define OBJECT					UCHAR_MAX

	// Constants for 3D connectivity:
	#define CONN6   711
	#define CONN18  712
	#define CONN26  713

	// Constants for GVF Skeletonization - critical point classification:
	#define	CPT_SADDLE				200
	#define CPT_ATTRACTING_NODE		225
	#define CPT_REPELLING_NODE		250

	// This is the size of the step used to advance to the next position 
	// when following the vector field (0.2 seems to be the best value
	// and it should be multiplied for the scale factor in case of down-
	// sampling):
	#define STEP_SIZE				0.2 

	// This is the distance used to check if current point is close to a 
	// skeleton point or a critical point (0.2 seems to be the best value
	// and it should be multiplied for the scale factor in case of down-
	// sampling):
	#define CLOSENESS				0.2

	#define MIN_LENGTH			    2.0

#endif

/*
	Macros:
*/

#ifndef P3D_GVF_MACROS
	#define P3D_GVF_MACROS

	#define I(i,j,k,N,M)    ( (j)*(N) + (i) + (k)*(N)*(M) ) 
	#define MIN(x,y)        (((x) < (y))?(x):(y))
	#define MAX(x,y)        (((x) > (y))?(x):(y))

	#define EPSILON			1E-3					/* Do not modify: 1E-3 is fair */
	#define EQUAL(n1, n2)	(IS_ZERO((n1) - (n2)))
	#define IS_ZERO(n)		(((n) < EPSILON) && ((n) > -EPSILON))

	#define SIGN(n)			(IS_ZERO(n) ? 0 : ((n) < 0.00 ? -1 : 1))

	
	/* A sort of TRY-CATCH constructor: */
	#define P3D_MEM_TRY( function ) if ( (function) == P3D_MEM_ERROR) { goto MEM_ERROR; }
	#define P3D_TRY( function ) if ( (function) == P3D_ERROR) { goto MEM_ERROR; }
	

#endif

/*
	Functions:
*/
int findNeighbor ( 
	unsigned char* im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,
	const int j,
	const int k,
	coords_t* coords
	);

int countNeighbors ( 
	unsigned char* im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,
	const int j,
	const int k
	);


int isSimplePoint(   
	unsigned char* im,		// IN: Input (binary) original volume
	const int dimx,
	const int dimy, 
	const int dimz,
	const int x,
	const int y,
	const int z
	);

int p3dCrop3D_uchar2uchar (	
    unsigned char* in_im,
	unsigned char* out_im,
	const int dimx, 
	const int dimy, 
	const int dimz, 
	const int size
	);

int p3dCrop3D_ushort2ushort (	
    unsigned short* in_im,
	unsigned short* out_im,
	const int dimx, 
	const int dimy, 
	const int dimz, 
	const int size
	);

int p3dCrop3D_float2float (	
    float* in_im,
	float* out_im,
	const int dimx, 
	const int dimy, 
	const int dimz, 
	const int size
	);

int p3dZeroPadding3D_uchar2float (	
    unsigned char* in_im,
	float* out_im,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size
	);


int p3dZeroPadding3D_uchar2uchar (	
    unsigned char* in_im,
	unsigned char* out_im,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size
	);

int p3dReplicatePadding3D_uchar2uchar (	
	unsigned char* in_im,
	unsigned char* out_im,
	const int dimx,
	const int dimy,
	const int dimz,
	const int size
	);

double interpolation ( 
	float* gvf,
	int dimx, 
	int dimy, 
	int dimz, 
	const double x, 
	const double y, 
	const double z	
	);

int isBoundary ( 
	unsigned char* in_im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,    
	const int j,    
	const int k
	);

int isFullNeighborhood ( 
	unsigned char* in_im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,    
	const int j,    
	const int k
	);

int p3dSpecialPadding3D_uchar2uchar (	
	unsigned char* in_im,
	unsigned char* out_im,
	const int dimx, 
	const int dimy, 
	const int dimz, 
	const int size
	);

#ifdef __cplusplus
    }
#endif