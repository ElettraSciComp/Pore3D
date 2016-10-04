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

#include <stdlib.h>
#include <omp.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dThinning.h"
#include "Common/p3dUtils.h"


int p3dThinningSkeletonization(   
	unsigned char* in_im, 
	unsigned char* out_im, 
	const int dimx,
	const int dimy, 
	const int dimz,
	int (*wr_log)(const char*, ...)
	)
{

	unsigned char* tmp_im;

	// Padding/cropping size and dimensions
	const int a_rad = 1;					
	const int a_dimx = dimx + a_rad*2;
	const int a_dimy = dimy + a_rad*2;
	const int a_dimz = dimz + a_rad*2;	

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dThinningSkeletonization");
    if (auth_code == '0') goto AUTH_ERROR;*/


	// Start tracking computational time:
	if (wr_log != NULL)
	{
		p3dResetStartTime(); 
		wr_log ("Pore3D - Performing thinning skeletonization..." );
	}



	// Init output voume with input volume values zero padded:
	P3D_TRY( tmp_im = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) ) );
	P3D_TRY( p3dZeroPadding3D_uchar2uchar ( in_im, tmp_im, dimx, dimy, dimz, a_rad ) );

	// Call in-place version:
	P3D_TRY( p3dThinning ( tmp_im, a_dimx, a_dimy, a_dimz ) );

	// Crop output:
	P3D_TRY( p3dCrop3D_uchar2uchar ( tmp_im, out_im, a_dimx, a_dimy, a_dimz, a_rad ) );




	// Print elapsed time (if required):
	if (wr_log != NULL)
	{			
		wr_log ("Pore3D - Thinning skeletonization performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}	
	
    // Release resources:		
	if ( tmp_im != NULL ) free(tmp_im);	

	// Return OK:
	return P3D_SUCCESS;
	

MEM_ERROR:

	if (wr_log != NULL)
	{
		wr_log ("Pore3D - Not enough (contiguous) memory. Program will exit.");
	}

    // Release resources:	
	if ( tmp_im != NULL ) free(tmp_im);	

	
	return P3D_MEM_ERROR;
        
      /*  AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

