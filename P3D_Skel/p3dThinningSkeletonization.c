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

