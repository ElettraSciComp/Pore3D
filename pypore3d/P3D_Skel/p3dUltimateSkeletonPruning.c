#include <stdlib.h>
#include <omp.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dThinning.h"
#include "Common/p3dUtils.h"


int p3dUltimateSkeletonPruning(   
	unsigned char* in_im,		// IN: Input (binary) skeleton
	unsigned char* out_im,		// OUT: Labeled skeleton
	const int dimx,
	const int dimy, 
	const int dimz,
	const int iterative,
	int (*wr_log)(const char*, ...)
	)
{
	int a,b,c;

	int noMoreChanges = P3D_FALSE; // Flag for iterative version of pruning

	// Dimensions for padded/cropped volumes:
	int a_dimx, a_dimy, a_dimz;
	int a_rad;

	// Padded and cropped temporary images:
	unsigned char* tmp_im;
	unsigned char* tmp_im2;

	
	// Start tracking computational time:
	if (wr_log != NULL)
	{
		p3dResetStartTime(); 
		wr_log ("Pore3D - Performing ultimate pruning..." );
	}

	// Create temporary input replicate-padded:
	a_rad = 1;

	// Compute dimensions of padded REV:
	a_dimx = dimx + a_rad*2;
	a_dimy = dimy + a_rad*2;
	a_dimz = dimz + a_rad*2;

	// Initialize input:
	P3D_MEM_TRY( tmp_im = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) ) );
	P3D_MEM_TRY( tmp_im2 = (unsigned char*) calloc( a_dimx*a_dimy*a_dimz,sizeof(unsigned char) ) );

	P3D_TRY( p3dZeroPadding3D_uchar2uchar ( in_im, tmp_im2, dimx, dimy, dimz, a_rad ) );


	// Loop for iterations:
	while ( noMoreChanges == P3D_FALSE)
	{	
		noMoreChanges = P3D_TRUE;

		// Perform skeleton labeling:
		P3D_TRY( p3dSkeletonLabeling ( tmp_im2, tmp_im, a_dimx, a_dimy, a_dimz, NULL ) );

		// Volume scanning:
		#pragma omp parallel for private(a, b) 
		for( c = a_rad; c < (a_dimz - a_rad); c++ )  
			for( b = a_rad; b < (a_dimy - a_rad); b++ )
				for( a = a_rad; a < (a_dimx - a_rad); a++ )
				{
					// If we're on a skeleton voxel different than one of a NODETOEND branch:
					if ( ( tmp_im[ I( a, b, c, a_dimx, a_dimy ) ] != BACKGROUND ) && 
						( tmp_im[ I( a, b, c, a_dimx, a_dimy ) ] != NODETOEND_LABEL ) ) 					
					{
						tmp_im2[ I( a, b, c, a_dimx, a_dimy) ] = OBJECT;
					}
					if ( tmp_im[ I( a, b, c, a_dimx, a_dimy ) ] == NODETOEND_LABEL ) 
					{
						noMoreChanges = P3D_FALSE;
					}
				}

		// Remove isolated end points created by previous pruning:
		#pragma omp parallel for private(a, b) 
		for( c = a_rad; c < (a_dimz - a_rad); c++ )  
			for( b = a_rad; b < (a_dimy - a_rad); b++ )
				for( a = a_rad; a < (a_dimx - a_rad); a++ )
				{
					// If we're on a skeleton voxel:
					if ( tmp_im2[ I( a, b, c, a_dimx, a_dimy ) ] == OBJECT ) 					
					{
						// Check Neighborhood:
						if ( countNeighbors ( tmp_im2, a_dimx, a_dimy, a_dimz, a, b, c )  == 0 )
						{
							tmp_im2[ I( a, b, c, a_dimx, a_dimy ) ] = BACKGROUND;
						}
					}
				}

		// Re-skeletonize:
		p3dThinning ( tmp_im2, a_dimx, a_dimy, a_dimz );

		// Avoid iterations if this option was undesidered:
		if ( iterative == P3D_FALSE )
			noMoreChanges = P3D_TRUE;

	}

	// Crop output:
	P3D_TRY( p3dCrop3D_uchar2uchar ( tmp_im2, out_im, a_dimx, a_dimy, a_dimz, a_rad ) );	

	
	// Print elapsed time (if required):
	if (wr_log != NULL)
	{
		wr_log ("Pore3D - Ultimate skeleton pruning performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}	
	
    // Release resources:		
	if ( tmp_im != NULL ) free(tmp_im);
	if ( tmp_im2 != NULL ) free(tmp_im2);

	// Return OK:
	return P3D_SUCCESS;
	

MEM_ERROR:

	if (wr_log != NULL)
	{
		wr_log ("Pore3D - Not enough (contiguous) memory. Program will exit.");
	}

    // Release resources:	
	if ( tmp_im != NULL )  free ( tmp_im );
	if ( tmp_im2 != NULL ) free ( tmp_im2 );

	
	return P3D_ERROR;
        
}



