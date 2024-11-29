#include <stdlib.h>
#include <omp.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dUtils.h"


int p3dSkeletonAnalysisFeasibility (
	 unsigned char*  in_im,	// IN: binary volume
	 unsigned char*	 sk_im,    // IN: skeleton mask
          double* ratio,				// OUT: percentage of total pore space
								// occupied by the medial axis.
	 const int dimx,
	 const int dimy, 
	 const int dimz,	
	 int (*wr_log)(const char*, ...)
	 )
{    
        int i,j,k;
        double porosityct;
	double skeletonct;
	
	
	// Start tracking computational time:
	if (wr_log != NULL)
	{
		p3dResetStartTime(); 
		wr_log ("Pore3D - Computing skeleton analysis feasibility..." );
	}

	// Initialize variables:
	porosityct = 0.0;
	skeletonct = 0.0;
	
	// Scan volumes:
	#pragma omp parallel for private(i, j)  reduction (+ : porosityct, skeletonct )
        for( k = 0; k < dimz; k++ )
		for( j = 0; j < dimy; j++ )
			for( i = 0; i < dimx; i++ )            
            {  
                if ( in_im[ I(i,j,k,dimx,dimy) ] == BACKGROUND )
                    porosityct = porosityct + 1.0;
                
                if ( sk_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND )
                    skeletonct = skeletonct + 1.0;
            }	


	// Return percentage of total pore space occupied by the medial axis:
	*ratio = skeletonct / ((double) (porosityct));

	if (wr_log != NULL)
	{
		wr_log ("\tPercentage of total space occupied by the skeleton [%]: %0.3f.", (*ratio * 100) );			
	}


	// Print elapsed time (if required):
	if (wr_log != NULL)
	{
		wr_log ("Pore3D - Skeleton analysis feasibility computed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}	

	// Return OK:
	return P3D_SUCCESS;

}

