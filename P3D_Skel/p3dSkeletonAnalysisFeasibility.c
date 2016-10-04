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

    /*char auth_code;

    //
    // Authenticate:
    //
    //auth_code = authenticate("p3dSkeletonAnalysisFeasibility");
    //if (auth_code == '0') goto AUTH_ERROR;*/

	
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

MEM_ERROR:

	if (wr_log != NULL)
	{
            wr_log ("Pore3D - Authentication error!. Program will exit.");
	}


	// Return error:
	return P3D_MEM_ERROR;
        
       /* AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }
	
    return P3D_AUTH_ERROR;*/
}

