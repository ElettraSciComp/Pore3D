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
#include <string.h>
#include <omp.h>
#include <math.h>
#include <limits.h>

#include "p3dGetHighDivPoints.h"
#include "p3dHighDivPointList.h"

#include "../Common/p3dUtils.h"

#define SEARCH_GRID		1
#define CELL_SIZE		1.00 / SEARCH_GRID


int  _p3dGetHighDivPoints_cmp (const void *a, const void *b)
{
    const highDiv_point_t *da = (const highDiv_point_t *) a;
    const highDiv_point_t *db = (const highDiv_point_t *) b;
     
    return (da->div > db->div) - (da->div < db->div);
}

int _p3dGetHighDivPoints_pointIsCloseToGroup ( highDiv_point_t pt, highDiv_point_list_t list) 
{
	highDiv_point_list_t tmp_list;
	highDiv_point_t      tmp_point;

	tmp_list = list;

	while ( highDiv_point_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current element:
		tmp_point = tmp_list->elem;
		
		if(	(fabs( pt.x - tmp_point.x ) <= 1) &&
			(fabs( pt.y - tmp_point.y ) <= 1) &&
			(fabs( pt.z - tmp_point.z ) <= 1) )
		{	
			// Input point is "close" to a point already present into list:
			return P3D_TRUE;
		}

		// The list will point on the next element:
		tmp_list = tmp_list->next;	
    }
  
	// The input point is not "close" to a point already present into list:
	return P3D_FALSE;
}


int p3dGetHighDivPoints (
	unsigned char* in_im,
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,	
	const double thresh,
	highDiv_point_list_t* highDiv_point_list	
	)
{
	int a_rad;
	int i, j, k;
	int a, b, c;
	double x, y, z;

	double* in_gvf_x;
	double* in_gvf_y;
	double* in_gvf_z;

	double maxDiv = LONG_MIN;
	double minDiv = LONG_MAX;
	double div;
	double threshold;

	double vdist = (CELL_SIZE) / 2.00;
  
	highDiv_point_t   highDivPoint;
	highDiv_point_t*  vtmp;
	int ct_HighDivPoints = 0;

	

			

	P3D_TRY ( in_gvf_x = (double*) malloc(6*sizeof(double)));
	P3D_TRY ( in_gvf_y = (double*) malloc(6*sizeof(double)));
	P3D_TRY ( in_gvf_z = (double*) malloc(6*sizeof(double)));
  	
	a_rad = 1; // An high divergence point could not occur on boundaries
	

	//
	// Compute divergence throughout the dataset:
	//
	for( k = a_rad; k < (dimz - a_rad); k++ )    
		for( j = a_rad; j < (dimy - a_rad); j++ )
			for( i = a_rad; i < (dimx - a_rad); i++ )
			{
			
				// Skip voxels on object boundary (6-connection used):
				if ( ( in_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND ) &&
					 ( isFullNeighborhood ( in_im, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) 
					 )
				{	
					for ( c = 0; c < SEARCH_GRID; c++ ) 
						for ( b = 0; b < SEARCH_GRID; b++ ) 
							for ( a = 0; a < SEARCH_GRID; a++ ) 
							{
								x = i + (a * CELL_SIZE);
								y = j + (b * CELL_SIZE);
								z = k + (c * CELL_SIZE);

								// Interpolate values in the neighborhood:
								in_gvf_x[0] = interpolation( gvf_x, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_x[1] = interpolation( gvf_x, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_x[2] = interpolation( gvf_x, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_x[3] = interpolation( gvf_x, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_x[4] = interpolation( gvf_x, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_x[5] = interpolation( gvf_x, dimx, dimy, dimz, x, y, z - vdist );
										
								in_gvf_y[0] = interpolation( gvf_y, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_y[1] = interpolation( gvf_y, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_y[2] = interpolation( gvf_y, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_y[3] = interpolation( gvf_y, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_y[4] = interpolation( gvf_y, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_y[5] = interpolation( gvf_y, dimx, dimy, dimz, x, y, z - vdist );

								in_gvf_z[0] = interpolation( gvf_z, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_z[1] = interpolation( gvf_z, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_z[2] = interpolation( gvf_z, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_z[3] = interpolation( gvf_z, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_z[4] = interpolation( gvf_z, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_z[5] = interpolation( gvf_z, dimx, dimy, dimz, x, y, z - vdist );

		      
								div = (	(in_gvf_x[0] - in_gvf_x[1]) + (in_gvf_y[2] - in_gvf_y[3]) + 
										(in_gvf_z[4] - in_gvf_z[5])) / (2.0 * vdist);

								// Add the point to the list
								/*highDivPoint.x = x;
								highDivPoint.y = y;
								highDivPoint.z = z;

								highDivPoint.div = div;

								P3D_TRY ( highDiv_point_list_push( highDiv_point_list, highDivPoint ) );

								ct_HighDivPoints++;*/
		      
						     
				      			if( div > maxDiv )	maxDiv = div;
								if( div < minDiv )	minDiv = div;
							}
				}
			}
  
	// case 1:
	// take <perc> percent of the lowest negative value:
	// !! have to change the comparison
	threshold = maxDiv - minDiv;
	threshold = thresh * threshold;
	threshold = minDiv + threshold;  


	//
	// Re-compute divergence throughout the dataset in order to save memory occupation:
	//
	for( k = a_rad; k < (dimz - a_rad); k++ )    
		for( j = a_rad; j < (dimy - a_rad); j++ )
			for( i = a_rad; i < (dimx - a_rad); i++ )
			{
			
				// Skip voxels on object boundary (6-connection used):
				if ( ( in_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND ) &&
					 ( isFullNeighborhood ( in_im, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) 
					 )
				{					
		
					for ( c = 0; c < SEARCH_GRID; c++ ) 
						for ( b = 0; b < SEARCH_GRID; b++ ) 
							for ( a = 0; a < SEARCH_GRID; a++ ) 
							{
								x = i + (a * CELL_SIZE);
								y = j + (b * CELL_SIZE);
								z = k + (c * CELL_SIZE);

								// Interpolate values in the neighborhood:
								in_gvf_x[0] = interpolation( gvf_x, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_x[1] = interpolation( gvf_x, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_x[2] = interpolation( gvf_x, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_x[3] = interpolation( gvf_x, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_x[4] = interpolation( gvf_x, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_x[5] = interpolation( gvf_x, dimx, dimy, dimz, x, y, z - vdist );
										
								in_gvf_y[0] = interpolation( gvf_y, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_y[1] = interpolation( gvf_y, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_y[2] = interpolation( gvf_y, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_y[3] = interpolation( gvf_y, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_y[4] = interpolation( gvf_y, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_y[5] = interpolation( gvf_y, dimx, dimy, dimz, x, y, z - vdist );

								in_gvf_z[0] = interpolation( gvf_z, dimx, dimy, dimz, x + vdist, y, z );
								in_gvf_z[1] = interpolation( gvf_z, dimx, dimy, dimz, x - vdist, y, z );
								in_gvf_z[2] = interpolation( gvf_z, dimx, dimy, dimz, x, y + vdist, z );
								in_gvf_z[3] = interpolation( gvf_z, dimx, dimy, dimz, x, y - vdist, z );
								in_gvf_z[4] = interpolation( gvf_z, dimx, dimy, dimz, x, y, z + vdist );
								in_gvf_z[5] = interpolation( gvf_z, dimx, dimy, dimz, x, y, z - vdist );

		      
								div = (	(in_gvf_x[0] - in_gvf_x[1]) + (in_gvf_y[2] - in_gvf_y[3]) + 
										(in_gvf_z[4] - in_gvf_z[5])) / (2.0 * vdist);

								if( div <= threshold ) 
								{
									// Add the point to the list
									highDivPoint.x = x;
									highDivPoint.y = y;
									highDivPoint.z = z;

									highDivPoint.div = div;

									P3D_TRY ( highDiv_point_list_push( highDiv_point_list, highDivPoint ) );

									ct_HighDivPoints++;
								}
							}
				}
			}
	
	//
	// Convert to array the list of high divergence points:
	//
	P3D_TRY ( vtmp = highDiv_point_list_toarray ( highDiv_point_list, ct_HighDivPoints ) );	
	
	//
	// Sort the array according to the divergence value;
	// 
	qsort(vtmp, ct_HighDivPoints, sizeof(highDiv_point_t), _p3dGetHighDivPoints_cmp);
	
  
	//
	// Remove duplicate points according to a "closeness" criteria:	
	//
	for( i = 0; i < ct_HighDivPoints; i++ ) 
	{
		
		/*if ( ( vtmp[i].div < threshold ) && ( _p3dGetHighDivPoints_pointIsCloseToGroup( 
			   vtmp[i], *highDiv_point_list ) == P3D_FALSE ) )*/
		if ( _p3dGetHighDivPoints_pointIsCloseToGroup( vtmp[i], *highDiv_point_list ) == P3D_FALSE )
		{
			// Add the point to the output list:
			P3D_TRY ( highDiv_point_list_push ( highDiv_point_list, vtmp[i] ) );
		}
	}
	
  
	// Release resources:	
	if ( vtmp != NULL ) free ( vtmp );
	if ( in_gvf_x != NULL ) free ( in_gvf_x );
	if ( in_gvf_y != NULL ) free ( in_gvf_y );
	if ( in_gvf_z != NULL ) free ( in_gvf_z );
	
	// Return success:
    return P3D_SUCCESS;

MEM_ERROR:

	// Release resources and return error:
	if ( vtmp != NULL ) free (vtmp);
		
	if ( in_gvf_x != NULL ) free( in_gvf_x );
	if ( in_gvf_y != NULL ) free( in_gvf_y );
	if ( in_gvf_z != NULL ) free( in_gvf_z );

	// Return error:
	return P3D_MEM_ERROR;
}