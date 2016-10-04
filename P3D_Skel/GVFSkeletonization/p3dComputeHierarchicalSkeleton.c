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

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "p3dComputeHierarchicalSkeleton.h"
#include "p3dHighDivPointList.h"

#include "../Common/p3dUtils.h"


//
// Procedures
//

int _p3dComputeHierarchicalSkeleton_stabilityReached ( fcoords_list_t segm_list, fcoords_t pt )
{
	fcoords_list_t tmp_list;
	fcoords_t      tmp_elem;

	// Scan list:
	tmp_list = segm_list;

	// If input point is already present into the current skeleton segment probably we are 
	// chasing our own tail:
	while ( fcoords_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current element:
		tmp_elem = tmp_list->elem;	

		if ( EQUAL( tmp_elem.x, pt.x) && 
			 EQUAL( tmp_elem.y, pt.y) &&
			 EQUAL( tmp_elem.z, pt.z) )
				return P3D_TRUE;

		// The list will point on the next element:
		tmp_list = tmp_list->next;	
	}

	return P3D_FALSE;
}

int _p3dComputeHierarchicalSkeleton_skelPointReached ( 	
	fcoords_list_t*   skel_point_list,
	fcoords_t         orig_point,
	fcoords_t         pt,
	const double      close_dist
	)
{
	fcoords_list_t tmp_list;
	fcoords_t tmp_point;
	double a, b, c;
	
  
	// Print critical points:
	tmp_list = (*skel_point_list);

	// Scan list:
	while ( fcoords_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current element:
		tmp_point = tmp_list->elem;

		// Skip current critical point (it would be better to test an ID...):
		if (!( EQUAL( orig_point.x, pt.x) && 
			   EQUAL( orig_point.y, pt.y) &&
			   EQUAL( orig_point.z, pt.z) ))
		{				
			a = fabs( pt.x - tmp_point.x );
			b = a + fabs( pt.y - tmp_point.y );
			c = a + b + fabs( pt.z - tmp_point.z );

			if( (a < close_dist ) && ( b < close_dist ) && ( c < close_dist ) )
			{				
				return P3D_TRUE;				
			}
		}

		// The list will point to the next element:
		tmp_list = tmp_list->next;		
	}

	return P3D_FALSE;
}



int _p3dComputeHierarchicalSkeleton_followStreams ( 
	fcoords_list_t*   skel_point_list,
	fcoords_t		  point, 
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,
	const double step,
	const double close_dist
	)
{	
	int flag;

	double out_x, out_y, out_z;
	double len;	
	
	fcoords_t fpoint1, fpoint2;

	fcoords_list_t tmp_list;


	// Init list:
	fcoords_list_init ( &tmp_list );

	// Initialize starting point and put it into the temporary list for current segment:
	fpoint2.x = point.x;
	fpoint2.y = point.y;
	fpoint2.z = point.z;	
	
	// Start following path until exit conditions:
	//   i. Skeleton point reached
	//  ii. Stability reached (safety condition)
	do {

		// Set current point and put it into list:
		fpoint1.x = fpoint2.x;
		fpoint1.y = fpoint2.y;
		fpoint1.z = fpoint2.z;		

		P3D_TRY ( fcoords_list_push( &tmp_list, fpoint1 ) );

		// Interpolate gradient vector flow:
		out_x = interpolation ( gvf_x, dimx, dimy, dimz, fpoint1.x, fpoint1.y, fpoint1.z );
		out_y = interpolation ( gvf_y, dimx, dimy, dimz, fpoint1.x, fpoint1.y, fpoint1.z );
		out_z = interpolation ( gvf_z, dimx, dimy, dimz, fpoint1.x, fpoint1.y, fpoint1.z );
		
		// Normalize:
		len = sqrt((out_x * out_x) + (out_y * out_y) + (out_z * out_z));

		if(len > 0.00) {
			out_x = out_x / len;
			out_y = out_y / len;
			out_z = out_z / len;
		}

		// Increment path following:		
		fpoint2.x = fpoint1.x + out_x*step;
		fpoint2.y = fpoint1.y + out_y*step;
		fpoint2.z = fpoint1.z + out_z*step;

		flag = _p3dComputeHierarchicalSkeleton_skelPointReached( skel_point_list, point, fpoint2, close_dist );
	}
	while ( ( flag == P3D_FALSE ) &&
			( _p3dComputeHierarchicalSkeleton_stabilityReached( tmp_list, fpoint2 ) == P3D_FALSE ) );
	

	// Current segment is added to skeleton if long enough (spurious voxels may occur) 
	// and connected to core skeleton:
	if ( ( flag == P3D_TRUE) && ( sqrt( (fpoint1.x - point.x) * (fpoint1.x - point.x) + 
		       (fpoint1.y - point.y) * (fpoint1.y - point.y) + 
			   (fpoint1.z - point.z) * (fpoint1.z - point.z)) > MIN_LENGTH ) )
	{
		while ( fcoords_list_isempty( tmp_list ) == P3D_FALSE )
		{
			P3D_TRY ( fcoords_list_push( skel_point_list, fcoords_list_pop( &tmp_list )));
		}
	}
	/*else
		printf("Qualcosa di diverso ho fatto!\n");*/

	// Return OK:
	return P3D_SUCCESS;

MEM_ERROR:

	// Release resources of temporary list:
	while ( fcoords_list_isempty( tmp_list ) == P3D_FALSE )
	{
		fcoords_list_pop( &tmp_list );
	}

	return P3D_MEM_ERROR;
}


int p3dComputeHierarchicalSkeleton (    
	highDiv_point_list_t*   highDiv_point_list,
	fcoords_list_t*         skel_point_list,
    float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,
	const double step,
	const double close_dist
	)
{
	highDiv_point_t  point;	
	fcoords_t        fpoint;

	//
	// Following the streamlines starting at an high divergence point:
	//  

	// Scan list (list is deleted at the end of the loop):
	while ( highDiv_point_list_isempty( *highDiv_point_list ) == P3D_FALSE )
	{
		// Get current critical point:
		point = highDiv_point_list_pop( highDiv_point_list );

		// "Convert" highDiv_point_t to fcoords_point:
		fpoint.x = point.x;
		fpoint.y = point.y;
		fpoint.z = point.z;

		P3D_TRY ( _p3dComputeHierarchicalSkeleton_followStreams ( skel_point_list, fpoint, gvf_x, gvf_y, gvf_z, 
			dimx, dimy, dimz, step, close_dist ) );		
	}


	// Return OK:
	return P3D_SUCCESS;

MEM_ERROR:

	// Release resources:
	while ( highDiv_point_list_isempty( *highDiv_point_list ) == P3D_FALSE )
	{
		highDiv_point_list_pop( highDiv_point_list );
	}

	// Return error:
	return P3D_MEM_ERROR;
}
