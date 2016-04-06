#include <math.h>
#include <string.h>

#include "p3dComputeCoreSkeleton.h"

#include "../Common/p3dCoordsT.h"
#include "../Common/p3dUtils.h"

//
// Procedures
//

int stabilityReached ( fcoords_list_t segm_list, fcoords_t pt )
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


int skelPointReached ( 	
	fcoords_list_t*   skel_point_list,
	fcoords_t         curr_crit_point,
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
		if (!( EQUAL( curr_crit_point.x, pt.x) && 
			   EQUAL( curr_crit_point.y, pt.y) &&
			   EQUAL( curr_crit_point.z, pt.z) ))
		{				
			a = fabs( pt.x - tmp_point.x );
			b = a + fabs( pt.y - tmp_point.y );
			c = a + b + fabs( pt.z - tmp_point.z );

			if( (a < close_dist ) && ( b < close_dist ) && ( c < close_dist ) )
				return P3D_TRUE;		
		}

		// The list will point on the next element:
		tmp_list = tmp_list->next;		
	}

	return P3D_FALSE;
}


int critPointReached ( 	
	crit_point_list_t crit_point_list,
	fcoords_t         curr_crit_point,
	fcoords_t         pt,
	const double      close_dist
	)
{
	crit_point_list_t tmp_list;
	crit_point_t crit_point;	  
	double a, b, c;
  
	// Print critical points:
	tmp_list = crit_point_list;

	// Scan list:
	while ( crit_point_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current element:
		crit_point = tmp_list->elem;

		// Skip current critical point (it would be better to test an ID...):
		if (!( EQUAL( curr_crit_point.x, crit_point.x) && 
			   EQUAL( curr_crit_point.y, crit_point.y) &&
			   EQUAL( curr_crit_point.z, crit_point.z) ))
		{		
			a = fabs( pt.x - crit_point.x );
			b = a + fabs( pt.y - crit_point.y );
			c = a + b + fabs( pt.z - crit_point.z );

			if( (a < close_dist ) && ( b < close_dist ) && ( c < close_dist ) )
				return P3D_TRUE;					
		}

		// The list will point on the next element:
		tmp_list = tmp_list->next;		
	}

	return P3D_FALSE;
}

int followStreams ( 
	fcoords_list_t*   skel_point_list,
	crit_point_list_t crit_point_list,
	fcoords_t		  crit_point, 
	fcoords_t         verse, 
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
	double out_x, out_y, out_z;
	double len;	
	
	fcoords_t fpoint1, fpoint2;

	fcoords_list_t tmp_list;


	// Init list:
	fcoords_list_init ( &tmp_list );

	// Initialize starting point and put it into the temporary list for current segment:
	fpoint1.x = crit_point.x;
	fpoint1.y = crit_point.y;
	fpoint1.z = crit_point.z;	

	P3D_TRY ( fcoords_list_push( &tmp_list, fpoint1 ) );

	// Initialize next point and put it into the temporary list for current segment:
	fpoint2.x = crit_point.x + verse.x*step;
	fpoint2.y = crit_point.y + verse.y*step;
	fpoint2.z = crit_point.z + verse.z*step;

	
	// Start following path until exit conditions:
	//    i. Critical point reached
	//   ii. Skeleton point reached
	//  iii. Stability reached (safety condition)
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
	}
	while ( ( skelPointReached( skel_point_list, crit_point, fpoint2, close_dist ) == P3D_FALSE ) &&
			( critPointReached( crit_point_list, crit_point, fpoint2, close_dist ) == P3D_FALSE ) &&
			( stabilityReached( tmp_list, fpoint2 ) == P3D_FALSE ) );
	
	// Current segment is added to skeleton:
	while ( fcoords_list_isempty( tmp_list ) == P3D_FALSE )
	{
		P3D_TRY ( fcoords_list_push( skel_point_list, fcoords_list_pop( &tmp_list )));
	}

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


int p3dComputeCoreSkeleton (
    crit_point_list_t crit_point_list,
	fcoords_list_t*   skel_point_list,
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
	crit_point_list_t  tmp_list;
	crit_point_t       crit_point;

	fcoords_t  verse;
	fcoords_t  point;	

	//
	// Following the streamlines starting at a saddle point in the direction of the
	// positive eigenvector(s):
	//  	
	tmp_list = crit_point_list;

	// Scan list:
	while ( crit_point_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current critical point:
		crit_point = tmp_list->elem;

		// If current critical point is a saddle point:
		if ( crit_point.type == CPT_SADDLE )
		{	
			
			point.x = crit_point.x;
			point.y = crit_point.y;
			point.z = crit_point.z;

			// Get the direction pointed by the positive first eigenvector:
			if( crit_point.eval0 > 0 ) 
			{

				// UP direction given by the eigenvector:
				verse.x = crit_point.evect0_x;
				verse.y = crit_point.evect0_y;
				verse.z = crit_point.evect0_z;

				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ) );
	
				// DOWN direction of the eigenvector:
				verse.x = - crit_point.evect0_x;
				verse.y = - crit_point.evect0_y;
				verse.z = - crit_point.evect0_z;
	
				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ) );
			}

			// Get the direction pointed by the positive second eigenvector:
			if( crit_point.eval1 > 0 ) 
			{	
				// UP direction given by the eigenvector:
				verse.x = crit_point.evect1_x;
				verse.y = crit_point.evect1_y;
				verse.z = crit_point.evect1_z;
	
				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ));

				// DOWN direction of the eigenvector:
				verse.x = - crit_point.evect1_x;
				verse.y = - crit_point.evect1_y;
				verse.z = - crit_point.evect1_z;
	
				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ));
			}

			// Get the direction pointed by the positive third eigenvector:
			if( crit_point.eval2 > 0 ) 
			{	
				// UP direction given by the eigenvector:
				verse.x = crit_point.evect2_x;
				verse.y = crit_point.evect2_y;
				verse.z = crit_point.evect2_z;
	
				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ));

				// DOWN direction of the eigenvector:
				verse.x = - crit_point.evect2_x;
				verse.y = - crit_point.evect2_y;
				verse.z = - crit_point.evect2_z;
	
				P3D_TRY ( followStreams ( skel_point_list, crit_point_list, point, verse, gvf_x, gvf_y, 
					gvf_z, dimx, dimy, dimz, step, close_dist ));
			}			
		}	

		// The list will point on the next element:
		tmp_list = tmp_list->next;
	}


	// Return OK:
	return P3D_SUCCESS;

MEM_ERROR:

	return P3D_MEM_ERROR;
}