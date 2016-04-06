#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "p3dCriticalPoints.h"
#include "p3dComputeEigenVal.h"

#include "../Common/p3dUtils.h"


//
// Constants:
//
#define CELL_SUBDIVISION_FACTOR	1048576.00
#define NR_OF_SIGN_CHANGES		3



int isChangeInSignNeighborhood ( 
	float* gvf, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,    
	const int j,    
	const int k
	)
{
	int a,b,c,s; 

	// Initializations:

	s = SIGN( gvf[ I(i,j,k,dimx,dimy) ] );

	// Check 8-neighborhood:

	c = k + 1; b = j; a = i;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k; b = j + 1; a = i;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k; b = j; a = i + 1;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k + 1; b = j + 1; a = i;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k + 1; b = j; a = i + 1;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k; b = j + 1; a = i + 1;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	c = k + 1; b = j + 1; a = i + 1;
	if ( SIGN( gvf[ I(a,b,c,dimx,dimy) ]) != s ) return P3D_TRUE;

	// No change in sign but sign is 0, return P3D_TRUE anyway:
	if( s == 0 ) return P3D_TRUE;     

	// No change in sign within neighborhood:
	return P3D_FALSE;
}


int changeInSign(
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	// ncols
	const int dimy,	// nrows
	const int dimz,	// nplanes
	const int i,
	const int j,
	const int k
	) 
{
  
	int ct = 0;
  
	// Voxels having both positive and negative values for every 
	// vector component in its neighborhood are assumed to be
	// critical points. If one component is 0 in all vectors, it 
	// should be considered as a change in sign. 

	if( isChangeInSignNeighborhood( gvf_x, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) ct++;
	if( isChangeInSignNeighborhood( gvf_y, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) ct++;
	if( isChangeInSignNeighborhood( gvf_z, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) ct++;
    
	if( ct == NR_OF_SIGN_CHANGES ) return P3D_TRUE;
      
	return P3D_FALSE;
}





int isChangeInSignNeighborhoodInterp ( 
	float* gvf, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const double x,    
	const double y,    
	const double z,
	const double cellsize
	)
{
	int s; 

	// Initializations:
	s = SIGN( interpolation ( gvf, dimx, dimy, dimz, x, y, z ) );

	// Check 8-neighborhood:
	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x + cellsize, y, z )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x, y + cellsize, z )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x, y, z + cellsize )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x + cellsize, y + cellsize, z )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x + cellsize, y, z + cellsize )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x, y + cellsize, z + cellsize )) != s ) return P3D_TRUE;

	if ( SIGN( interpolation ( gvf, dimx, dimy, dimz, x + cellsize, y + cellsize, z + cellsize )) != s ) return P3D_TRUE;

	// No change in sign but sign is 0, return P3D_TRUE anyway:
	if( s == 0 ) return P3D_TRUE;     

	// No change in sign within neighborhood:
	return P3D_FALSE;
}


int changeInSignInterp(
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	// ncols
	const int dimy,	// nrows
	const int dimz,	// nplanes
	const double i,
	const double j,
	const double k,
	const double cellsize
	) 
{
  
	int ct = 0;
  
	// Voxels having both positive and negative values for every 
	// vector component in its neighborhood are assumed to be
	// critical points. If one component is 0 in all vectors, it 
	// should be considered as a change in sign. 

	if( isChangeInSignNeighborhoodInterp( gvf_x, dimx, dimy, dimz, i, j, k, cellsize ) == P3D_TRUE ) ct++;
	if( isChangeInSignNeighborhoodInterp( gvf_y, dimx, dimy, dimz, i, j, k, cellsize ) == P3D_TRUE ) ct++;
	if( isChangeInSignNeighborhoodInterp( gvf_z, dimx, dimy, dimz, i, j, k, cellsize ) == P3D_TRUE ) ct++;
    
	if( ct == NR_OF_SIGN_CHANGES ) return P3D_TRUE;
      
	return P3D_FALSE;
}








int getCriticalPointInFloatCell ( 
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	// ncols
	const int dimy,	// nrows
	const int dimz,	// nplanes
	double x, 
	double y, 
	double z, 
	double cellsize,
	crit_point_list_t* crit_point_list
	)
{
	crit_point_t crit_point;
	double in_gvf_x, in_gvf_y, in_gvf_z;
	int err_code;
	  
	// Get the first sub-cell:
	in_gvf_x = interpolation ( gvf_x, dimx, dimy, dimz, x, y, z );
	in_gvf_y = interpolation ( gvf_y, dimx, dimy, dimz, x, y, z );
	in_gvf_z = interpolation ( gvf_z, dimx, dimy, dimz, x, y, z );

	
	// If first vertex vector (coresponding to (x, y, z)) is (0, 0, 0)
	// then return (x, y, z) as the critical point:
	if ( (IS_ZERO(in_gvf_x)) && (IS_ZERO(in_gvf_y)) && (IS_ZERO(in_gvf_z)) )
    {
		// Add critical point into list:
		crit_point.x = x;
		crit_point.y = y;
		crit_point.z = z;

		P3D_TRY( crit_point_list_push ( crit_point_list, crit_point));

		return P3D_TRUE;
    }
	else
	{	  
		// The cell is a candidate cell if there is a change of sign in one of the vector 
		// components among all eight vertices of the cell:
		if ( changeInSignInterp( gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, x, y, z, cellsize) == P3D_TRUE ) 
		{
			// Stop here if maximum subdivision reached and assume
			// that the critical point is in the center of the cell:
			if( cellsize <= (1.00 / CELL_SUBDIVISION_FACTOR) ) 
			{
				// Add critical point into list:
				crit_point.x = x + cellsize/ 2.0;
				crit_point.y = y + cellsize/ 2.0;
				crit_point.z = z + cellsize/ 2.0;

				P3D_TRY ( crit_point_list_push ( crit_point_list, crit_point ));

				return P3D_TRUE;
			}
			else
			{    
				// We are on a candidate cell but it's not small enough, so we are
				// dividing the cell in 8 subcells and try to find the critical point 
				// in one of the subcells (recursive programming used):
				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x, y, z, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x + cellsize/2.0, y, z, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x, y + cellsize/2.0, z, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x + cellsize/2.0, y + cellsize/2.0, z, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x, y, z + cellsize/2.0, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x + cellsize/2.0, y, z + cellsize/2.0, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x, y + cellsize/2.0, z + cellsize/2.0, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

				err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, 
					x + cellsize/2.0, y + cellsize/2.0, z + cellsize/2.0, cellsize/2.0, crit_point_list );
				if ( err_code == P3D_TRUE ) return P3D_TRUE;
				if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

			}
		}
	}	

	return P3D_FALSE;

MEM_ERROR:

	return P3D_MEM_ERROR;
}



int p3dGetCriticalPoints (
	unsigned char* in_im,
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,	
	crit_point_list_t* crit_point_list	
	)
{
	int a_rad;
	int i,j,k;

	crit_point_t crit_point;

	int err_code;



	a_rad = 1; // A critical point could not occur on boundaries.
	
	//
	// Determine critical points according to the sign of GVF.
	//

	for( k = a_rad; k < (dimz - a_rad); k++ )    
		for( j = a_rad; j < (dimy - a_rad); j++ )
			for( i = a_rad; i < (dimx - a_rad); i++ )
			{
				// Skip voxels on object boundary (6-connection used) :
				if ( ( in_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND ) &&
					 ( isFullNeighborhood ( in_im, dimx, dimy, dimz, i, j, k ) == P3D_TRUE ) 
					 )
				{
					if(  ( IS_ZERO( gvf_x[ I(i,j,k,dimx,dimy) ] )) &&
						 ( IS_ZERO( gvf_y[ I(i,j,k,dimx,dimy) ] )) &&
						 ( IS_ZERO( gvf_z[ I(i,j,k,dimx,dimy) ] )) )
					{
						// Add critical point into list:
						crit_point.x = (double) i;
						crit_point.y = (double) j;
						crit_point.z = (double) k;
									
						P3D_TRY ( crit_point_list_push ( crit_point_list, crit_point) );
					}
					else if ( changeInSign( gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i, j, k ) == P3D_TRUE )
					{
						// Explore subcells:		
						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i, j, k, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i, j, k + 0.5, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i, j + 0.5, k, 0.5, crit_point_list ); 
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i + 0.5, j, k, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i, j + 0.5, k + 0.5, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i + 0.5, j, k + 0.5, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i + 0.5, j + 0.5, k, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;

						err_code = getCriticalPointInFloatCell (gvf_x, gvf_y, gvf_z, dimx, dimy, dimz, i + 0.5, j + 0.5, k + 0.5, 0.5, crit_point_list );
						if ( err_code == P3D_TRUE ) continue;
						if ( err_code == P3D_MEM_ERROR ) return P3D_MEM_ERROR;
					}
				}
			}
	
	
	// Return success code:
	return P3D_SUCCESS;

MEM_ERROR:

	return P3D_MEM_ERROR;
}


int p3dClassifyCriticalPoints (
    crit_point_list_t crit_point_list,
    float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz	
	)
{
	crit_point_list_t tmp_list;
	crit_point_t crit_point;

	double i,j,k;

	double* in_gvf_x;
	double* in_gvf_y;
	double* in_gvf_z;

	
	double* jacob;
	double* eigvectors;
	double* eigvals_re;
	double* eigvals_im;
	

	double vecLength;
	double vdist = 1.00 / CELL_SUBDIVISION_FACTOR;

	
	P3D_TRY ( jacob = (double*) malloc(9*sizeof(double)));
	P3D_TRY ( eigvectors = (double*) malloc(9*sizeof(double)));
	P3D_TRY ( eigvals_re = (double*) malloc(3*sizeof(double)));
	P3D_TRY ( eigvals_im = (double*) malloc(3*sizeof(double)));
	

	P3D_TRY ( in_gvf_x = (double*) malloc(6*sizeof(double)));
	P3D_TRY ( in_gvf_y = (double*) malloc(6*sizeof(double)));
	P3D_TRY ( in_gvf_z = (double*) malloc(6*sizeof(double)));


	// Scan list:
	tmp_list = crit_point_list;

	while ( crit_point_list_isempty( tmp_list ) == P3D_FALSE )
	{
		// Get current element:
		crit_point = tmp_list->elem;

		i = crit_point.x;
		j = crit_point.y;
		k = crit_point.z;			
						
		//
		// Classify the critical points as atracting/repelling nodes or saddles
		// based on the real and imaginary part of the eigenvalues of the 
		// Jacobian matrix.
		//

		// NOTE:	same array used as input for Jacobian matrix and as output for
		//			eigenvalues
		
		// Interpolate values in the neighborhood:
		in_gvf_x[0] = interpolation( gvf_x, dimx, dimy, dimz, i + vdist, j, k );
		in_gvf_x[1] = interpolation( gvf_x, dimx, dimy, dimz, i - vdist, j, k );
		in_gvf_x[2] = interpolation( gvf_x, dimx, dimy, dimz, i, j + vdist, k );
		in_gvf_x[3] = interpolation( gvf_x, dimx, dimy, dimz, i, j - vdist, k );
		in_gvf_x[4] = interpolation( gvf_x, dimx, dimy, dimz, i, j, k + vdist );
		in_gvf_x[5] = interpolation( gvf_x, dimx, dimy, dimz, i, j, k - vdist );
				
		in_gvf_y[0] = interpolation( gvf_y, dimx, dimy, dimz, i + vdist, j, k );
		in_gvf_y[1] = interpolation( gvf_y, dimx, dimy, dimz, i - vdist, j, k );
		in_gvf_y[2] = interpolation( gvf_y, dimx, dimy, dimz, i, j + vdist, k );
		in_gvf_y[3] = interpolation( gvf_y, dimx, dimy, dimz, i, j - vdist, k );
		in_gvf_y[4] = interpolation( gvf_y, dimx, dimy, dimz, i, j, k + vdist );
		in_gvf_y[5] = interpolation( gvf_y, dimx, dimy, dimz, i, j, k - vdist );

		in_gvf_z[0] = interpolation( gvf_z, dimx, dimy, dimz, i + vdist, j, k );
		in_gvf_z[1] = interpolation( gvf_z, dimx, dimy, dimz, i - vdist, j, k );
		in_gvf_z[2] = interpolation( gvf_z, dimx, dimy, dimz, i, j + vdist, k );
		in_gvf_z[3] = interpolation( gvf_z, dimx, dimy, dimz, i, j - vdist, k );
		in_gvf_z[4] = interpolation( gvf_z, dimx, dimy, dimz, i, j, k + vdist );
		in_gvf_z[5] = interpolation( gvf_z, dimx, dimy, dimz, i, j, k - vdist );

		// Determine Jacobian matrix:
		jacob[0] = (in_gvf_x[0] - in_gvf_x[1]) / (2.0 * vdist);
		jacob[1] = (in_gvf_x[2] - in_gvf_x[3]) / (2.0 * vdist);
		jacob[2] = (in_gvf_x[4] - in_gvf_x[5]) / (2.0 * vdist);
	    
		jacob[3] = (in_gvf_y[0] - in_gvf_y[1]) / (2.0 * vdist);
		jacob[4] = (in_gvf_y[2] - in_gvf_y[3]) / (2.0 * vdist);
		jacob[5] = (in_gvf_y[4] - in_gvf_y[5]) / (2.0 * vdist);
	    
		jacob[6] = (in_gvf_z[0] - in_gvf_z[1]) / (2.0 * vdist);
		jacob[7] = (in_gvf_z[2] - in_gvf_z[3]) / (2.0 * vdist);
		jacob[8] = (in_gvf_z[4] - in_gvf_z[5]) / (2.0 * vdist);

		//
		// Calculate the eigenvalues and eigenvectors of the Jacobian matrix:
		//
		p3dComputeEigenVal( jacob, eigvectors, eigvals_re, eigvals_im, 3 );
	
					
		if(	(eigvals_re[0] < 0.00) && (eigvals_re[1] < 0.00) && (eigvals_re[2] < 0.00) )
		{
			// All real parts of the eigenvalues are negative => attracting node:
			crit_point.type = CPT_ATTRACTING_NODE;
		}
		else 
		{							
			if(	(eigvals_re[0] > 0.00) && (eigvals_re[1] > 0.00) &&	(eigvals_re[2] > 0.00))
			{
				// All real parts of the eigenvalues are positive => repelling node:
				crit_point.type = CPT_REPELLING_NODE;
			}
			else 
			{
				// Two real parts of the eigenvalues are of one sign and the other 
				// has the opposite sign => this is a saddle point
				crit_point.type = CPT_SADDLE;
			}		
		}

		// Set eigenvalues:
		crit_point.eval0 = eigvals_re[0];
		crit_point.eval1 = eigvals_re[1];
		crit_point.eval2 = eigvals_re[2];
	    
		// Set normalized eigenvectors:
		vecLength = sqrt(	eigvectors[0]*eigvectors[0] + 
							eigvectors[3]*eigvectors[3] + 
							eigvectors[6]*eigvectors[6] );
		if ( vecLength == 0.00 ) 
		{
			vecLength = 1.00;
		}

		crit_point.evect0_x = eigvectors[0] / vecLength;
		crit_point.evect0_y = eigvectors[3] / vecLength;
		crit_point.evect0_z = eigvectors[6] / vecLength;
	    


		vecLength = sqrt(	eigvectors[1]*eigvectors[1] + 
							eigvectors[4]*eigvectors[4] + 
							eigvectors[7]*eigvectors[7] );
		if(vecLength == 0.00) 
		{
			vecLength = 1.00;
		}

		crit_point.evect1_x = eigvectors[1] / vecLength;
		crit_point.evect1_y = eigvectors[4] / vecLength;
		crit_point.evect1_z = eigvectors[7] / vecLength;
	    
		vecLength = sqrt(	eigvectors[2]*eigvectors[2] + 
							eigvectors[5]*eigvectors[5] + 
							eigvectors[8]*eigvectors[8] );
		if(vecLength == 0.00) 
		{
			vecLength = 1.00;
		}

		crit_point.evect2_x = eigvectors[2] / vecLength;
		crit_point.evect2_y = eigvectors[5] / vecLength;
		crit_point.evect2_z = eigvectors[8] / vecLength;
	  
		// Set crit_point to current element:
		tmp_list->elem = crit_point;

		// The list will point on the next element:
		tmp_list = tmp_list->next;		
	}

	// Release resources and return success:		

	if ( in_gvf_x != NULL ) free( in_gvf_x );
	if ( in_gvf_y != NULL ) free( in_gvf_y );
	if ( in_gvf_z != NULL ) free( in_gvf_z );

	if ( jacob != NULL ) free( jacob );
	if ( eigvectors != NULL ) free( eigvectors );
	if ( eigvals_re != NULL ) free( eigvals_re );
	if ( eigvals_im != NULL ) free( eigvals_im );

	return P3D_SUCCESS;

MEM_ERROR:

	// Release resources and return error:		

	if ( in_gvf_x != NULL ) free( in_gvf_x );
	if ( in_gvf_y != NULL ) free( in_gvf_y );
	if ( in_gvf_z != NULL ) free( in_gvf_z );

	if ( jacob != NULL ) free( jacob );
	if ( eigvectors != NULL ) free( eigvectors );
	if ( eigvals_re != NULL ) free( eigvals_re );
	if ( eigvals_im != NULL ) free( eigvals_im );

	return P3D_MEM_ERROR;
}