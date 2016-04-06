#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "p3dGVF.h"

#include "../Common/p3dUtils.h"


//
// In-line gradient computation:
//
double grad_x ( unsigned char* in_im, 
			   const int dimx, const int dimy, const int dimz,
			   const int i,    const int j,    const int k )
{
	return (0.5 * ( in_im[ I(i+1,j,k,dimx,dimy) ] - in_im[ I(i-1,j,k,dimx,dimy) ] ));
}
double grad_y ( unsigned char* in_im, 
			   const int dimx, const int dimy, const int dimz,
			   const int i,    const int j,    const int k )
{
	return (0.5 * ( in_im[ I(i,j+1,k,dimx,dimy) ] - in_im[ I(i,j-1,k,dimx,dimy) ] ));
}
double grad_z ( unsigned char* in_im, 
			   const int dimx, const int dimy, const int dimz,
			   const int i,    const int j,    const int k )
{
	return (0.5 * ( in_im[ I(i,j,k+1,dimx,dimy) ] - in_im[ I(i,j,k-1,dimx,dimy) ] )); 
}
//
// In-line gradient magnitude computation:
//
double grad_mag ( unsigned char* in_im, 
				 const int dimx, const int dimy, const int dimz,
				 const int i,    const int j,    const int k )
{
	return ( grad_x(in_im,dimx,dimy,dimz,i,j,k)*grad_x(in_im,dimx,dimy,dimz,i,j,k) +
			 grad_y(in_im,dimx,dimy,dimz,i,j,k)*grad_y(in_im,dimx,dimy,dimz,i,j,k) +
			 grad_z(in_im,dimx,dimy,dimz,i,j,k)*grad_z(in_im,dimx,dimy,dimz,i,j,k) );		
}

//
// In-line laplacian computation:
//
double lapl (  float* in_im, 
			   const int dimx, const int dimy, const int dimz,
			   const int i,    const int j,    const int k)
{
	return ( (1.0/6.0) * ( in_im[ I(i-1,j,k,dimx,dimy) ] + 
					in_im[ I(i+1,j,k,dimx,dimy) ] + 
					in_im[ I(i,j-1,k,dimx,dimy) ] + 
					in_im[ I(i,j+1,k,dimx,dimy) ] + 
					in_im[ I(i,j,k-1,dimx,dimy) ] + 
					in_im[ I(i,j,k+1,dimx,dimy) ] ) - 
					in_im[ I(i,j,k,dimx,dimy) ] ); 
}


int p3dGVF (	
    unsigned char* in_im,
	float* out_x,
	float* out_y,
	float* out_z,
	const int dimx,	
	const int dimy,	
	const int dimz,	
	const double  mu,
	const double  eps
	)
{
	int rad = 1;
	int ct  = 0;

	int i,j,k;

	double r;
	float  tmp_x, tmp_y, tmp_z;
	double eps_x, eps_y, eps_z;
	double loc_eps_x, loc_eps_y, loc_eps_z;
	int    stable_x, stable_y, stable_z;

	
	//
	// PHASE 1: Compute gradient (central difference used).
	//
		   
	#pragma omp parallel for private(i, j)
	for( k = rad; k < (dimz - rad); k++ )    
		for( j = rad; j < (dimy - rad); j++ )
			for( i = rad; i < (dimx - rad); i++ )
			{		
				if ( in_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND )
				{
					out_x[ I(i,j,k,dimx,dimy) ] = (float) (0.5 * ( (float) (in_im[ I(i+1,j,k,dimx,dimy) ])
						- (float) (in_im[ I(i-1,j,k,dimx,dimy) ]) ));
					out_y[ I(i,j,k,dimx,dimy) ] = (float) (0.5 * ( (float) (in_im[ I(i,j+1,k,dimx,dimy) ])
						- (float) (in_im[ I(i,j-1,k,dimx,dimy) ]) )); 	 
					out_z[ I(i,j,k,dimx,dimy) ] = (float) (0.5 * ( (float) (in_im[ I(i,j,k+1,dimx,dimy) ])
						- (float) (in_im[ I(i,j,k-1,dimx,dimy) ]) )); 	
				}
			}

	//
	// PHASE 2: Iteratively solve GVF:
	//	

	
	// Solve GVF iteratively:
	//
	//   gvf_x = gvf_x + mu*4*lapl( gvf_x ) - SqrMagf.*( gvf_x - grad_x );
    //   gvf_y = gvf_y + mu*4*lapl( gvf_y ) - SqrMagf.*( gvf_x - grad_y );
    //   gvf_z = gvf_z + mu*4*lapl( gvf_z ) - SqrMagf.*( gvf_z - grad_z );
	//
    //   mu = 0.15 with 500 iterations (Hassouna, 2009)
	
	stable_x = P3D_FALSE;
	stable_y = P3D_FALSE;
	stable_z = P3D_FALSE;


	while ( (stable_x == P3D_FALSE) || (stable_y == P3D_FALSE) || (stable_z == P3D_FALSE) )
	{    
		eps_x = 0.0;
		eps_y = 0.0;
		eps_z = 0.0;

		#pragma omp parallel for private(i, j, tmp_x, tmp_y, tmp_z, loc_eps_x, loc_eps_y, loc_eps_z)		
		for( k = rad; k < (dimz - rad); k++ )
		{
			loc_eps_x = 0.0;
			loc_eps_y = 0.0;
			loc_eps_z = 0.0;

			for( j = rad; j < (dimy - rad); j++ )
				for( i = rad; i < (dimx - rad); i++ )
				{
					if (in_im[ I(i,j,k,dimx,dimy) ] != BACKGROUND )
					{
						tmp_x = (float) (out_x[ I(i,j,k,dimx,dimy) ] +
							mu*4.0*lapl(out_x,dimx,dimy,dimz,i,j,k) - grad_mag(in_im,dimx,dimy,dimz,i,j,k) *
							(out_x[ I(i,j,k,dimx,dimy) ] - grad_x(in_im,dimx,dimy,dimz,i,j,k) ));

						tmp_y = (float) (out_y[ I(i,j,k,dimx,dimy) ] +
							mu*4.0*lapl(out_y,dimx,dimy,dimz,i,j,k) - grad_mag(in_im,dimx,dimy,dimz,i,j,k) *
							(out_y[ I(i,j,k,dimx,dimy) ] - grad_y(in_im,dimx,dimy,dimz,i,j,k) ));

						tmp_z = (float) (out_z[ I(i,j,k,dimx,dimy) ] +
							mu*4.0*lapl(out_z,dimx,dimy,dimz,i,j,k) - grad_mag(in_im,dimx,dimy,dimz,i,j,k) *
							(out_z[ I(i,j,k,dimx,dimy) ] - grad_z(in_im,dimx,dimy,dimz,i,j,k) ));		

						// Compute differences:
						loc_eps_x = MAX(loc_eps_x,fabs(tmp_x - out_x[ I(i,j,k,dimx,dimy) ]));
						loc_eps_y = MAX(loc_eps_y,fabs(tmp_y - out_y[ I(i,j,k,dimx,dimy) ]));
						loc_eps_z = MAX(loc_eps_z,fabs(tmp_z - out_z[ I(i,j,k,dimx,dimy) ]));

						// Assign values for this step:
						out_x[ I(i,j,k,dimx,dimy) ] = tmp_x;
						out_y[ I(i,j,k,dimx,dimy) ] = tmp_y;
						out_z[ I(i,j,k,dimx,dimy) ] = tmp_z;
					}
				}
			
			#pragma omp critical
			{
				eps_x = MAX(eps_x,loc_eps_x); 
				eps_y = MAX(eps_y,loc_eps_y); 
				eps_z = MAX(eps_z,loc_eps_z); 
			}
		}

		// Check for stability:		
		
		// If at least one value has a difference from the 
		// previous step which is greater than EPS continue 
		// the iterations:
		stable_x = (eps_x > eps) ? P3D_FALSE : P3D_TRUE;
		stable_y = (eps_y > eps) ? P3D_FALSE : P3D_TRUE;
		stable_z = (eps_z > eps) ? P3D_FALSE : P3D_TRUE;

		ct++;
	}


	//
	// Normalize vectors flow:
	//

	#pragma omp parallel for private(r)
	for ( i = 0; i < (dimx*dimy*dimz); i++)
	{
		if( in_im[ i ] != BACKGROUND ) 
		{	  
			r = out_x[ i ]*out_x[ i ] + 
				out_y[ i ]*out_y[ i ] + 
				out_z[ i ]*out_z[ i ];
		    
			if( r > 0.00 ) 
			{
				r = sqrt(r);
		    
				out_x[ i ] = (float) ( out_x[ i ] / r );
				out_y[ i ] = (float) ( out_y[ i ] / r );
				out_z[ i ] = (float) ( out_z[ i ] / r );
			}
		}
	}

	// Return the number of iterations:
	return ct;
}

