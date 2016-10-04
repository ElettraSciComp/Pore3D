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

/* 
* p3dSquaredEuclideanDistanceTransform  
* 
* Computes the squared Euclidean distance transform of the input volume. For 
* each non-zero voxel in input image, the distance transform assigns a number 
* that is the squared distance between that voxel and the nearest zero pixel 
* (background). The squared Euclidean distance transform is a non-approximated
* Euclidean distance metric.
*
* Remarks
* -------
* Computational cost is O(N^3) for a NxNxN input volume. Memory 
* requirement is two internal volumes of a NxNxN and two linear arrays of N 
* elements, i.e. O(N^3) of type int. Current implementation uses int elements,
* therefore (for most compilers) memory requirement is 8 x O(N^3) the size of
* input volume.
*
* References
* ----------
* [1] T. Hirata. "A unified linear-time algorithm for computing distance 
* maps", Information Processing Letters, 58(3):129-133, May 1996.
*
* [2] A. Meijster, J.B.T.M. Roerdink and W. H. Hesselink. "A general 
* algorithm for computing distance transforms in linear time", 
* Mathematical Morphology and its Applications to Image and Signal 
* Processing, pp. 331-340. Kluwer, 2000.
*

* 
* Copyright 2008, SYRMEP Group - Sincrotrone Trieste S.C.p.A.
*
* Author:	Brun Francesco
* Version:	1.0
* Date:		19 june 2008
*
*/
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#include "p3dSquaredEuclideanDT.h"
#include "p3dUtils.h"

// Set a signed integer type able to get over 65.535 and the 
// maximum value representable (INF) with this specified type 
// (see <limits.h>). These settings depends on the specific 
// compilator used. 
typedef int int_type;	
#define INF INT_MAX	

/* -------------------
 * Mathematical operators with INF management.
 * ------------------- */

/** 
 **************************************************
 * @b sum
 * @param a int_type number with INF
 * @param b int_type number with INF
 * @return The sum of a and b handling INF
 **************************************************/
int_type sum( int_type a, int_type b ) 
{
    if ((a == INF) || (b == INF))     
        return INF;    
    else 
        return a + b;
}

/** 
 **************************************************
 * @b prod
 * @param a long number with INF
 * @param b long number with INF
 * @return The product of a and b handling INF
 **************************************************/
int_type prod( int_type a, int_type b ) 
{
    if ((a == INF) || (b == INF)) 
        return INF;  
    else 
        return a * b;
}
/** 
 **************************************************
 * @b opp
 * @param a long number with INF
 * @return The opposite of a  handling INF
 **************************************************/
int_type opp( int_type a ) 
{
    if (a == INF) 
        return INF;
    else 
        return -a;
}

/** 
 **************************************************
 * @b intdivint
 * @param divid long number with INF
 * @param divis long number with INF
 * @return The division (integer) of divid out of divis handling INF
 **************************************************/
int_type intdivint( int_type divid, int_type divis ) 
{
    if (divis == 0) 
        return  INF;
    if (divid == INF) 
        return  INF;
    else 
        return  divid / divis;
}

/* -------------------
 * Mathematical operators for the parabola definition.
 * ------------------- */

/** 
 **************************************************
 * @b f
 * @param x 
 * @param i 
 * @param gi2 
 * @return Definition of a parabola
 **************************************************/
int_type f( int x, int i, int_type gi2 )
{
    return sum((x-i)*(x-i), gi2);
}

/** 
 **************************************************
 * @b sep
 * @param i 
 * @param u 
 * @param gi2 
 * @param gu2 
 * @return The abscissa of the intersection point between two parabolas
 **************************************************/
int_type sep( int i, int u, int_type gi2, int_type gu2 ) 
{
    return intdivint(sum( sum((long) (u*u - i*i),gu2), opp(gi2) ), 2*(u-i));
}



/* -------------------
 * Steps of distance transform algorithm.
 * ------------------- */


/** 
 **************************************************
 * @b stepX
 * @param voi input volume
 * @param sdt_x DT along the x-direction
 **************************************************/
int stepX (	const unsigned char* voi, 
			int_type* sdt_x, 
			const int NRows, 
			const int NCols, 
			const int NPlanes
			)
{
    int x,y,z;
    
    for (z = 0; z < NPlanes; z++) 	    
        for (y = 0; y < NCols ; y++) 
        {
        	if (voi[ I(0,y,z,NRows,NCols) ] == 0)
                sdt_x[ I(0,y,z,NRows,NCols) ] = 0;
        	else 	    
                sdt_x[ I(0,y,z,NRows,NCols) ] = INF;  
	  
        	// Forward scan
        	for (x = 1; x < NRows; x++) 	    
                if (voi[ I(x,y,z,NRows,NCols) ] == 0)
            	    sdt_x[ I(x,y,z,NRows,NCols) ] = 0;
            	else 
            	    sdt_x[ I(x,y,z,NRows,NCols) ]  = sum(1, sdt_x[ I(x-1,y,z,NRows,NCols)] );	  
	
        	//Backward scan
        	for (x = NRows-2; x >= 0; x--)      
                if (sdt_x[ I(x+1,y,z,NRows,NCols) ] < sdt_x[ I(x,y,z,NRows,NCols) ]) 
            	    sdt_x[ I(x,y,z,NRows,NCols) ] = sum(1, sdt_x[ I(x+1,y,z,NRows,NCols) ]);
        }

	// Return success code:
	return P3D_SUCCESS;
}

/** 
 **************************************************
 * @b stepY
 * @param sdt_x the DT along the x-direction
 * @param sdt_xy the DT in the xy-slices
 **************************************************/
int stepY (	const int_type* sdt_x, 
			int_type* sdt_xy, 
			const int NRows, 
			const int NCols, 
			const int NPlanes
			)
{
  
    int* s;  // Center of the upper envelope parabolas
    int* t;  // Separating index between 2 upper envelope 
                            // parabolas 
    int q; 
    int w;
    
    int x,z,u;

    
    P3D_TRY ( s = (int*) calloc(NCols, sizeof(int)) );
    P3D_TRY ( t = (int*) calloc(NCols, sizeof(int)) );
	

    for (z = 0; z < NPlanes; z++) 	    
        for (x = 0; x < NRows; x++) 
        {
            q = 0;    
            s[0] = 0;
            t[0] = 0;
	
            //Forward scanning:
            for (u=1; u < NCols ; u++) 
            {
                while ((q >= 0) &&
                    ( f(t[q],s[q],prod(sdt_x[ I(x,s[q],z,NRows,NCols) ],sdt_x[ I(x,s[q],z,NRows,NCols) ])) > 
                    f(t[q],u,prod(sdt_x[ I(x,u,z,NRows,NCols) ],sdt_x[ I(x,u,z,NRows,NCols) ]))) ) 
                    q--;
	    
                    if (q<0) 
                    {
                        q=0;
                        s[0]=u;
                    }
                else 
                {
                    w = 1 + sep(s[q],u,
                    prod(sdt_x[ I(x,s[q],z,NRows,NCols) ],sdt_x[ I(x,s[q],z,NRows,NCols) ]),
                    prod(sdt_x[ I(x,u,z,NRows,NCols) ],sdt_x[ I(x,u,z,NRows,NCols) ]));

                    if (w < NCols) 
                    {
                        q++;
                        s[q]=u;
                        t[q]=w;
                    }
                }
            }

            //Backward scanning:
            for (u = NCols-1; u >= 0; --u) 
            {
                sdt_xy[ I(x,u,z,NRows,NCols) ] = f(u,s[q],prod(sdt_x[ I(x,s[q],z,NRows,NCols) ],sdt_x[ I(x,s[q],z,NRows,NCols) ]));	      
                if (u==t[q]) 
                    q--;
            }
        }
    
	// Free memory:
    if ( s != NULL ) free(s);
    if ( t != NULL ) free(t);

	// Return success code:
	return P3D_SUCCESS;

MEM_ERROR:

    // Free memory:
    if ( s != NULL ) free(s);
    if ( t != NULL ) free(t);

	// Return error code:
	return P3D_MEM_ERROR;
}

/** 
 **************************************************
 * @b stepZ
 * @param sdt_xy the DT in the xy-slices
 * @param sdt_xyz the final DT
 **************************************************/
int stepZ (	const int_type* sdt_xy, 
			int_type* sdt_xyz, 
			const int NRows, 
			const int NCols, 
			const int NPlanes
			)
{
  
    int* s;  // Center of the upper envelope parabolas
    int* t;  // Separating index between 2 upper envelope 
                            // parabolas 
    int q; 
    int w;
    
    int x,y,u;

    
    P3D_TRY ( s = (int*) calloc(NPlanes, sizeof(int)) );
    P3D_TRY ( t = (int*) calloc(NPlanes, sizeof(int)) );


    for (y = 0; y < NCols; y++) 	    
        for (x = 0; x < NRows; x++) 
        {
            q = 0;
            s[0] = 0;
            t[0] = 0;
	
            //Forward scanning:
            for ( u = 1; u < NPlanes ; u++ ) 
            {
                while ((q >= 0) &&
                    (f(t[q],s[q], sdt_xy[ I(x,y,s[q],NRows,NCols) ]) > 
                    f(t[q],u,sdt_xy[ I(x,y,u,NRows,NCols) ])) ) 
                    q--;
	    
                if (q < 0) 
                {
                    q = 0;
                    s[0] = u;
                }
                else 
                {
                    w = 1 + sep(s[q], u,
                    sdt_xy[ I(x,y,s[q],NRows,NCols) ],
                    sdt_xy[ I(x,y,u,NRows,NCols) ]);
	
                    if (w < NPlanes) 
                    {
                        q++;
                        s[q]=u;
                        t[q]=w;
                    }
                }
            }

            //Backward scanning:
            for (u = NPlanes-1; u >= 0; --u) 
            {
                sdt_xyz[ I(x,y,u,NRows,NCols) ] = f(u,s[q],sdt_xy[ I(x,y,s[q],NRows,NCols) ]);	      
                if ( u == t[q] ) 
                    q--;
            }
        }
    
	// Free memory:
    if ( s != NULL ) free(s);
    if ( t != NULL ) free(t);

	// Return success code:
	return P3D_SUCCESS;

MEM_ERROR:

    // Free memory:
    if ( s != NULL ) free(s);
    if ( t != NULL ) free(t);

	// Return error code:
	return P3D_MEM_ERROR;
}


/************************************************************************
 *  DT3 - Entry point:
 *
 *
 ************************************************************************/
int p3dSquaredEuclideanDT(   
	  unsigned char* in_rev,
	  unsigned short* out_rev,
	  const int dimx,
	  const int dimy, 
	  const int dimz
	  )
{   
	// Temporary values and volume:
	int_type* tmp_rev;
	int_type* tmp_out_rev;

	// Counter:
	int i;


	// Allocate temporary output:
	P3D_TRY ( tmp_out_rev = (int_type*) malloc(dimx*dimy*dimz*sizeof(int_type)) );
	P3D_TRY ( tmp_rev = (int_type*) malloc( dimx*dimy*dimz*sizeof(int_type) ) );

    // Compute transform:
    P3D_TRY ( stepX ( in_rev, tmp_out_rev, dimx, dimy, dimz ) );
    P3D_TRY ( stepY ( tmp_out_rev, tmp_rev, dimx, dimy, dimz ) );
    P3D_TRY ( stepZ ( tmp_rev, tmp_out_rev, dimx, dimy, dimz ) ); 

	// Cast output from int_type to unsigned short:
	for ( i = 0; i < (dimx*dimy*dimz); i++)
		out_rev[i] = (tmp_out_rev[i] > USHRT_MAX) ? USHRT_MAX : (unsigned short) ( tmp_out_rev[i] );

    
    // Free allocated memory:
    if ( tmp_rev != NULL ) free(tmp_rev);
	if ( tmp_out_rev != NULL ) free(tmp_out_rev);

	// Return OK:
	return P3D_SUCCESS;

MEM_ERROR:	

	// Free allocated memory:
    if ( tmp_rev != NULL ) free(tmp_rev);
	if ( tmp_out_rev != NULL ) free(tmp_out_rev);

	// Return error:
	return P3D_MEM_ERROR;
}

