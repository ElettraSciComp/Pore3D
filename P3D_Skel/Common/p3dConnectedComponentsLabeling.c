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

/************************************************************************
 * CCL3 Performs connected component labeling of a 3D binary volume.
 *   CCL3 performs connected component labeling of a 3D binary volume 
 *   using 26-connectivity for object voxels and 6-connectivity for 
 *   background voxels. The algorithm uses a flexible tradeoff between
 *   speed and memory by means of parameter SIZE.
 *
 *   LB = CCL3(VOI) apply the fastest connected component labeling (i.e. 
 *   SIZE = 3).
 *
 *   LB = CCL3(VOI, SIZE) apply connected component labeling using the 
 *   specified SIZE for the sub-volume scanning (see references).
 *
 *   Class Support
 *   -------------      
 *   The input volume VOI must be logical. SIZE must be an odd value greater
 *   than 3.
 *
 *   Remarks
 *   -------
 *   The SIZE parameters controls the tradeoff between computational time
 *   and memory occupation. The fastest execution requires SIZE = 3.
 *
 *   Example
 *   ------------
 *
 *   lbl = ccl3(voi);
 *
 *   References
 *   ----------
 *   [1] Q. Hu et al. Fast connected-component labeling in
 *   three-dimensional binary images based on iterative recursion.
 *   Computer Vision and Image Understanding, 99(2005):414-434, 2005.


 *   Author: Brun Francesco - Universit� degli studi di Trieste
 *     Date: 10 may 2008

 ************************************************************************/
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#include "p3dCoordsQueue.h"
#include "p3dUIntList.h"
#include "p3dBoundingBoxList.h"
#include "p3dUtils.h"


#define ON_LABEL 1     // ON Value:
#define FIRST_LABEL 3  // First label:
#define TEMP_LABEL 2   // Temporary label:


void putCoordsInQueue (
	unsigned short* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int i,
	const int j,
	const int k,
	const int ct,
	const int win_size,
	unsigned short m,
	coords_queue_t* queue				// FIFO data structure for coords storage
	)
{
	coords_t tmp_coords;		// Temporary coords


	// If we're on an object voxel:
	if ( out_rev[ I(i,j,k,dimx,dimy) ] == ON_LABEL )
	{
		// Check if this is the last step:
		if ( ct == (win_size - 2))
		{
			// Border voxels are set to m:
			out_rev[ I(i,j,k,dimx,dimy) ] = m;		

			// A border voxel is pushed into queue:
			tmp_coords.x = i;
			tmp_coords.y = j;
			tmp_coords.z = k;
			coords_queue_push(queue,tmp_coords);
		}
		else
		{
			// Mark voxel for further scan:
			out_rev[ I(i,j,k,dimx,dimy) ] = TEMP_LABEL;
		}
	}
}



void conn_fun6 ( 
	unsigned short* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int i,
	const int j,
	const int k,
	const int ct,
	const int win_size,
	unsigned short m,
	coords_queue_t* queue				// FIFO data structure for coords storage
	)
{
	int a,b,c;   

	c = k - 1; b = j; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j - 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j + 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		
	
}


void conn_fun18 ( 
	unsigned short* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int i,
	const int j,
	const int k,
	const int ct,
	const int win_size,
	unsigned short m,
	coords_queue_t* queue				// FIFO data structure for coords storage
	)
{
	int a,b,c;   

	// Perform 6-connectivity:
	conn_fun6 (out_rev, dimx, dimy, dimz, i, j, k, ct, win_size, m, queue );

	// Do other 12 tests:

	// On k-1 plane:
	c = k - 1; b = j - 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j + 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j + 1; a = i;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		


	// On k plane:
	c = k; b = j - 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j - 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j + 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k; b = j + 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );	

	
}


void conn_fun26 ( 
	unsigned short* out_rev,
	const int dimx,
	const int dimy,
	const int dimz,
	const int i,
	const int j,
	const int k,
	const int ct,
	const int win_size,
	unsigned short m,
	coords_queue_t* queue				// FIFO data structure for coords storage
	)
{
	int a,b,c;  

	// Perform 18-connectivity:
	conn_fun18 (out_rev, dimx, dimy, dimz, i, j, k, ct, win_size, m, queue );

	// Do other 8 tests:
	
	// On k-1 plane:
	c = k - 1; b = j - 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j - 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j + 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k - 1; b = j + 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );	


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j - 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j + 1; a = i - 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );		

	c = k + 1; b = j + 1; a = i + 1;
	putCoordsInQueue( out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue );	

}


int first_unlabeled (
	unsigned short* lbl, 
	const int nRows, 
	const int nCols, 
	const int nPlanes, 
	const int winsize, 
	coords_t* coords
	)
{	
	int i,j,k;
	int off;

	off = winsize / 2; // integer division:

	// Volume scanning:
    for( k = off; k < (nPlanes - off); k++ )  
        for( j = off; j < (nCols - off); j++ )
			for( i = off; i < (nRows - off); i++ )
			{
				if (lbl[ I(i,j,k,nRows,nCols) ] == ON_LABEL)
				{
					coords->x = i;
					coords->y = j;
					coords->z = k;
					return P3D_TRUE;
				}
			}
	
	return P3D_FALSE;
}

int p3dConnectedComponentsLabeling (
	 unsigned char* in_rev,
	 unsigned short* out_rev,	 
	 int* numOfConnectedComponents,	// OUT: dim of arrays
	 unsigned int** volumes,        // OUT: array of sizes (voxel counting)
	 bb_t** boundingBoxes,          // OUT: array of bounding boxes
	 const int dimx,
	 const int dimy, 
	 const int dimz,	
	 const int conn,
	 const int skip_borders
	 )
{    
	// Size of local window (3 in current implementation).
	// Modify this value with 5,7,... if memory problem
	// occurs. Uncomment also some code (see further).
	int win_size;
		
	coords_queue_t queue;		// FIFO data structure for coords storage
	coords_t coords;			// Current coords

	uint_list_t vol_list;		// List for region volume counting
	bb_list_t bb_list;			// List for bounding boxes
	bb_t curr_bb;

	unsigned short m;			// Current label
	int ct, offset;   	// Counter for local subvolume size
	int a,b,c;
	int i,j,k;
	int m_ct;			// Counter for volume counting
	int lbl_ct;

	// Dimensions for padded/cropped volumes:
	int a_dimx, a_dimy, a_dimz;
	int a_rad;

	// Padded and cropped temporary input and output:
	unsigned char* tmp_in_rev;
	unsigned short* tmp_out_rev;
 
	// Pointer to a function for kind of connectivity:
	void (*conn_fun) ( 
		unsigned short* out_rev,
		const int dimx,
		const int dimy,
		const int dimz,
		const int i,
		const int j,
		const int k,
		const int ct,
		const int win_size,
		unsigned short m,
		coords_queue_t* queue				// FIFO data structure for coords storage
		);

	// The win_size parameters controls the tradeoff between computational time
    // and memory occupation. The fastest execution requires win_size = 3.
    win_size = 3;

	// Set the correct type of connectivity based on input parameter:
	if ( conn == CONN6 )
		conn_fun = conn_fun6;
	else if ( conn == CONN18 )
		conn_fun = conn_fun18;
	else // default:
		conn_fun = conn_fun26;    

    // Apply algorithm:   

	// Create temporary input replicate-padded:
	a_rad = 1;

	// Compute dimensions of padded REV:
	a_dimx = dimx + a_rad*2;
	a_dimy = dimy + a_rad*2;
	a_dimz = dimz + a_rad*2;

	// Initialize input:
	P3D_MEM_TRY ( tmp_in_rev = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) ) );
	p3dZeroPadding3D_uchar2uchar ( in_rev, tmp_in_rev, dimx, dimy, dimz, a_rad);


	// Initialize output label volume with ON_LABEL on non-zero values
	// of input volume:
	P3D_MEM_TRY ( tmp_out_rev = (unsigned short*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned short) ) );
	
	#pragma omp parallel for
	for (ct = 0; ct < (a_dimx*a_dimy*a_dimz); ct++)
		tmp_out_rev[ct] = ( tmp_in_rev[ct] ) ? ON_LABEL : 0;			

	

	// Initialize variables:
	coords_queue_init(&queue);
	uint_list_init(&vol_list);
	bb_list_init(&bb_list);
	

	m = FIRST_LABEL;
	lbl_ct = 0;
	

	// While not all voxels are labeled:
	while( first_unlabeled(tmp_out_rev, a_dimx, a_dimy, a_dimz, win_size, &coords) == P3D_TRUE )
	{		
		// Re-init counter for current connected component:
		m_ct = 0;

		// Init bounding box:
		curr_bb.max_x = 0;
		curr_bb.max_y = 0;
		curr_bb.max_z = 0;
		curr_bb.min_x = INT_MAX;
		curr_bb.min_y = INT_MAX;
		curr_bb.min_z = INT_MAX;

		// Push the first unlabeled object voxel into queue:
		coords_queue_push(&queue, coords);

		// While the queue is not empty:
		while( coords_queue_isempty(queue) == P3D_FALSE )
		{
			// Pop the first element of the queue:
			coords = coords_queue_pop(&queue);
			
			// Mark the extracted element:
			tmp_out_rev[ I(coords.x,coords.y,coords.z,a_dimx,a_dimy) ] = TEMP_LABEL;			

			// Perform sub-volume iterative scanning:
			ct = 1;
			while (ct < win_size)
			{
				offset = ct / 2; // integer division

				for (k = (coords.z - offset); k <= (coords.z + offset); k++)
					for (j = (coords.y - offset); j <= (coords.y + offset); j++)
						for (i = (coords.x - offset); i <= (coords.x + offset); i++)
						{	
					
							if (tmp_out_rev[ I(i,j,k,a_dimx,a_dimy) ] == TEMP_LABEL)
							{
								// Perform 6-, 18-, or 26-connectivity:
								conn_fun ( tmp_out_rev, a_dimx, a_dimy, a_dimx, i,
									j, k, ct, win_size, m, &queue );
								
								// A voxels previously marked is set to m:
								tmp_out_rev[ I(i,j,k,a_dimx,a_dimy) ] = m;	

								// Define bounding box:
								curr_bb.max_x = MAX( curr_bb.max_x, i - a_rad );
								curr_bb.max_y = MAX( curr_bb.max_y, j - a_rad );
								curr_bb.max_z = MAX( curr_bb.max_z, k - a_rad );
								curr_bb.min_x = MIN( curr_bb.min_x, i - a_rad );
								curr_bb.min_y = MIN( curr_bb.min_y, j - a_rad );
								curr_bb.min_z = MIN( curr_bb.min_z, k - a_rad );
								
								// Increment counter for region volume:
								m_ct++;
							}
						}
				// Increase local sub-volume window:
				ct = ct + 2;
			}
		}


		// Skip object if required:
		if ( ( skip_borders == P3D_TRUE ) && ( (curr_bb.min_x == 0) || (curr_bb.min_y == 0) 
			|| (curr_bb.min_z == 0) || (curr_bb.max_x == (dimx - 1))
			|| (curr_bb.max_y == (dimy - 1)) ||  (curr_bb.max_z == (dimz - 1)) ) )
		{
			// Restore OBJECT label instead of current label:
			#pragma omp parallel for private(a,b)
			for ( c = curr_bb.min_z; c <= curr_bb.max_z; c++ )
				for ( b = curr_bb.min_y; b <= curr_bb.max_y; b++ )
					for ( a = curr_bb.min_x; a <= curr_bb.max_x; a++ )
					{
						if( tmp_out_rev[ I(a+a_rad,b+a_rad,c+a_rad,a_dimx,a_dimy) ] == m )
						{ 
							tmp_out_rev[ I(a+a_rad,b+a_rad,c+a_rad,a_dimx,a_dimy) ] = USHRT_MAX; // Object
						}
					}	
		}
		else
		{
			// Save connected component volume:
			uint_list_add( &vol_list, m_ct );
			bb_list_add( &bb_list, curr_bb );

					
			// Increment counter (it is equal to (m - 3) if labels are NOT randomly assigned):
			m++;
			lbl_ct++;
		}
	} 

	// Return number of connected components labeled:
	if ( numOfConnectedComponents != NULL )
		*numOfConnectedComponents = lbl_ct;

	
	// Return the array of volumes for each connected component.
	// (The list is deleted after this operation):
	if ( volumes != NULL )
		(*volumes) = uint_list_toarray ( &vol_list, lbl_ct );
	else
		uint_list_clear ( &vol_list );

	if ( boundingBoxes != NULL )
		(*boundingBoxes) = bb_list_toarray ( &bb_list, lbl_ct );
	else
		bb_list_clear ( &bb_list );

	
	// Crop output:
	p3dCrop3D_ushort2ushort ( tmp_out_rev, out_rev, a_dimx, a_dimy, a_dimz, a_rad );



	// Release resources:
	if ( tmp_in_rev != NULL ) free(tmp_in_rev);
	if ( tmp_out_rev != NULL ) free(tmp_out_rev);

	// Return OK:
	return P3D_SUCCESS;

MEM_ERROR: 

	// Release resources:
	if (tmp_in_rev != NULL) free(tmp_in_rev);
	if (tmp_out_rev != NULL) free(tmp_in_rev);

	// Return error code:
	return P3D_ERROR;
}



int p3dGetMaxVolumeRegion (   
	unsigned char* in_rev, 
	unsigned char* out_rev, 
	const int dimx,
	const int dimy, 
	const int dimz,
	int conn
	)
{	
	unsigned short* lbl_rev;

	unsigned int* volumes;
	unsigned short lbl;

	int i,j,k, num_el;

	int lbl_max;
	unsigned int vol_max;

	int err_code;

	// Initialize output by cloning input:
	memcpy(out_rev, in_rev, dimx*dimy*dimz*sizeof(unsigned char));

	// Allocate memory for labels REV:
	lbl_rev = (unsigned short*) malloc (dimx*dimy*dimz*sizeof(unsigned short));

	if ( lbl_rev == NULL )
	{
		// Return error code and exit:
		return P3D_ERROR;
	}

	// Perform connected component labeling:
	err_code = p3dConnectedComponentsLabeling ( in_rev, lbl_rev, &num_el, &volumes, NULL, dimx, 
		dimy, dimz, conn, P3D_FALSE);

	if ( err_code == P3D_ERROR )
	{		
		// Free memory if previous malloc were successfully:
		if ( lbl_rev != NULL ) free ( lbl_rev );
		if ( volumes != NULL ) free ( volumes );

		// Return error code and exit:
		return P3D_ERROR;
	}

	
	lbl_max = 0;
	vol_max = 0;

	// For each connected component labeled:
	for ( lbl = 0; lbl < num_el; lbl++ )
	{
		if (volumes[lbl] > vol_max)
		{
			lbl_max = lbl;
			vol_max = volumes[lbl];
		}
	}

	// Remove connected component:
	for( k = 0; k < dimz; k++ )
		for( j = 0; j < dimy; j++ )
			for( i = 0; i < dimx; i++ )            
			{  
				// Labels start from 3:
				if (lbl_rev[ I(i,j,k,dimx,dimy) ] != (lbl_max + 3) )
				{
					out_rev[ I(i,j,k,dimx,dimy) ] = 0;
				}
			}

	// Free memory:
	if ( lbl_rev != NULL ) free ( lbl_rev );
	if ( volumes != NULL ) free ( volumes );

	// Return OK:
	return P3D_SUCCESS;
}
