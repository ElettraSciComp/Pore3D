
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


 *   Author: Brun Francesco 
 *     Date: 10 may 2008

 ************************************************************************/
#include <limits.h>

#include "p3dBoundingBoxT.h"
#include "p3dUIntList.h"

int p3dConnectedComponentsLabeling_ushort (
	 unsigned char* in_rev,
	 unsigned short* out_rev,	 
	 unsigned int* numOfConnectedComponents,	// OUT: dim of arrays
	 unsigned int** volumes,	    // OUT: array of sizes (voxel counting)
	 bb_t** boundingBoxes,          // OUT: array of bounding boxes
	 const int dimx,
	 const int dimy, 
	 const int dimz,	
	 const int conn,
     const int random_lbl,
	 const int skip_borders
	 );

int p3dConnectedComponentsLabeling_uint (
	 unsigned char* in_rev,
	 unsigned int* out_rev,	 
	 unsigned int* numOfConnectedComponents,	// OUT: dim of arrays
	 unsigned int** volumes,	// OUT: array of sizes (voxel counting)
	 bb_t** boundingBoxes,          // OUT: array of bounding boxes
	 const int dimx,
	 const int dimy, 
	 const int dimz,	
	 const int conn,
     const int random_lbl,
	 const int skip_borders
	 );


