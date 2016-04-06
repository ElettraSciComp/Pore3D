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


