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
* requirement is three extra volumes of a NxNxN and two linear arrays of N 
* elements, i.e. O(N^3). Current implementation uses unsigned short elements
* (see compiler specifications).
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



/************************************************************************
 *  DT3 - Entry point:
 *
 *
 ************************************************************************/
int p3dSquaredEuclideanDT(   
	  unsigned char* in_im,
	  unsigned short* out_im,
	  const int dimx,
	  const int dimy, 
	  const int dimz
	  );
