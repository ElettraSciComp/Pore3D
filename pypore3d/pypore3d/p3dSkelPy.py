import pypore3d.p3d_common_lib
from pypore3d.p3d_common_lib import *

import pypore3d.p3dSkel
from pypore3d.p3dSkel import *
 



################## Skeletonization functions
    
def py_p3dSkeletonLabeling(image_data, dimx, dimy, dimz = 0):
	"""
 
  Scan the input skeleton assigning a different label to topological relevant voxels. On the assumption that input skeleton has  the thinness property, i.e. is one-voxel wide, a labeling in terms of nodes and branches is possible. A more detailed labeling  of skeleton branches is also possible identifying node-to-node, node-to-end and end-to-end branches. 
  
  Syntax:
  ------
  Result = py_p3dSkeletonLabeling ( image_data, dimx, dimy [, dimz = value])
  
  Return Value:
  ------------
  Returns an image of type BYTE with the same dimensions of input image in which each skeleton element is assigned a different  gray-level according to the following codes:
 Node voxel
 
 1. A node voxel has more than two voxels in its neighborhood (End voxel)
 2. An end voxel has exactly one voxel in its neighborhood. (Isolated voxel)
 3. An isolated voxel has no voxels in its neighborhood. In principle, an isolated voxel occurs in the skeletonization of a  perfect sphere however, practically, it should be interpreted as a spurious voxel. (Node-to-node branch)
 4. A node-to-node branch connects two node voxels. (Node-to-end branch)
 5. A node-to-end connects a node voxel with an end voxel. (End-to-end branch)
 6. An end-to-end branch connects two end voxels ("isolated" branch).
  
  Arguments:
  ---------
  image_data: A 3D matrix of type BYTE representing the binary skeleton image to label. The image should be organized having  value 255 on skeleton voxels and zero elsewhere.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
	"""
	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dSkeletonLabeling(image_data,out_image,dimx,dimy,dimz, None)
	py_printErrorMessage(err_code)
	return out_image
   
def py_p3dGVFSkeletonization(image_data, dimx, dimy, dimz = 0, mu= 0.15, eps= 1E-4, hierarc = 1.0, scale= 0.0):
	"""
 
  Computes the skeleton of a 3D binary image using the Brun and Dreossi algorithm [1].
  
  Syntax:
  ------
  Result = py_p3dGVFSkeletonization ( image_data, dimx, dimy [, dimz = value] [, mu=value ] [, eps=value ] [, hierarc=value ] [, scale=value  ])
  
  Return Value:
  ------------
 Returns a volume with the same dimensions of input volume having value 255 on skeleton voxels and 0 elsewhere.
 
  Arguments:
  ---------
  image_data: A 3D matrix of type BYTE.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  mu: A decimal value in the range [0.0 0.5] (default 0.15) that controls the computation of the Gradient Vector Flow (GVF) (see  [1] for details). Changes in the MU value affect the "smoothness" of the output skeleton but poorly affects computational  requirements.
 
  eps: A decimal value in the range [0.0 0.1] (default 1E-4) that controls the iterative process of the computation of the  Gradient Vector Flow (GVF) (see [1] for details). An increase of this value (for instance EPS = 1E-3) greatly improves the  computational performances but may affect the "quality" of the output skeleton.
 
  hierarc: By specifying a decimal value in the range [0.0 0.5] (default 0.0) additional node-to-end branches will be added in  the output skeleton. The addition of more node-to-end branches requires additional computational time.
 
  scale: Since the algorithm internally operates on a continuous model, by specifying a SCALE value (default = 1.0) the discrete  input volume can be downsampled or upsampled as an attempt to, respectively, simplify or enrich the output skeleton (which,  however, will still be one voxel wide and with the same dimensions of the non-resampled input volume). The downsampling  improves performances, while upsampling of the input volume might be useful in a limited resolution regime.
  
  Remarks:
  -------
  Computational time depends on the "complexity" of the structure to skeletonize and not only on volume dimensions.  Skeletonization with this method might require long computational time, especially in the last step of the process.
 
  References
  -------
 [1] F. Brun and D. Dreossi, Efficient curve-skeleton computation for the analysis of biomedical 3D images, Biomedical Sciences  Instrumentation, Vol. 46, pp. 475-80, 2010.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dGVFSkeletonization(image_data, out_image, dimx, dimy, dimz, mu, eps, hierarc, scale, None)
	py_printErrorMessage(err_code)
	return out_image

def py_p3dLKCSkeletonization(image_data, dimx, dimy, dimz = 0):
	"""
 
  Computes the skeleton of a 3D binary image using the Palágyi and Kuba algorithm [1].
  
  Syntax:
  ------
  Result = p3dPKSkeletonization ( image_data, dimx, dimy [, dimz = value] )
  
  Return Value:
  ------------
  Returns a volume with the same dimensions of input volume having value 255 on skeleton voxels and 0 elsewhere.
 
  Arguments:
  ---------
  image_data: A 3D matrix of type BYTE.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  Remarks:
  -------
  The result is a homotopic one voxel thin skeleton without assurance of medialness.
 
  References
  -------
 [1] [1] K. Palágyi and A. Kuba, A Parallel 3D 12-Subiteration Thinning Algorithm, Graphical Models and Image Processing, Vol.  61, pp. 199-221, 1999.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dLKCSkeletonization(image_data,out_image,dimx,dimy,dimz, None)
	py_printErrorMessage(err_code)
	return out_image

def py_p3dSkeletonAnalysis(image_data, skeleton_image, dimx, dimy, dimz = 0, nodes_im = None, pores_im= None, ends_im= None, throats_im= None, merging_factor= 0.85, tortuosity_depth= 3, resolution= 1.0, skel_stats_file = "skeleton_stats.txt"):
	"""
  Performs a series of analysis on the input volume based on its skeleton.
 
 By scanning the skeleton it is possible to extract the number of nodes and branches, length and thickness measures based on the  concept of maximal inscribed sphere [1] as well as connectivity indexes. This approach is typically suitable for the analysis  of porous media having an interconnected porous space. Conceptually, the nodes correspond to pore bodies and the branches of  the pore space skeleton correspond to the channels (or paths) connecting the pores. Unfortunately, since operatively a node is  a skeleton voxel with more than two voxels in its neighborhood (see p3dSkeletonLabeling), this 1-to-1 correspondence does not  exist. In fact, while every channel has its corresponding branch in the skeleton, not every branch corresponds to a channel.  Furthermore, while each pore body is represented by some nodes of the skeleton, several nodes may occur in the same pore body  and these nodes are usually connected through very small (and therefore spurious) branches. A merging criterion is therefore  required in order to merge nodes that occur within the same pore body and in order to avoid to consider the related spurious  branches.
 
 The algorithm used in p3dSkeletonAnalysis merges two or more nodes if there is any overlap among the maximal spheres centered  at the nodes. The set of these overlapped spheres (cluster of spheres) determines a subvolume within the pore body. The largest  maximal sphere is then searched within the previously determined subvolume and it is assumed as the center of a pore. In order  to tune the amount of merging, the user can control the size of the cluster of spheres via the MERGING_FACTOR parameter.
  
  Syntax:
  ------
  Result = p3dSkeletonAnalysis ( image_data, skeleton_image, dimx, dimy [, dimz = value] [, nodes_im=bytearray] [, pores_im=bytearray] [,  throats_im=bytearray] [, pores_im=bytearray][, merging_factor=value] [, tortuosity_depth=value] [, resolution=value,] [skel_stats_file = value] )
  
  Return Value:
  ------------
  Returns a struct of parameters. The name of the fields are:
 CONNECTIVITY_DENSITY [mm-3]: A scalar value representing the number of redundant connections normalized to the total volume V.  It is computed as (1 - ΧV )/V where ΧV = (n - b), being n the number of pores and b the number of node-to-node branches.
 
 COORDINATION_NUMBER [-]: An array of length PORES_COUNT containing the number of branches that spread out from each node.
 
 PORES_COUNT [-]: An integer value representing the number of pores determined after the application of the merging criterion.  Therefore, it does not necessarly correspond to the number of skeleton nodes.
 
 PORES_WIDTH [mm]: An array of length PORES_COUNT containing the pore-size distribution computed as diameter of the maximal  inscribed sphere for each pore. The center of the maximal sphere is affected by the merging criterion.
 
 ENDPOINTS_COUNT [-]: An integer value representing the number of skeleton end points.
 
 ENDPOINTS_WIDTH [mm]: An array of length ENDPOINTS_COUNT containing the width of each end point computed as the diameter of the  maximal sphere centered on the end point.
 
 ENDTOEND_COUNT [-]: An integer value representing the number of end-to-end branches.
 
 ENDTOEND_LENGTH [mm]: An array of length ENDTOEND_COUNT containing the length of each end-to-end branch computed from the  surface to the maximal sphere of an end point to the surface of the maximal sphere of the other end point.
 
 ENDTOEND_MEANWIDTH [mm]: An array of length ENDTOEND_COUNT containing the mean width of each endToEndBranches. The width is  computed averaging the diameter of the maximal spheres of each branch voxel.
 
 ENDTOEND_MINWIDTH [mm]: An array of length ENDTOEND_COUNT containing the minimum width of each end-to-end branch. This value is  the diameter of the smallest maximal spheres among all the maximal spheres centered on each branch voxel.
 
 ENDTOEND_MAXWIDTH [mm]: An array of length ENDTOEND_COUNT containing the maximum width of each end-to-end branch. This value is  the diameter of the largest maximal spheres among all the maximal spheres centered on each branch voxel.
 
 NODETOEND_COUNT [-]: An integer value representing the number of node-to-end branches.
 
 NODETOEND_LENGTH [mm]: An array of length NODETOEND_COUNT containing the length of each node-to-end branch computed from the  surface to the maximal sphere of the node point to the surface of the maximal sphere of the end point.
 
 NODETOEND_MEANWIDTH [mm]: An array of length NODETOEND_COUNT containing the mean width of each node-to-end branch. The width is  computed averaging the diameter of the maximal spheres of each branch voxel.
 
 NODETOEND_MINWIDTH [mm]: An array of length NODETOEND_COUNT containing the minimum width of each node-to-end branch. This value  is the diameter of the smallest maximal spheres among all the maximal spheres centered on each branch voxel.
 
 NODETOEND_MAXWIDTH [mm]: An array of length NODETOEND_COUNT containing the maximum width of each node-to-end branch. This value  is the diameter of the largest maximal spheres among all the maximal spheres centered on each branch voxel.
 
 NODETONODE_COUNT [-]: An integer value representing the number of node-to-node branches.
 
 NODETONODE_LENGTH [mm]: An array of length NODETONODE_COUNT containing the length of each node-to-node branch computed from the  surface of the maximal sphere inscribed within the pore to the surface of the maximal sphere of the other pore.
 
 NODETONODE_MEANWIDTH [mm]: An array of length NODETONODE_COUNT containing the mean width of each node-to-node branch. The width  is computed averaging the diameter of the maximal spheres of each branch voxel.
 
 NODETONODE_MINWIDTH [mm]: An array of length NODETONODE_COUNT containing the minimum width of each node-to-node branch. This  value is the diameter of the smallest maximal spheres among all the maximal spheres centered on each branch voxel. The smallest  thickness along a node-to-node branch is usually defined as throat.
 
 NODETONODE_MAXWIDTH [mm]: An array of length NODETONODE_COUNT containing the maximum width of each node-to-node branch. This  value is the diameter of the largest maximal spheres among all the maximal spheres centered on each branch voxel.
 
  Arguments:
  ---------
  image_data: A 3D matrix of type BYTE.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  SkelImage: A 3D matrix of type BYTE representing the skeleton of OrigImage. OrigImage and SkelImage should have same  dimensions.
 
  nodes_im: Optionally, a 3D matrix of type BYTE with the same dimensions of input images having the filled and sometimes  overlapping balls centered on each skeleton nodes can be returned as output. The size of these balls is tuned with the  MERGING_FACTOR parameter. This volume is the starting point for the determination of the PORES_IM image.
 
  pores_im: Optionally, a 3D matrix of type BYTE with the same dimensions of input images having the maximal balls on each of  what have been considered as pores can be returned as output. The diameter of these balls are used for the pore size  distribution assessment. The center of each of these balls lies within the cluster of overlapped balls of nodes_im.
 
  throats_im: Optionally, a 3D matrix of type BYTE with the same dimensions of input images having the filled maximal balls on  skeleton throats (i.e. the minimum thickness along a node-to-node branch) can be returned as output.
 
  throats_im: Optionally, a 3D matrix of type BYTE with the same dimensions of input images having the filled maximal balls on  skeleton throats (i.e. the minimum thickness along a node-to-node branch) can be returned as output.
 
  merging_factor: A decimal value in the range [0.0,1.0] for reducing the size of the maximal balls to use for merging adjacent  nodes (default 0.85). If the value 1.0 is specified too many skeleton nodes might be merged leading to inaccurate results. On  the other hand, if a too small value is specified, spurious branches might be considered.
 
  resolution: A decimal value representing the resolution of input image. If this value is not specified a voxelsize of 1.0 is  assumed that means output values are expressed in voxel unit.
  
  skel_stats_file: tab_delimeted/categorized file where output statistics are written. If no file address is specifed, the statustics will be saved in the default folder under the name "skeleton_stats.txt"
 
  Remarks:
  -------
  The result is a homotopic one voxel thin skeleton without assurance of medialness.
 
  References
  -------
 [1] K. Palágyi and A. Kuba, A Parallel 3D 12-Subiteration Thinning Algorithm, Graphical Models and Image Processing, Vol.  61,  pp. 199-221, 1999.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	skeleton_stats = malloc_PSkeletonStats()
	err_code=p3dSkeletonAnalysis(image_data,skeleton_image,skeleton_stats,nodes_im,pores_im,ends_im,throats_im,dimx,dimy,dimz,merging_factor,tortuosity_depth, resolution,None)
	py_printErrorMessage(err_code)
	PrintSkelStruct(skeleton_stats, skel_stats_file)
	return skeleton_stats

    
def py_p3dSkeletonPruning(image_data, dimx, dimy, dimz = 0, thresh=3, ultimate = False, iterative = False):
	"""
  Prunes the skeleton by removing node-to-end branches.
 
 Skeleton pruning can be performed by simply fixing a threshold on the length (maximum number of voxels) of a node-to-end  branch, therefore a node-to-end branch having a number of voxels below the specified value is removed from output skeleton.  After first step of pruning, new node-to-end branches could be created and they could be considered spurious as well. If the  ITERATIVE keyword is specified, the pruning algorithm iteratively scans for node-to-end branches after a first removal and it  stops when no more changes are possibile. If the ULTIMATE keyword is specified all the node-to-end branches are removed from  the skeleton independently from their length. If both ITERATIVE and ULTIMATE keywords are specified, only closed paths or end- to-end branches are preserved in the output skeleton. If the input skeleton does not feature closed paths or end-to-end  branches, an application of the ULTIMATE and/or ITERATIVE keywords results in a blank output.
  
  Syntax:
  ------
  Result = py_p3dSkeletonPruning ( image_data, dimx, dimy [, dimz = value] [,thresh=value] [,ultimate = bool] [,iterative = bool])
  
  Return Value:
  ------------
  Returns a volume with the same dimensions of input volume having value 255 on skeleton voxels and 0 elsewhere.
 
  Arguments:
  ---------
  image_data: A 3D matrix of type BYTE.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	if iterative == True and ultimate == False:
		print("p3dIterativeSkeletonPruning")
		err_code = p3dIterativeSkeletonPruning (image_data, out_image,dimx,dimy,dimz, thresh, None )
		py_printErrorMessage(err_code)
		return out_image
	if ultimate == True:
		print("p3dUltimateSkeletonPruning")
		err_code =  p3dUltimateSkeletonPruning (image_data, out_image, dimx,dimy,dimz, iterative, None )
		py_printErrorMessage(err_code)
		return out_image
    # Default: iterative == False and ultimate = False
	print("p3dSimpleSkeletonPruning")
	err_code = p3dSimpleSkeletonPruning(image_data, out_image,dimx,dimy,dimz, thresh, None )
	py_printErrorMessage(err_code)
	return out_image
    

  
  
