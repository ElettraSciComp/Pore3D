
import pypore3d.p3d_common_lib
from pypore3d.p3d_common_lib import *

import pypore3d.p3dBlob
from pypore3d.p3dBlob import *


#################
  
#Basic analysis
def py_p3dBasicAnalysis(image_data, dimx, dimy, dimz = 0, resolution = 1.0):
	"""
 
  Performs a series of direct 3D basic analysis for randomly distributed porous media. 
  
  Syntax:
  ------
  Result = py_p3dBasicAnalysis ( image_data, dimx, dimy, [dimz = value] [, resolution = value] )
  
  Return Value:
  ------------
  Performs a series of direct 3D basic analysis for randomly distributed porous media. Computed parameters are [1]:
 
  VV [-]: Density (VV). A measure of density based on the number of object voxels with respect to the total number of volume  voxels.
 
  SV [mm-1]: Specific surface area (SV). A measure of the surface of the object with respect to the total volume. Tipically it     is related to the mechanical properties of the object.
 
  MV [mm-2]: Integral of mean curvature (MV). A positive value implies the dominance of convex structures, while MV < 0 occurs  in the case of predominance of concave structures.
 
  CV [mm-3]: Euler characteristic (ΧV). This is an index of connectivity of the object network.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  resolution: A decimal value representing the resolution of input image. If this value is not specified a voxelsize of 1.0 is  assumed that means output values are expressed in voxel unit.
 
  References:
  ---------
 
  [1] J. Ohser and F. Mücklich, Statistical Analysis of Microstructures in Materials Science. Wiley & Sons, 2000.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_BasicStats = BasicStats()
	err_code = p3dBasicAnalysis(image_data,out_BasicStats,dimx,dimy,dimz, resolution, None)
	py_printErrorMessage(err_code)
	return out_BasicStats

#Anisotrtopy analysis
def py_p3dAnisotropyAnalysis(image_data, dimx, dimy, dimz = 0, resolution = 1.0):
	"""
 
  Returns a struct of parameters computed from input binary image. Computed parameters are [1]:
  
  Syntax:
  ------
  Result = py_p3dAnisotropyAnalysis ( image_data, dimx, dimy, [dimz = value] [, resolution = value] [, details = boolean] )
  
  Return Value:
  ------------
  Performs a series of direct 3D basic analysis for randomly distributed porous media. Computed parameters are [1]:
 
  I [-]: sotropy index. It measures the similarity of a fabric to a uniform distribution and varies between 0 (all observation  confined to a single plane or axis) and 1 (perfect isotropy)
 
  E [-]: Elongation index. It measures the preferred orientation of a fabric in the u1/u2 plane and varies between 0 (no  preferred orientation) and 1 (a perfect preferred orientation with all observations parallel).
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  resolution: A decimal value representing the resolution of input image. If this value is not specified a voxelsize of 1.0 is  assumed that means output values are expressed in voxel unit.
 
 
  References:
  ----------
 
  [1] S.C. Cowin and A.J. Laborde, The relationship between the elasticity tensor and the fabric tensor. Mechanics of Materials,   Vol. 4, No. 22, pp. 137-147, 1985.
 
  [2] S.C. Cowin, Wolff's law of trabecular architecture at remodeling equilibrium. Journal of Biomechanical Engineering, Vol.  108, No. 1, pp. 83-88, 1986.
 
  [3] W.J. Whitehouse, The quantitative morphology of anisotropic trabecular bone. Journal of Microscoscopy, Vol. 101, pp. 153-168, 1974.
 
  [4] T.P. Harrigan and R.W. Mann, Characterization of microstructural anisotropy in orthotropic materials using a second rank  tensor. Journal of Material Science, Vol. 19, No. 3, pp. 761-767, 1984.
 
  [5] D. I. Benn, Fabric shape and the interpolation of sedimentary fabric data. Journal of Sedimentary Research, Vol. 64, No.    2, pp. 910-915, 1994.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_AnisotropyStats = AnisotropyStats()
	im_dummy_mask = malloc_uchar(dimx*dimy*dimz)
	im_dummy_mask = image_data
	details = False
	err_code = p3dAnisotropyAnalysis(image_data, None, out_AnisotropyStats, dimx, dimy, dimz, resolution, details , None)
	py_printErrorMessage(err_code)
	return out_AnisotropyStats

#Blob labeling
def py_p3dBlobLabeling(image_data, dimx, dimy, dimz = 0, conn3D = CONN6 ):
	"""
 
  Performs connected components (blobs) labeling using the algorithm proposed in [1].
 
  Connected components labeling scans an image and groups its voxels into components based on voxel connectivity. Once all
 groups have been determined, each voxels is labeled with a graylevel according to the component it was assigned to.
  
  Syntax:
  ------
  Result = py_p3dBlobLabeling ( image_data, dimx, dimy, dimz [, conn3D = value] )
  
  Return Value:
  ------------
  Returns an image of type UINT with the same dimensions of input image in which each connected component is assigned a
 different gray-level.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  conn3D: Specify the desired connectivity: 6 (default), 18 or 26 for 3D images.
 
  Remarks
  -------
  No casting to 8-bit format is performed if less than 256 connected components are determined in an image. Returned output
 remains in 16-bit format in any case. Obviously, this implies that no more than 65535 connected components may be recognized in
 input image.
 
  References:
  ----------
 
  [1] Q. Hu et al., Fast connected-component labeling in three-dimensional binary images based on iterative recursion. Computer   Vision and Image Understanding, Vol. 99, pp. 414-434, 2005.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_ushort(dimx*dimy*dimz)
	skip_borders = 1
	rand = 0
	err_code = p3dBlobLabeling_ushort(image_data, out_image, dimx,dimy,dimz,conn3D, rand , skip_borders, None)
	py_printErrorMessage(err_code)
	return out_image

#Blob Get max
def py_p3dGetMaxVolumeBlob3D(image_data, dimx, dimy, dimz = 0, conn3D = CONN6):
	"""

 Returns an image containing only the connected component (blob) of input binary image having the greatest volume.

 The algorithm performs a connected components labeling with contextually computation of the volume (or area for 2D images) of   each component. An image containing the connected component having only the greatest volume is returned as output.

 Syntax:
 ------
 Result = py_p3dGetMaxVolumeBlob3D ( image_data, dimx, dimy, [dimz = value] [, conn3D = value] )
 
 Return Value:
 ------------
 Returns a filtered image with the same dimensions and type of input imaged.
 
 Arguments:
 ---------

 image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.

 dimx,dimy,dimz: three variables representing the dimensions of image to read. 

 conn3D: The desired connectivity, i.e. 6, 18 or 26 for 3D images. The default connectivity is 6 for three dimensional images.

	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dGetMaxVolumeBlob3D(image_data, out_image, dimx,dimy,dimz,conn3D, None)
	py_printErrorMessage(err_code)
	return out_image

#Blob Get min
def py_p3dGetMinVolumeBlob3D(image_data, dimx, dimy, dimz = 0, conn3D = CONN6):
	"""
     
 Returns an image containing only the connected component (blob) of input binary image having the minimum volume.
 
 The algorithm performs a connected components labeling with contextually computation of the volume (or area for 2D images) of  each component. An image containing the connected component having only the greatest volume is returned as output.
 
  Syntax:
  ------
  Result = py_p3dGetMinVolumeBlob3D ( image_data, dimx, dimy, [dimz = value] [, conn3D = value] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  conn3D: The desired connectivity, i.e. 6, 18 or 26 for 3D images. The default connectivity is 6 for three dimensional images.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dGetMinVolumeBlob3D(image_data, out_image, dimx,dimy,dimz,conn3D, None)
	py_printErrorMessage(err_code)
	return out_image

#ChamferDT
def py_p3dChamferDT(image_data, dimx, dimy, dimz = 0, w1 = 3, w2 = 4, w3 = 5):
	"""
 In a distance transformed image, each object voxel has a value measuring the distance to the nearest background voxel.
 
 The distance transform (DT) is sometimes called "burn number distribution". Imagine that foreground regions in the input binary  image are made of some uniform slow burning inflammable material. Then consider simultaneously starting a fire at all points on  the boundary of a foreground region and letting the fire burn its way into the interior. Labeling each point in the interior  with the amount of time that the fire took to first reach that point, means effectively computing the distance transform of  that region.
 
 There are several different sorts of distance transform, depending upon which distance metric is being used to determine the  distance between voxels. The different distance measures are achieved by using different sets of weights in a sequential  scanning algorithm, as described in [1]. Not all combinations of local distances <w1,w2,w3> result in useful distance  transforms. The most interesting weighted (or chamfer) distance transforms implemented in Pore3D are:
  Cityblock (or "Manhattan") :<1,-,->
  Chessboard: <1,1,1>
  Quasi-Euclidean: <3,4,5>
 
  Syntax:
  ------
  Result = py_p3dChamferDT ( image_data, dimx, dimy, [dimz = value] [, w1 = value] [, w2 = value] [, w3 = value])
  
  Return Value:
  ------------
  Returns the distance transform image with the same dimensions of input image and type BYTE.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  w1, w2, w3: three integer values (default: [3,4,5]) representing the weights of the chamfer distance transform.
 
  References:
  ----------
  [1] G. Borgefors, On Digital Distance Transforms in Three Dimensions, Computer Vision and Image Understanding, Vol. 64, No. 3,   pp. 368-376, 1996.
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_ushort(dimx*dimy*dimz)
	err_code = p3dChamferDT(image_data, out_image, dimx,dimy,dimz,w1, w2, w3, None)
	py_printErrorMessage(err_code)
	return out_image

#SquaredEuclideanDT
def py_p3dSquaredEuclideanDT(image_data, dimx, dimy, dimz = 0):
	"""
  In a distance transformed image, each object voxel has a value measuring the distance to the nearest background voxel.
 
 The distance transform (DT) is sometimes called "burn number distribution". Imagine that foreground regions in the input binary  image are made of some uniform slow burning inflammable material. Then consider simultaneously starting a fire at all points on  the boundary of a foreground region and letting the fire burn its way into the interior. Labeling each point in the interior  with the amount of time that the fire took to first reach that point, means effectively computing the distance transform of  that region.
 
 There are several different sorts of distance transform, depending upon which distance metric is being used to determine the  distance between voxels. The euclidean distance is the one adopted in p3dSquaredEuclideanDT and in order to avoid decimal  number, squared values are returned.
 
  Syntax:
  ------
  Result = py_p3dSquaredEuclideanDT ( image_data, dimx, dimy, [dimz = value])
  
  Return Value:
  ------------
  Returns the distance transform image with the same dimensions of input image and type UINT.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  References:
  ----------
 
  [1] T. Hirata. A unified linear-time algorithm for computing distance maps. Information Processing Letters, 58(3):129-133, May   1996.
 
  [2] A. Meijster, J.B.T.M. Roerdink and W. H. Hesselink. A general algorithm for computing distance transforms in linear time.   Mathematical Morphology and its Applications to Image and Signal Processing, pp. 331-340. Kluwer, 2000.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_ushort(dimx*dimy*dimz)
	err_code = p3dSquaredEuclideanDT(image_data, out_image, dimx, dimy, dimz, None)
	py_printErrorMessage(err_code)
	return out_image

# Morphometric Analysis
def py_p3dMorphometricAnalysis(image_data, dimx, dimy, dimz, resolution = 1.0):
	"""
  Performs a series of direct 3D analysis suitable for trabecular-like porous media.
 
  Syntax:
  ------
  Result = py_p3dMorphometricAnalysis ( image_data, dimx, dimy [, dimz = value] [,resolution = 1.0])
  
  Return Value:
  ------------
  Returns a struct of parameters computed from input binary image. The name of the fields of the struct are the ones according  to "Bone ASBMR" [1] as follows:
 
  BVTV [-]: Bone Volume / Total Volume (BV/TV). The ratio of object voxels and the total number of voxels in the considered  Volume of Interest (VOI). If no irregular VOI is defined, i.e. no mask parameter is applied, the output is identical to the  parameter Density of the command p3dBasicAnalysis.
 
  TBTH [mm]: Trabecular thickness (Tb.Th). A measure of thickness of the solid phase objects obtained with a variation of the  directed secant method exposed in [2], assuming the parallel plate model (see [3]). A similar but model-independent parameter  can be extracted performing skeleton analysis.
 
  TBSP [mm]: Trabecular separation (Tb.Sp). A measure of separation of the solid phase objects (i.e. a measure of "thickness" of  the void phase), obtained with a variation of the directed secant method exposed in [2], assuming the parallel plate model (see  [3]). A similar but model-independent parameter can be extracted performing skeleton analysis.
 
  TBN [mm-1]: Trabecular number (Tb.N). A measure related to the number of traversals across a solid structure, tipically  interpreted as a measure of architectural complexity.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  resolution: A decimal value representing the resolution of input image. If this value is not specified a voxelsize of 1.0 is  assumed that means output values are expressed in voxel unit.
 
  References:
  ----------
 
  [1] A.M. Parfitt et al., Bone histomorphometry: Standardization of nomenclature, symbols and units, Journal of Bone and  Mineral Research, Vol. 2, pp. 595-610, 1987.
 
  [2] C.A. Simmons and J.A. Hipp, Method-based differences in the automated analysis of the three-dimensional morphology of  trabecular bone, Journal of Bone and Mineral Research, Vol. 12, No. 6, pp. 942-947, 1997
 
  [3] A.P. Accardo et al., Medical imaging analysis of the three dimensional (3D) architecture of trabecular bone: techniques  and their applications. Medical imaging system technology: analysis and computational methods, C.T. Leondes, World Scientific,  2005.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_MorphometricStat = MorphometricStats()
	#mask_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dMorphometricAnalysis(image_data, None, out_MorphometricStat, dimx,dimy,dimz, resolution, None)
	py_printErrorMessage(err_code)
	return out_MorphometricStat

# Texture Analysis
def py_p3dTextureAnalysis(image_data, dimx, dimy, dimz = 0):
	"""
  Computes a series of direct 3D textural measures.
 
  Syntax:
  ------
  Result = py_p3dTextureAnalysis ( image_data, dimx, dimy [, dimz = value])
  
  Return Value:
  ------------
  Returns a struct of parameters computed from input binary image. The name of the fields of the struct is:
 
  FB [-]: Fractal dimension. A textural measure based on the fractal theory [1].
 
  References:
  ----------
 
  [1] B. Mandelbrot The fractal geometry of nature. W.H. Freeman and Company, 1982.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_TextureStat = TextureStats()
	err_code = p3dTextureAnalysis(image_data, out_TextureStat, dimx,dimy,dimz, None)
	py_printErrorMessage(err_code)
	return out_TextureStat

# MinVolumeFilter
def py_p3dMinVolumeFilter3D(image_data, dimx, dimy, dimz = 0, min_vol = 5, conn = 6):
	"""
  Removes from a binary image all the connected components (blob) that have fewer than the specified number of voxels (or  pixels).
 
 The algorithm starts performing a connected components labeling (see p3dBlobLabeling) with contextually computation of the  volume of each component. Connected components having volume lower than the specified threshold are removed from the image.
 
  Syntax:
  ------
  Result = py_p3dMinVolumeFilter3D ( image_data, dimx, dimy [, dimz = value] [, resolution = value])
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.
 
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE  representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  min_vol: All connected components having volume (i.e. number of voxels) below this value will be removed (default = 5).
 
  conn: Specify the desired connectivity: 6 (default), 18 or 26.
 
  Remarks:
  ---------
 
  Input binary image is assumed to have value 255 for object voxels (the blobs to characterize) and 0 for background.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_image = malloc_uchar(dimx*dimy*dimz)
	if conn == 6 :
		conn = CONN6;
	if conn == 18 :
		conn = CONN18;
	if conn == 26 :
		conn = CONN26;
	err_code = p3dMinVolumeFilter3D(image_data, out_image, dimx,dimy,dimz, min_vol, conn, None)
	py_printErrorMessage(err_code)
	return out_image


# Blob Analysis
def py_p3dBlobAnalysis(image_data, dimx, dimy, dimz = 0, blob_im = None, star_im = None, resolution = 1.0, conn = 6, blob_analysis_file = "blob_analysis.txt"):
	"""
  Performs a series of direct 3D analysis suitable for porous media having isolated pores. If the pore space is formed by an  isolated set of “blobs” (connected components), a series of descriptors for size and shape of each “blob” can be computed. The  analysis is based on the concept of connected com- ponents and their labeling (see py_p3dBlobLabeling_8).
  
  Syntax:
  ------
  Result = py_p3dBlobAnalysis ( image_data, dimx, dimy [, dimz = value] [, blob_im = bytearray] [, star_im = bytearray] [, resolution =  value] [, conn = value,] [blob_analysis_file = value])
  
  Return Value:
  ------------
  Returns a struct of parameters computed from input binary image. The fields are:
 COUNT [-]: The number of identified blobs.
 
 VOLUME [mm3]: An array of length COUNT with the volume of each identified blob computed as the number of voxels rescaled  according to the specified voxel size.
 
 MAX_SPHERE [mm]: An array of length COUNT with the diameter of the maximum inscribed sphere of each identified blob. It is  computed as two times the maximum value of the Euclidean distance transform within the blob.
 
 EQ_SPHERE [mm]: An array of length COUNT with the diameter of the equivalent sphere, i.e. the diameter of a sphere with the  same volume as the blob. It is computed exploiting the inverse formula of the volume of a sphere.
 
 MIN_AXIS [mm]: An array of length COUNT with the minor axis length, i.e. the length of the shortest segment among all the  segments fully included into the blob and passing through its center of mass. The so-called “star” of segments from which  selecting the shortest is generated using random orientations. The "star" image can be optionally returned as output in order  to determine if more random segments have to be computed.
 
 MAX_AXIS [mm]: An array of length COUNT with the major axis length, i.e. the length of the longest segment among all the  segments fully included into the blob and passing through its center of mass.
 
 SPHERICITY [-]: An array of length COUNT with the ratio of MAX_SPHERE and EQ_SPHERE for each blob.
 
 ASPECT_RATIO [-]: An array of length COUNT with the ratio of MIN_AXIS and MAX_AXIS for each blob.
 
 EXTENT [-]: An array of length COUNT with the ratio between the volume of the blob and the volume of the minimum bounding box,  i.e. the smallest parallelepiped oriented according to image axis containing the blob.
  
  Arguments:
  ---------
 
  image_data: A 2D or 3D matrix of type BYTE representing the input image to write to disk.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  blob_im: A 3D matrix of type BYTE with the same dimensions of input image with the maximal ball on each identified blob can be  returned as output. The diameter of these balls are used for the assessment of Max_Sphere distribution.
 
  star_im: A 3D matrix of type BYTE with the same dimensions of input image with the minor and major axis for each blob can be  returned as output. This image is labeled with gray level = 1 for the center of mass (center of the "star"), gray level = 2 for  the minor axis and gray level = 3 for the major axis.
 
  resolution: A decimal value representing the voxel size of input image. If this value is not specified a voxelsize of 1.0 is  assumed that means output values are expressed in voxel unit.
 
  conn: Specify the desired connectivity: 6 (default), 18 or 26.
  
  blob_analysis_file: tab-delimeted file where the results of the blob statiscs are saved. If no name is specified, the file will be saved under the name "blob_analysis.txt" in the default folder location. 
 
  Remarks:
  ---------
 
  Input binary image is assumed to have value 255 for object voxels (the blobs to characterize) and 0 for background.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	out_BlobStat = BlobStats()
	max_rot = 1024  
	skip_borders = 0
    
	if conn == 6 :
		conn = CONN6;
	if conn == 18 :
		conn = CONN18;
	if conn == 26 :
		conn = CONN26;
        
	if skip_borders == 0 :
		borders = P3D_FALSE;
	else :
		borders = P3D_TRUE;
	err_code = p3dBlobAnalysis(image_data, out_BlobStat, blob_im, star_im, dimx,dimy,dimz,resolution,conn,max_rot,borders, None)
	py_printErrorMessage(err_code)
	PrintBlobStruct(out_BlobStat, blob_analysis_file)
	return out_BlobStat

    
    
