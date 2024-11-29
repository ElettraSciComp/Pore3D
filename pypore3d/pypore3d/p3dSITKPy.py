import argparse
import SimpleITK as sitk
import os
import tempfile


import pypore3d.p3d_SITK_common_lib
from pypore3d.p3d_SITK_common_lib import *

import pypore3d.p3d_SITK_read_raw
from pypore3d.p3d_SITK_read_raw import *

import pypore3d.p3dFiltPy
from pypore3d.p3dFiltPy import py_p3dReadRaw8,py_p3dWriteRaw8


########################################## Filters ##########################################
###### median filter
def py_p3d_SITK_Median(img, dimx,dimy,dimz, kWidth = 1):
    img = p3d_to_sitk_file_format(img, dimx,dimy,dimz)
    median_filter = sitk.MedianImageFilter()
    median_filter.SetRadius(kWidth)
    outImg = median_filter.Execute(img)   
    outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
    return outImg


###### Binary Dilation
def py_p3d_Dilate(input_image, dimx, dimy, dimz = 0, kWidth = 3):
	"""
 
  Dilation is one of the two basic operators in the area of mathematical morphology. The basic effect of the operator on a  binary image is to gradually enlarge the boundaries of regions of foreground voxels. Thus areas of foreground voxels grow in  size while holes within those regions become smaller.
  
  Syntax:
  ------
  Result = py_p3d_sitk_Dilate ( input_image, dimx, dimy [, dimz=0] [, kWidth =value ])
  
  Return Value:
  ------------
  Returns a dilated image with the same dimensions and type of input imaged.                  
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  kWidth: An odd integer value in the range [3, 51] (default: 3) representing the diameter of the spherical (or circular in case  of 2D input) structuring element.
 
  Remarks:
  -------
  Often a 3x3x3 structuring element (or a 3x3 for 2D images) is used, although larger kernels (e.g. a 5x5x5 element) can be used  for more severe dilation. Note that a small kernel can be applied more than once in order to produce a similar but not  identical effect as a single pass with a large kernel.
 
	"""
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	dilation_filter = sitk.BinaryDilateImageFilter()
	dilation_filter.SetKernelRadius(kWidth)
	dilation_filter.SetForegroundValue (255)
    
	outImg = dilation_filter.Execute(outImg)    
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### Binary Erosion
def py_p3d_Erode(input_image, dimx, dimy, dimz = 0, kWidth = 3):
	"""
 
  Erosion is one of the two basic operators in the area of mathematical morphology. The basic effect of the operator on a binary  image is to erode away the boundaries of regions of foreground voxels. Thus areas of foreground voxels shrink in size, and  holes within those areas become larger.
  
  Syntax:
  ------
  Result = py_p3d_sitk_Erode ( input_image, dimx, dimy[, dimz = 0 ] [, kWidth =value ])
  
  Return Value:
  ------------
  Returns a dilated image with the same dimensions and type of input imaged.                  
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the binary input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  kWidth: An odd integer value in the range [3, 51] (default: 3) representing the diameter of the spherical (or circular in case  of 2D input) structuring element.
 
  Remarks:
  -------
  Often a 3x3x3 structuring element (or a 3x3 for 2D images) is used, although larger kernels (e.g. a 5x5x5 element) can be used  for more severe erosion. Note that a small kernel can be applied more than once in order to produce a similar but not  identical effect as a single pass with a large kernel.
 
	"""   
	outImg = p3d_to_sitk_file_format(input_image, dimx, dimy, dimz )
	erosion_filter = sitk.BinaryErodeImageFilter()
	erosion_filter.SetKernelRadius(kWidth)
	erosion_filter.SetForegroundValue (255)
    
	outImg = erosion_filter.Execute(outImg)   
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### HMinimaImageFilter
def py_p3d_HMinimaFilter(input_image, dimx, dimy, dimz = 0, height = 3): 
	"""
 
  Perform H-minima transform of input image, i.e. all minima in input image whose depth is less than the specified threshold are  suppressed [1].
 
 H-minima transformation is a mere application of the concept of morphological reconstruction [2]. By performing the  morphological reconstruction using the input image as the mask image and the result of subtraction of the threshold value from the input image as the marker image the H-maxima transform is performed. The H-minima transformed image can be simply obtained by complementing both input and output image of the process. 
  
  Syntax:
  ------
  Result = py_p3d_HMinimaFilter ( input_image, dimx, dimy[, dimz = 0 ] [, threshold =value ])
  
  Return Value:
  ------------
  Returns the H-minima transformed image with the same dimensions and type of input image.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  height: Set the height that a local maximum must be above the local background (local contrast) in order to survive the      processing. Local maxima below this value are replaced with an estimate of the local background. Default = 3
 
  References
  ----------
 
 [1] P. Soille, Morphological Image Analysis: Principles and Applications, Springer-Verlag, 1999, pp. 170-171.
 
 [2] L. Vincent, Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms, IEEE  Transactions on Image Processing, Vol. 2, No. 2, pp. 176-201, 1993.
 
	"""    
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.HMinimaImageFilter()  
	Filter.SetHeight(threshold)
    
	outImg = Filter.Execute(outImg)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### Mutli Thresholding
def py_p3d_MultiThresholding(input_image, dimx, dimy, dimz = 0, regionsNum = 2): 
	"""
 
  Performs multithresholding,
  
  Syntax:
  ------
  Result = py_p3d_MultiThresholding ( input_image, dimx, dimy[, dimz = 0 ] [, regionsNum =value ])
  
  Return Value:
  ------------
  Returns segmented image with the same dimensions and type of input image.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  regionsNum: A number representing segments number.
 
	"""    
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	thresholdFilter = sitk.OtsuMultipleThresholdsImageFilter()
	thresholdFilter.SetNumberOfThresholds(regionsNum)
	outImg = thresholdFilter.Execute(outImg)
    
	rescaler = sitk.RescaleIntensityImageFilter()
	rescaler.SetOutputMinimum(0)
	rescaler.SetOutputMaximum(255)
    
	outImg = rescaler.Execute(outImg)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### GradientMagnitudeImageFilter
def py_p3d_GradientMagnitudeImageFilter(input_image, dimx, dimy, dimz = 0): 
	"""
 
  Performs GradientMagnitudeImageFilter,
  
  Syntax:
  ------
  Result = py_p3d_GradientMagnitudeImageFilter ( input_image, dimx, dimy[, dimz = 0 ] )
  
  Return Value:
  ------------
  Returns filtered image with the same dimensions and type of input image.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
	"""
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.GradientMagnitudeImageFilter()  
	#Filter.SetUseImageSpacing (False)
 
	outImg = Filter.Execute(outImg)    
	#outImg =  apply_rescaler(outImg, dimx,dimy,dimz)   
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
    
	return outImg


###### CurvatureFlowImageFilter
def py_p3d_CurvatureFlowImageFilter(input_image, dimx, dimy, dimz = 0, time_step = 3.0, iterations= 3): 
	"""
 
  The sitk::CurvatureFlowImageFilter performs edge-preserving smoothing in a similar fashion to the classical anisotropic  diffusion.
 
  The filter uses a level set formulation where the iso-intensity contours in a image are viewed as level sets, where pixels of  a  particular intensity form one level set. The level set function is then evolved under the control of a diffusion equation  where the speed is proportional to the curvature of the contour. Areas of high curvature will diffuse faster than areas of low  curvature. Hence, small jagged noise artifacts will disappear quickly, while large scale interfaces will be slow to evolve,  thereby preserving sharp boundaries between objects. However, it should be noted that although the evolution at the boundary is  slow, some diffusion still occur.
  
  Syntax:
  ------
  Result = py_p3d_CurvatureFlowImageFilter ( input_image, dimx, dimy[, dimz = 0 ] [, iterations =value ])
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  iterations: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the domain gaussian.
 
	"""
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.CurvatureFlowImageFilter()  
	Filter.SetNumberOfIterations(iterations)
	Filter.SetTimeStep (time_step)
    
	outImg = Filter.Execute(outImg)
	#outImg =  apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg


###### IsolatedConnectedImageFilter 
def py_p3d_IsolatedConnectedImageFilter (input_image, dimx, dimy, dimz = 0): 
	"""
 
  Pore3D python wrapper for sitk.IsolatedConnectedImageFilter()
  Details: https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1IsolatedConnectedImageFilter.html#details
 
  Syntax:
  ------
  Result = py_p3d_IsolatedConnectedImageFilter ( input_image, dimx, dimy[ , dimz ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
	"""   
	input_image = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.IsolatedConnectedImageFilter()  
    
	#Filter = sitk.SetFindUpperThreshold (True)
	#Filter = sitk.SetIsolatedValueTolerance (double IsolatedValueTolerance)
	#Filter = sitk.SetLower (double Lower) 
	#Filter = sitk.SetReplaceValue (100)
	#Filter = sitk.SetSeed1 (std::vector< unsigned int > Seed1)
	#Filter = sitk.SetSeed2 (std::vector< unsigned int > Seed2)
	#Filter = sitk.SetUpper (double Upper)
    
	outImg = Filter.Execute(input_image)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### NeighborhoodConnectedImageFilter
def py_p3d_NeighborhoodConnectedImageFilter(input_image, dimx, dimy, dimz = 0): 
	"""
 
  Pore3D python wrapper for sitk.NeighborhoodConnectedImageFilter()  
 
  Syntax:
  ------
  Result = py_p3d_NeighborhoodConnectedImageFilter ( input_image, dimx, dimy[ , dimz ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
	"""     
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.NeighborhoodConnectedImageFilter()  
    
	outImg = Filter.Execute(outImg)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### ConnectedThresholdImageFilter
def py_p3d_ConnectedThresholdImageFilter(input_image, dimx, dimy, dimz = 0): 
	"""
 
  Pore3D python wrapper for sitk.ConnectedThresholdImageFilter()  
 
  Syntax:
  ------
  Result = py_p3d_ConnectedThresholdImageFilter ( input_image, dimx, dimy[ , dimz ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
	"""       
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.ConnectedThresholdImageFilter()  
    
	outImg = Filter.Execute(outImg)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

###### sitk_WatershedSegmentation_8
def py_p3d_WatershedSegmentation(input_image, dimx, dimy, dimz = 0, level = 0, markWatershedLine = False, connected = False): 
	"""
 
  Pore3D python wrapper for sitk.MorphologicalWatershedImageFilter()  
 
  Syntax:
  ------
  Result = py_p3d_WatershedSegmentation ( input_image, dimx, dimy[ , dimz ] [, level = value] [, markWatershedLine = bool] [,  connected = bool])
  
  Return Value:
  ------------
  Returns a segmented image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: AA 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  level:
  
  markWatershedLine:
 
  connected:
 
	"""           
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.MorphologicalWatershedImageFilter() 
	Filter.SetFullyConnected(connected)
	Filter.SetLevel(level)
	Filter.SetMarkWatershedLine(markWatershedLine)
    
	outImg = Filter.Execute(outImg)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg


###### k means
def py_p3d_KMeansClustering(input_image, dimx, dimy, dimz = 0, segment_num = 1): 
	"""
 
  Performs K-means clustering multiphase segmentation.
 
 Classification includes a broad range of decision-theoretic approaches to the identification of images (or parts thereof). All  classification algorithms are based on the assumption that the image in question depicts one or more features and that each of  these features belongs to one of several distinct and exclusive classes. The classes may be specified a priori by the user (as  in supervised classification) or automatically clustered (as in unsupervised classification) into sets of prototype classes,  where the user merely specifies the number of desired categories. Within image processing context, classification should be  thought as a multiphase segmentation process: in fact the classes are then used for generating a labeled image in which every  voxel is assigned to one of the classes.
 
 The K-means algorithm converts an input image into vectors of equal size and then determines the k prototype mean vectors by  minimizing of the sum of the squared distances from all points in a class to the class center.
   
  Syntax:
  ------
  Result = py_p3d_KMeansClustering ( input_image, dimx, dimy[ , dimz ])
  
  Return Value:
  ------------
  Returns a segmented image of type BYTE with the same dimensions of input image.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  segment_num: number of segments
 
	"""   
	Img = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	#Img = Read_Raw8(file_name, [dimx,dimy,dimz])
    
	#Filter1 = sitk.RelabelComponentImageFilter()
	#Img = Filter1.Execute(Img)
      
	Filter = sitk.ScalarImageKmeansImageFilter()
	Filter. SetUseNonContiguousLabels (True)
    
	if segment_num <=0:
		print ("Error: segment_num must be greater than 0")
		return
    
	shift_val = 250.0/segment_num
	mean_val = 0.0
	initial_mean_vec =  []
	for i in range(segment_num):
		initial_mean_vec.append(mean_val)
		mean_val = mean_val+shift_val
	for i in range(segment_num):
		print (initial_mean_vec[i])
	#initial_mean_vec =  [0.0,50.0,100.0,150.0,200.0,250.0]
	Filter.SetClassWithInitialMean(initial_mean_vec)
        
	outImg = Filter.Execute(Img)
    
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg


###### CurvatureAnisotropicDiffusionImageFilter (16 bit only)
def py_p3d_CurvatureAnisotropicDiffusionImageFilter(input_image, dimx, dimy, dimz = 0, conductance=7.0, iterations=1): 
	"""
 
 Performs anisotropic diffusion on an image using a modified curvature diffusion equation (MCDE).
 
 Qualitatively, MCDE compares well with other non-linear diffusion techniques. It is less sensitive to contrast than classic  anisotropic diffusion and preserves finer detailed structures in images. Each iteration of the filter requires more  computational time than the classic anisotropic diffusion, however. fewer iterations may be required to reach an acceptable  solution. 
 
  Syntax:
  ------
  Result = py_p3d_CurvatureAnisotropicDiffusionImageFilter ( input_image, dimx, dimy[ , dimz ] [, conductance = value ] [,  iterations = value ])
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  conductance: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the domain gaussian.
 
  iterations: An integer value greater than 0 (default: 3.0) representing the the diffusion time, i.e. the number of times the  algorithm is iterated.
 
	"""       
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	#feature_img = sitk.GradientMagnitude(img)
	Filter = sitk.CurvatureAnisotropicDiffusionImageFilter()  
	Filter.SetNumberOfIterations(iterations)
	Filter.SetConductanceParameter (conductance)

	outImg = Filter.Execute(outImg)    
	outImg =  apply_rescaler(outImg, dimx,dimy,dimz)   
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
    
	return outImg

###### GradientAnisotropicDiffusionImageFilter(16 bit only)
def py_p3d_GradientAnisotropicDiffusionImageFilter(input_image, dimx, dimy, dimz = 0, conductance = 7.0, iterations = 1): 
	"""
 
  Pore3D python wrapper for sitk.GradientAnisotropicDiffusionImageFilter
 
  Syntax:
  ------
  Result = py_p3d_GradientAnisotropicDiffusionImageFilter ( input_image, dimx, dimy[ , dimz ] [, conductance = value ] [, iterations   = value] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  conductance: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the domain gaussian.
 
  iterations: An integer value greater than 0 (default: 3.0) representing the the diffusion time, i.e. the number of times the  algorithm is iterated.
 
	"""     
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.GradientAnisotropicDiffusionImageFilter()  
	Filter.SetNumberOfIterations(iterations)
	Filter.SetConductanceParameter (conductance)
    
	outImg = Filter.Execute(outImg)
	outImg =  apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg


###### MinMaxCurvatureFlowImageFilter (16 bit only)
def py_p3d_MinMaxCurvatureFlowImageFilter(input_image, dimx, dimy, dimz = 0, iterations = 3, width = 3): 
	"""
 
  The sitk::CurvatureFlowImageFilter performs edge-preserving smoothing applying a variant of the curvature flow algorithm where  diffusion is turned on or off depending on the scale of the noise.
 
 The minimum-maximum variant of the curvature flow filter results in sharper edges than the application with the simple  curvature flow with similar parametrization.
 
  Syntax:
  ------
  Result = py_p3d_MinMaxCurvatureFlowImageFilter ( input_image, dimx, dimy[ , dimz ] [, iterations = value ] [, width = value] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input imaged.                 
  
  Arguments:
  ---------
  input_image: A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  iterations: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the domain gaussian.
 
  width: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the domain gaussian.
 
	"""      
	outImg = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
	Filter = sitk.MinMaxCurvatureFlowImageFilter()  
	Filter.SetNumberOfIterations(iterations)
	Filter.SetStencilRadius(width)
    
	outImg = Filter.Execute(outImg)
	outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	return outImg

##### MinMaxCurvatureFlowImageFilter (16 bit only)
#def py_p3d_MinMaxCurvatureFlowImageFilter_16(input_image, dimx,dimy,dimz = 0, ITERATIONS=3, WIDTH = 3): 
#    
#    input_image = p3d_to_sitk_file_format(input_image, dimx,dimy,dimz)
#    Filter = sitk.MinMaxCurvatureFlowImageFilter()  
#    Filter.SetNumberOfIterations(ITERATIONS)
#    Filter.SetStencilRadius(WIDTH)
#    
#    outImg = Filter.Execute(input_image)
#    outImg = apply_rescaler(outImg, dimx,dimy,dimz)
#    outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
#    return outImg



#FUNCTION        P3DCONFIDENCEREGIONGROWING				1       1   KEYWORD#S