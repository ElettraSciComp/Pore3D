import argparse
import SimpleITK as sitk
import os
import tempfile


import pypore3d.p3d_SITK_common_lib
from pypore3d.p3d_SITK_common_lib import *
from pypore3d.p3d_SITK_common_lib_16 import *


import pypore3d.p3d_SITK_read_raw
from pypore3d.p3d_SITK_read_raw import *

import pypore3d.p3dFiltPy
from pypore3d.p3dFiltPy import py_p3dReadRaw8,py_p3dWriteRaw8

import pypore3d.p3dFiltPy_16
from pypore3d.p3dFiltPy_16 import py_p3dReadRaw16,py_p3dWriteRaw16


########################################## Filters ##########################################


def py_p3d_SITK_Median_16(img, dimx,dimy,dimz, kWidth = 1):
	"""
	"""
	median_filter = sitk.MedianImageFilter()
	median_filter.SetRadius(kWidth)
	outImg = median_filter.Execute(img)   
	return outImg




###### sitk_WatershedSegmentation
###### sitk_WatershedSegmentation_16
def py_p3d_WatershedSegmentation(dist_file_path, watershed_image_path, dimx, dimy, dimz = 0, level = 0, markWatershedLine = True, connected = False): 
	"""
	Syntax:
	------
	Result = py_p3d_WatershedSegmentation ( dist_file_path, watershed_image_path, dimx, dimy [, dimz = value] [, level = value] [, markWatershedLine = value] [, connected = value])
  
	Return Value:
	------------
	Returns the H-py_p3d_WatershedSegmentation transformed image with the same dimensions and type of input image.                 
  
	Arguments:
	---------
	dist_file_path: File path storing A 2D or 3D matrix representing the input image. Usually it is a distance field image
 
	dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
	markWatershedLine: Default = True
    
	connected: Default = True
  
	watershed_image_path: the path where the output image is written
	""" 
	dist_img = Read_Raw16(dist_file_path, [dimx,dimy,dimz])

	Filter = sitk.MorphologicalWatershedImageFilter() 
	Filter.SetFullyConnected(connected)
	Filter.SetLevel(level)
	Filter.SetMarkWatershedLine(markWatershedLine)
	WS_Img = Filter.Execute(dist_img)
    
	out_file = watershed_image_path + '.mhd'
	sitk.WriteImage(WS_Img,out_file)
	return WS_Img


###### HMinimaImageFilter
def py_p3d_HMinimaFilter(input_file_path, output_image_path, dimx, dimy, dimz = 0, threshold = 3): 
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
  input_file_path: File path storing A 2D or 3D matrix of type BYTE representing the input image.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  threshold: Set the height that a local maximum must be above the local background (local contrast) in order to survive the      processing. Local maxima below this value are replaced with an estimate of the local background. Default = 3
  
  output_image_path: the path where the output image is written
 
  References
  ----------
 
 [1] P. Soille, Morphological Image Analysis: Principles and Applications, Springer-Verlag, 1999, pp. 170-171.
 
 [2] L. Vincent, Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms, IEEE  Transactions on Image Processing, Vol. 2, No. 2, pp. 176-201, 1993.
 
	"""    
	input_img = Read_Raw16(input_file_path, [dimx,dimy,dimz])
	Filter = sitk.HMinimaImageFilter()  
	Filter.SetHeight(threshold)
    
	outImg = Filter.Execute(input_img)
	#outImg = apply_rescaler(outImg, dimx,dimy,dimz)
	#outImg = sitk_to_p3d_file_format (outImg, dimx,dimy,dimz)
	out_file = output_image_path + '.mhd'
	sitk.WriteImage(outImg,out_file)
	return outImg

