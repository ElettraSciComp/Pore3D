import os.path

import pypore3d.p3d_common_lib
from pypore3d.p3d_common_lib import *

import pypore3d.p3dFilt
from pypore3d.p3dFilt import *


#read
def py_p3dReadRaw8(filename, dimx, dimy , dimz = 0):
	"""
 
  Read a RAW image in 8-bit format from the specified file.
  
  Syntax:
  ------
  Result = py_p3dReadRaw8 ( filename, dimx, dimy [, dimz = value])
  
  Return Value:
  ------------
  Returns a matrix of type BYTE with the dimensions specified in input representing the image read from disk.                  
  
  Arguments:
  ---------
  filename: A string with the filename with full path included.
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  Remarks:
  -------
  Trying to read a RAW image stored in SIGNED format is in principle an error. 
  Also, using p3dReadRaw8 for a 16-bit RAW image results in a completely wrong output but,again, no error could 
 be detected.
  Attention must be paid when using this function. Pore3D internally deals only with UNSIGNED format.
 
	"""
    
	if os.path.isfile(filename) == False:
		py_printErrorMessage(-5)
		return 
    
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return   
    
	image_data = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dReadRaw8(filename, image_data,dimx,dimy,dimz, None, None);
	py_printErrorMessage(err_code)
	return image_data

#write

def py_p3dWriteRaw8(image_data, filename, dimx, dimy, dimz = 0):
	"""
 
  Write a RAW image to the specified path using the specified file name (8-bit).
  
  Syntax:
  ------
  Result = py_p3dWriteRaw8 ( image_data, filename, dimx, dimy [, dimz = value])
  
  Return Value:
  ------------
  No return value.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to write to disk.
  filename: A string with the filename with full path included.
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return
    
	err_code = p3dWriteRaw8(image_data, filename,dimx,dimy,dimz, None, None);
	py_printErrorMessage(err_code)
	return image_data

#Gaussian Filter
def py_p3dGaussianFilter8(image_data, dimx, dimy, dimz = 0, width =3 ,sigma = 1.0):
	"""
 
  The gaussian smoothing operator is a 2-D convolution operator that is used to "blur" images and remove detail and noise. In  this sense it is similar to the mean filter, but it uses a different kernel that represents the shape of a gaussian ("bell- shaped") hump.
 
 The effect of gaussian smoothing is to blur an image, in a similar fashion to the mean filter. The degree of smoothing is  determined by the standard deviation of the gaussian. The gaussian outputs a "weighted average" of each pixel's neighborhood,  with the average weighted more towards the value of the central pixels. This is in contrast to the mean filter's uniformly  weighted average. Because of this, a gaussian provides gentler smoothing and preserves edges better than a similarly sized mean  filter.
  
  Syntax:
  ------
  Result = py_p3dGaussianFilter8 ( image_data, dimx, dimy [, dimz = value] [, width=value] [, sigma=value] )
  
  Return Value:
  ------------
  Returns a smoothed image with the same dimensions and type of input images.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  width : An odd integer value in the range [3, 51] (default: 3) representing the side of the square (or cubic in case of 3D  input) kernel.
 
  sigma=value : A decimal value greater than 0 (default: 1.0) representing the standard deviation of gaussian kernel.
 
  Remarks:
  -------
  Often a 3x3 square kernel (or a 3x3x3 cubic kernel for 3D images) is used, although larger kernels (e.g. 5x5 squares or 5x5x5  cubes) can be used for more severe smoothing. Note that a small kernel can be applied more than once in order to produce a  similar but not identical effect as a single pass with a large kernel.
 
"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return 

	if width < 3 or width > 51 or (width % 2) == 0:
		py_printErrorMessage(-4)
		return   
    
	if sigma < 0 :
		py_printErrorMessage(-4)
		return 
    
	out_image = malloc_uchar(dimx*dimy*dimz)
	err_code = p3dGaussianFilter3D_8(image_data,out_image,dimx,dimy,dimz, 3,1.0, None, None)
	py_printErrorMessage(err_code)
	return out_image			

# Mean Filter
def py_p3dMeanFilter8(image_data, dimx, dimy, dimz = 0, width =3):
	"""
 
  Mean filtering is a simple way to smooth an image, i.e. reducing the amount of intensity variation between one pixel and its  neighbors. It is often used to reduce noise in images.
 
 The idea of mean filtering is simply to replace each pixel value in an image with the average value of its neighbors, including  itself. This has the effect of eliminating pixel values which are unrepresentative of their surroundings. Mean filtering is  usually thought of as a convolution filter. Like other convolutions it is based around a kernel, which represents the size of  the neighborhood to be taked into account when calculating the mean.
  
  Syntax:
  ------
  Result = p3dMeanFilter ( image_data, dimx, dimy [, dimz = value] [, width=value] )
  
  Return Value:
  ------------
  Returns a smoothed image with the same dimensions and type of input images.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  width : An odd integer value in the range [3, 51] (default: 3) representing the side of the square (or cubic in case of 3D  input) kernel.
  
  Remarks:
  -------
  Often a 3x3 kernel (or a 3x3x3 kernel for 3D images) is used, although larger kernels (e.g. 5x5 squares or 5x5x5 cubes) can be  used for more severe smoothing. Note that a small kernel can be applied more than once in order to produce a similar but not  identical effect as a single pass with a large kernel.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	if width < 3 or width > 51 or (width % 2) == 0:
		py_printErrorMessage(-4)
		return   
    
	out_image =  malloc_uchar(dimx*dimy*dimz)
	err_code = p3dMeanFilter3D_8(image_data,out_image, dimx,dimy,dimz, width,None, None)				
	py_printErrorMessage(err_code)
	return out_image
	
# Median Filter
def py_p3dMedianFilter8(image_data, dimx, dimy, dimz = 0, width=3):
	"""
 
  The median filter is normally used to reduce noise in an image, somewhat like the mean filter. However, it often does a better  job than the mean filter of preserving useful detail in the image.
 
 Like the mean filter, the median filter considers each pixel in the image in turn and looks at its nearby neighbors to decide  whether or not it is representative of its surroundings. Instead of simply replacing the pixel value with the mean of  neighboring pixel values, it replaces it with the median of those values. In principle, the median is calculated by first  sorting all the pixel values from the surrounding neighborhood into numerical order and then replacing the pixel being  considered with the middle pixel value.
 
 The median filter has two main advantages over the mean filter: 1) the median is a more robust average than the mean and so a  single very unrepresentative pixel in a neighborhood will not affect the median value significantly; 2) since the median value  must actually be the value of one of the pixels in the neighborhood, the median filter does not create new unrealistic pixel  values when the filter straddles an edge. For this reason the median filter is much better at preserving sharp edges than the  mean filter. In general, the median filter allows a great deal of high spatial frequency detail to pass while remaining very  effective at removing noise on images where less than half of the pixels in a smoothing neighborhood have been effected. As a  consequence of this, median filtering can be less effective at removing noise from images corrupted with Gaussian noise.
  
  Syntax:
  ------
  Result = p3dMedianFilter ( image_data, dimx, dimy [, dimz = value] [, width=value ] )
  
  Return Value:
  ------------
  Returns a smoothed image with the same dimensions and type of input images.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx, dimy, dimz: three variables representing the dimensions of image to read. 
 
  width: An odd integer value in the range [3, 51] (default: 3) representing the side of the square (or cubic in case of 3D  input) kernel.
  
  Remarks:
  -------
  Trying to read a RAW image stored in SIGNED format is in principle an error. 
  Also, using p3dReadRaw8 for a 16-bit RAW image results in a completely wrong output but, again, no error could be detected.
  Attention #must be paid when using this function. Pore3D internally deals only with UNSIGNED format.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	if width < 3 or width > 51 or (width % 2) == 0:
		py_printErrorMessage(-4)
		return   
    
	out_image =  malloc_uchar(dimx*dimy*dimz)
	err_code = p3dMedianFilter3D_8(image_data,out_image, dimx,dimy,dimz, width,None, None)				
	py_printErrorMessage(err_code)
	return out_image

# AnisotropicDiffusionFilter
def py_p3dAnisotropicDiffusionFilter8(image_data, dimx, dimy, dimz = 0, m = 1, lambdaP = 0.01, sigma = 0.01, iterP = 10):
	"""
 
  Performs edge preserving smoothing using edge-enhancing anisotropic diffusion [1].
 
 Anisotropic diffusion is an edge preserving smoothing filter similar to bilateral filter. When used with very small diffusion  parameters and long diffusion time, it is less sensitive to contrast but better preserves finer detailed structures in images.
  
  Syntax:
  ------
  Result = py_p3dAnisotropicDiffusionFilter8 ( image_data, dimx, dimy [, dimz = value] [, m = value ] [, lambdaP = value ] [, sigma=value ]  [, iterP=value ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input images.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  m: A integer value greater than 0 (default: 1) representing the standard deviation of the domain gaussian
 
  lambdaP: A decimal value greater than 0 (default: 0.01) representing the in the diffusion equation. Higher values imply faster  diffusion with more smoothing of tiny details.

  sigma: A decimal value greater than 0 (default: 0.01) representing the in the diffusion equation. Higher values imply faster  diffusion with more smoothing of tiny details.
  
  iterP: An integer value greater than 0 (default: 10) representing the diffusion time, i.e. the number of times the algorithm  is iterated.
  
  Remarks:
  -------
  Intensive memory occupation.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	if m < 0:
		py_printErrorMessage(-4)
		return   

	if lambdaP < 0:
		py_printErrorMessage(-4)
		return   
    
	if sigma < 0:
		py_printErrorMessage(-4)
		return 

	if iterP < 0:
		py_printErrorMessage(-4)
		return 
    
	out_image =  malloc_uchar(dimx*dimy*dimz)
	err_code = p3dAnisotropicDiffusionFilter3D_8(image_data,out_image, dimx,dimy,dimz,m,lambdaP,sigma,iterP,None, None)				
	py_printErrorMessage(err_code)
	return out_image

# BilateralFilter
def py_p3dBilateralFilter8(image_data, dimx, dimy, dimz = 0, size = 3, sigma_d = 1.0, sigma_r = 3, iterations = 10):
	"""
 
  Bilateral filtering smooths images while preserving edges, by means of a nonlinear combination of nearby image values..
 
 The idea uderlying bilateral filtering is to do in the range of an image what traditional filter do in its domain. Two pixels  (or voxels) can be close to one another, that is, occupy nearby spatial location, or they can be similar to one another, that  is, have nearby values. Closeness refers to the vicinity in the domain, similarity to vicinity in the range. Traditional  filtering is domain filtering and enforces closeness by weighing pixel (or voxel) values with coefficients that fall off with  distance. Similarly, range filters average image values with weights that decay with dissimilarity. Range filters are nonlinear  because their weights depend on image intensity. The combination of domain and range filtering is denoted as bilateral  filtering [1].
  
  Syntax:
  ------
  Result = py_p3dBilateralFilter8 ( image_data, dimx, dimy [, dimz = value] [, size=value ] [,sigma_d =value ] [, sigma_r=value ] [,  iterations=value ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input images.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
  
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  size: An odd integer value in the range [3, 51] (default: 3) representing the side of the cubic kernel.
  
  sigma_d: A decimal value greater than 0 (default: 1.0) representing the standard deviation of the domain gaussian.
  
  sigma_r: A decimal value greater than 0 (default: 3.0) representing the standard deviation of the range gaussian.
  
  iterations: The number of times the algorithm is iterated (default: 10).
  
  Remarks:
  -------
 Often a 3x3 square kernel (or a 3x3x3 cubic kernel for 3D images) is used, although larger kernels (e.g. 5x5 squares or 5x5x5  cubes) can be used for more severe smoothing. Note that a small kernel can be applied more than once in order to produce a  similar but not identical effect as a single pass with a large kernel.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	if size < 3 or size > 51 or (size % 2) == 0:
		py_printErrorMessage(-4)
		return   

	if sigma_d < 0:
		py_printErrorMessage(-4)
		return   
    
	if sigma_r < 0:
		py_printErrorMessage(-4)
		return 

	if iterations < 0:
		py_printErrorMessage(-4)
		return 

	out_image =  malloc_uchar(dimx*dimy*dimz)
	err_code = p3dBilateralFilter3D_8(image_data,out_image, dimx,dimy,dimz,size,sigma_d,sigma_r,iter,None, None)				
	py_printErrorMessage(err_code)
	return out_image
	
# BoinHaibelRingRemover 2D only
def py_p3dBoinHaibelRingRemover8(image_data, dimx, dimy, centerX=0, centerY=0, winsize=0, iterations=0, precision=0):
	"""
  
  Syntax:
  ------
  Result = py_p3dBoinHaibelRingRemover8 ( image_data, dimx, dimy [, centerX=value ] [,centerY =value ] [, winsize=value ]     [,  iterations=value ] )
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input images.  
 
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  centerX,centerY:  
 
  winsize:  
 
  iterations: 
  
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return
    
	out_image = malloc_uchar(dimx*dimy)
	err_code = p3dBoinHaibelRingRemover2D_8(image_data,out_image, dimx,dimy,centerX,centerY,winsize,iterations, precision,None)
	py_printErrorMessage(err_code)
	return out_image

# SijbersPostnovRingRemover 2D only 
def py_p3dSijbersPostnovRingRemover8(image_data, dimx, dimy, centerX = 0, centerY = 0, winsize = 5, thresh = 1, iterations = 1, precision = 1.5, mask = None):
	"""
 
  Reduce ring artifacts from reconstructed images using a modified J. Sijbers and A. Postnov algorithm [1].
 
 The method is based on the observation that ring artifacts become straight vertical lines trasforming input image in polar  coordinates where the center of the ring artifacts is assumed as the center of polar transformation. Within a sliding window a  set of homogeneous rows of polar image is detected. Working on this set an artifact template is generated and used for the  correction of the image. At the end, the image is transformed back into cartesian coordinates.
  
  Syntax:
  ------
  Result = py_p3dSijbersPostnovRingRemover8 ( image_data, dimx, dimy [, mask= 2D image ] [, centerX=value ] [,centerY =value ] [, winsize=value  ] [, thresh=value ][, iterations=1 ][, precision=value ] )
  
  Return Value:
  ------------
  Returns an image with the same dimensions and type of input image having ring artifacts removed.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  centerX,centerY: Coordinates of the center of ring artifacts. If this keyword is not setted the image center will be assumed  as ring artifacts center (optimal condition). 
  
  winsize: An odd integer value in the range [3, 201] (default: 51). representing the size of the sliding window used to scan  the image. The purpose of the window is to detect homogeneous row segments. Therefore, the width should be chosen of adequate  size in order to detect a sufficiently large number of homogeneous rows.
 
  tbresh: A decimal value greater than 0 representing the threshold value that classifies row segments as homogeneous or  inhomogeneous.
 The choice of Thresh depends on the severity of the line artifacts. The less pronounced the line artifacts, the smaller the  value of threshold can be chosen, with a lower bound determined by the image noise variance. A good choice for Thresh is around  3 times the image noise variance (see further remarks for the estimation of noise variance).  
 
  iterations: Filter could be applied iteratively the specified number of times (default: 1) greater than 0. It is equivalent to a iterative  invocation of the filter but execution time is lower (filter remains in polar coordinates and performs the conversion to  cartesian coordinates only at the end of last iteration).
 
  precision: The discrete polar-cartesian conversion could require a denser grid in order to avoid compromising structures far  from ring artifacts center. Increasing the precision parameter (between 1.0 and 5.0) results in a more precise processing but  more computational time is required. If no parameter is specified, the fastest filtering is performed (default: 1.5).
 
   mask :A 2D matrix (coherent with Image argument) of type BYTE having value 255 on pixel or voxel representing the region of interest and 0 elsewhere. Filtering will be applied only within the specified Mask resulting in a faster filter execution with better results in the vicinity of mask boundary. Mask could be omitted for objects having circular simmetry with baricenter that fits neatly the center of ring artifacts. Mask could be omitted also for less pronounced artifacts and for high porous objects.

  
  Remarks:
  -------
  Trying to read a RAW image stored in SIGNED format is in principle an error. 
  Also, using p3dReadRaw8 for a 16-bit RAW image results in a completely wrong output but, again, no error could be detected.
  Attention #must be paid when using this function. Pore3D internally deals only with UNSIGNED format.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return
    
	if winsize < 3 or winsize > 201 or (winsize % 2) == 0:
		py_printErrorMessage(-4)
		return   

	if tbresh < 0:
		py_printErrorMessage(-4)
		return 

	if iterations < 0:
		py_printErrorMessage(-4)
		return 

	if precision < 1.0 or precision > 5.0:
		py_printErrorMessage(-4)
		return   
    
	if centerX == 0:
		centerX = dimx/2
	if centerY == 0:
		centerY = dimy/2
	out_image =  malloc_uchar(dimx*dimy)
	err_code = p3dSijbersPostnovRingRemover2D_8(image_data,out_image, dimx,dimy,centerX,centerY,winsize,thresh,iterations, precision, mask, None,None)				
	py_printErrorMessage(err_code)
	return out_image

# ClearBorderFilter
def py_p3dClearBorderFilter8(image_data, dimx, dimy, dimz = 0, connIn = 6):
	"""
 
  This filter suppresses connected components (blobs) connected to image margins.
 
 The algorithm starts performing a modified connected components labeling (see py_p3dBlobLabeling) starting from image margins.  Each connected components found is "labeled" with background value, i.e. it is removed from image.
  
  Syntax:
  ------
  Result = py_p3dClearBorderFilter8 (image_data, dimx, dimy [, dimz = value][,connIn = value])
  
  Return Value:
  ------------
  Returns a filtered image with the same dimensions and type of input image.
  
  Arguments:
  ---------
  image_data: A 2D or 3D matrix of type BYTE representing the binary input image to filter.
 
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
 
  connIn: The desired connectivity, i.e. 6, 18 or 26 for 3D images. The default connectivity is 6 for three dimensional images.
  
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return
    
	out_image =  malloc_uchar(dimx*dimy*dimz)
	conn3D = 711
	if connIn == 6:
		conn3D = 711       
	elif connIn == 18:
		conn3D = 712      
	elif connIn == 26:
		conn3D = 713     
	else:
		py_printErrorMessage(-4)


        
	err_code = p3dClearBorderFilter3D(image_data,out_image, dimx,dimy,dimz,conn3D,None)
	py_printErrorMessage(err_code)
	return out_image
	
# AutoThresholding
def py_p3dAutoThresholding8(image_data, dimx, dimy, dimz = 0, methodNum = 1, thresh_val = 7):
	"""
 
  Performs segmentation by thresholding using an automatically determined threshold according to one of the following methods:
 
  1. Kittler and Illingworth: The method consists in arbitrarily dividing the histogram into two parts (the foreground and the  background), modeling each part with a normal distribution, comparing the resulting model based on a mixture of these normal  distribution with the original histogram and assuming as optimal the threshold that minimizes a criterion function based on the  classification error probability [1].
 
  2. Otsu: The intra-class variance is computed, i.e., the weighted sum of the variances of each class (the background and the  foreground) adopting the number of voxels in the class as a weight. Among all the possible thresholds, the optimal is the one  that minimizes this intra-class variance [2].
 
  3.Pun: This method is based on entropic thresholding: an entropy-thresholded image is the one that preserves (as much as  possible) the information contained in the original unthresholded image in terms of entropy [3].
 
  4.Ridler and Calvard: A unique threshold is assumed to be the average of the foreground and background class means. The means  of the two parts can be evaluated only after the threshold is determined, but the threshold needs to be computed from the two  means. Therefore, an iterative algorithm was suggested: first, an initial threshold is selected (the mean of the entire  histogram is a sufficient starting point), then the two means for the two distributions on either side of the threshold are  calculated. A new threshold is obtained by averaging these means and the process continues until the value of the threshold  converges [4].
 
  5.Kapur et al.Kapur et al. improved the Pun's approach considering the image foreground and background as two different  classes of events. It first measures the class entropies, which is interpreted as a measure of class compactness. When the sum  of the two class entropies reaches its maximum, the image is said to be optimally thresholded [5].
 
 6.Tsai: The Tsai's method first computes gray-level moments from the input's histogram, and then obtains the threshold  according to the principle [6].
  
  Syntax:
  ------
  Result = py_p3dAutoThresholding8 ( image_data, dimx, dimy [, dimz = value] [, methodNum = value] [, thresh_val = variable])
  
  Return Value:
  ------------
  Returns a volume with the same dimensions of input volume having value 255 on skeleton voxels and 0 elsewhere.
  
  Arguments:
  ---------
  image_data: A string with the filename with full path included.
  
  dimx,dimy,dimz: three variables representing the dimensions of image to read. 
  
  methodNum: The method (default: 1) for the automatic determination of the threshold.
  
  thresh_val: An output parameter that will contain the automatically computed threshold. This parameter will be in the range 
 [0, 255] for a byte image.
  
  References
  ---------
  [1] Kittler J, Illingworth J. Minimum error thresholding. Pattern Recogn. 1986;19:41-7.
  [2] Otsu N. A threshold selection method from gray-level histograms. IEEE Trans Syst Man Cybern. 1979;9:62-6.
  [3] Kapur JN, Sahoo PK, Wong AKC. A new method for gray-level picture thresholding using the entropy of the histogram. Graph    Models Image Process. 1985;29:273-85.
  [4] Tsai WH. Moment-preserving thresholding: a new approach. Graph Models Image Process. 1985;19:377-93.
  [5] Ridler TW, Calvard S. Picture thresholding using an iterative selection method. IEEE Trans Syst Man Cybern. 1978;8:630-2.
  [6] Pun T. Entropic thresholding: a new approach. Comp Graph Image Process. 1981;16:210-39.
 
	"""
	if dimx == 0 or dimy == 0:
		py_printErrorMessage(-3)
		return

	if thresh_val < 6 or thresh_val > 255:
		py_printErrorMessage(-4)
		return

	thresh_val = 0
	out_image =  malloc_uchar(dimx*dimy*dimz)
	thresh_val = malloc_uchar(1)
	if methodNum == 1:
		err_code = p3dKittlerThresholding_8(image_data, out_image,dimx,dimy,dimz, thresh_val, None, None)
	elif methodNum == 2:
		err_code = p3dOtsuThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)	
	elif methodNum == 3:
		err_code = p3dPunThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)
	elif methodNum == 4:
		err_code = p3dRidlerThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)
	elif methodNum == 5:
		err_code = p3dKapurThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)
	elif methodNum == 6:
		err_code = p3dJohannsenThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)
	elif methodNum == 7:
		err_code = p3dHuangYagerThresholding_8(image_data, out_image,dimx,dimy,dimz,thresh_val, None, None)
	else:
		err_code = p3dOtsuThresholding_8(image_data, out_image,dimx,dimy,dimz, thresh_val, None, None)
	
	py_printErrorMessage(err_code)
	return out_image
		
		