import os
import setuptools
from setuptools import setup
from setuptools.extension import Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

this_dir = os.path.dirname(os.path.abspath(__file__))


ext_modules = Extension(
    name='pypore3d._p3dBlob',
    include_dirs=[os.path.join(this_dir, 'pypore3d')],
    sources=["pypore3d/p3dBlob.i",
    "pypore3d/P3D_Blob/_p3dTime.c",
    "pypore3d/P3D_Blob/p3dBasicAnalysis.c",
    "pypore3d/P3D_Blob/p3dAnisotropyAnalysis.c",
    "pypore3d/P3D_Blob/p3dBlobAnalysis.c" ,
    "pypore3d/P3D_Blob/p3dBlobLabeling.c", 
    "pypore3d/P3D_Blob/p3dChamferDT.c", 
    "pypore3d/P3D_Blob/p3dGetMaxVolumeBlob.c", 
    "pypore3d/P3D_Blob/p3dGetMinVolumeBlob.c", 
    "pypore3d/P3D_Blob/p3dMinVolumeFilter.c", 
    "pypore3d/P3D_Blob/p3dMorphometricAnalysis.c", 
    "pypore3d/P3D_Blob/p3dREVEstimation.c", 
    "pypore3d/P3D_Blob/p3dTextureAnalysis.c", 
    "pypore3d/P3D_Blob/Common/p3dBoundingBoxList.c", 
    "pypore3d/P3D_Blob/Common/p3dConnectedComponentsLabeling_uint.c", 
    "pypore3d/P3D_Blob/Common/p3dConnectedComponentsLabeling_ushort.c", 
    "pypore3d/P3D_Blob/Common/p3dCoordsList.c", 
    "pypore3d/P3D_Blob/Common/p3dCoordsQueue.c", 
    "pypore3d/P3D_Blob/Common/p3dDoubleList.c", 
    "pypore3d/P3D_Blob/Common/p3dFCoordsList.c", 
    "pypore3d/P3D_Blob/p3dSquaredEuclideanDT.c", 
    "pypore3d/P3D_Blob/Common/p3dUIntList.c", 
    "pypore3d/P3D_Blob/Common/p3dUtils.c" 

]) 


ext_modules1 = Extension(
    name='pypore3d._p3dSkel',
    include_dirs=[os.path.join(this_dir, 'pypore3d')],
    sources=["pypore3d/p3dSkel.i",
    "pypore3d/P3D_Skel/Common/p3dBoundingBoxList.c",
    "pypore3d/P3D_Skel/Common/p3dConnectedComponentsLabeling.c",
    "pypore3d/P3D_Skel/Common/p3dCoordsList.c",
    "pypore3d/P3D_Skel/Common/p3dCoordsQueue.c",
    "pypore3d/P3D_Skel/Common/p3dFCoordsList.c",
    "pypore3d/P3D_Skel/Common/p3dSquaredEuclideanDT.c",
    "pypore3d/P3D_Skel/Common/p3dThinning.c",
    "pypore3d/P3D_Skel/Common/p3dUIntList.c",
    "pypore3d/P3D_Skel/Common/p3dUtils.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dComputeCoreSkeleton.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dComputeEigenVal.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dComputeHierarchicalSkeleton.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dCriticalPoints.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dCritPointList.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dGetHighDivPoints.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dGVF.c",
    "pypore3d/P3D_Skel/GVFSkeletonization/p3dHighDivPointList.c",
    "pypore3d/P3D_Skel/_p3dTime.c",
    "pypore3d/P3D_Skel/p3dGVFSkeletonization.c",
    "pypore3d/P3D_Skel/p3dIterativeSkeletonPruning.c",
    "pypore3d/P3D_Skel/p3dLKCSkeletonization.c",
    "pypore3d/P3D_Skel/p3dSimpleSkeletonPruning.c",
    "pypore3d/P3D_Skel/p3dSkeletonAnalysis.c",
    "pypore3d/P3D_Skel/p3dSkeletonAnalysisFeasibility.c",
    "pypore3d/P3D_Skel/p3dSkeletonLabeling.c",
    "pypore3d/P3D_Skel/p3dThinningSkeletonization.c",
    "pypore3d/P3D_Skel/p3dUltimateSkeletonPruning.c"]) 


ext_modules2 = Extension(
    name='pypore3d._p3dFilt',
    include_dirs=[os.path.join(this_dir, 'pypore3d')],
    sources=["pypore3d/p3dFilt.i",
    "pypore3d/P3D_Filt/p3dMedianFilter.c", 
    "pypore3d/P3D_Filt/_p3dTime.c",
    "pypore3d/P3D_Filt/p3dMeanFilter.c", 
    "pypore3d/P3D_Filt/p3dBilateralFilter.c", 
    "pypore3d/P3D_Filt/p3dIORaw.c",
    "pypore3d/P3D_Filt/p3dAnisotropicDiffusionFilter.c",
    "pypore3d/P3D_Filt/p3dBoinHaibelRingRemover.c",
    "pypore3d/P3D_Filt/p3dClearBorderFilter.c",
    "pypore3d/P3D_Filt/p3dCreateBinaryShapes.c",
    "pypore3d/P3D_Filt/p3dCrop.c",
    "pypore3d/P3D_Filt/p3dFrom16To8.c",
    "pypore3d/P3D_Filt/p3dGaussianFilter.c",
    "pypore3d/P3D_Filt/p3dGetRegionByCoords.c",
    "pypore3d/P3D_Filt/p3dHuangYagerThresholding.c",
    "pypore3d/P3D_Filt/p3dJohannsenThresholding.c",
    "pypore3d/P3D_Filt/p3dKapurThresholding.c",
    "pypore3d/P3D_Filt/p3dKittlerThresholding.c",
    "pypore3d/P3D_Filt/p3dOtsuThresholding.c",
    "pypore3d/P3D_Filt/p3dPadding.c",
    "pypore3d/P3D_Filt/p3dPunThresholding.c",
    "pypore3d/P3D_Filt/p3dRidlerThresholding.c",
    "pypore3d/P3D_Filt/p3dSijbersPostnovRingRemover.c", 
    "pypore3d/P3D_Filt/Common/p3dRingRemoverCommon.c",
    "pypore3d/P3D_Filt/Common/p3dCoordsQueue.c"]) 

setup(
name = 'pypore3d',
version = '0.0.1',
author = 'Amal Aboulhassan',
author_email = "aboulhassan.amal@gmail.com",
description = "Pore3d with python wrappers",
long_description = long_description,
long_description_content_type = "text/markdown",
url = "https://gitlab.elettra.eu/amal.abouelhassan/pore3d_py",
packages=setuptools.find_packages(),
classifiers =
    ["Programming Language :: Python :: 3",
     "License :: OSI Approved :: MIT License",
     "Operating System :: OS Independent"],
ext_modules=[ext_modules, ext_modules1, ext_modules2]
)


#package_dir = 'src',
#ext_modules=[src],
#zip_safe=False,
#include_package_data=True,