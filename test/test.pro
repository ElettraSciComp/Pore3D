CD, 'C:\Users\Franz\Documents\Visual Studio 2012\Projects\Pore3D\test'

im = p3dReadRaw8('test_274x274x326.raw',[274,274,326])

flt1 = p3dCurvatureAnisotropicDiffusionFilter(im)

th = p3dAutoThresholding(flt1)

th = p3dMinVolumeFilter(th)

blb = p3dBlobLabeling(th)

stat = p3dBlobAnalysis(th)

th = BYTE(255-th)

stat2 = p3dLKCSkeletonization(th)

end