/* File: p3dMedianFilter.i */
%module p3dFilt

%include typemaps.i

%{
   #define SWIG_FILE_WITH_INIT
   
#include <omp.h>    
/*#include "p3dTime.h"*/
#include "P3D_Filt/p3dTime.h"   
#include "P3D_Filt/p3dFilt.h"
#include "P3D_Filt/Common/p3dCoordsQueue.h"
#include "P3D_Filt/Common/p3dCoordsT.h"
#include "P3D_Filt/Common/p3dRingRemoverCommon.h"

%}

%include <cmalloc.i>
%include <cdata.i>

%allocators(unsigned char, uchar)
%allocators(unsigned short, ushort)

/*#include "p3dTime.h"*/
%include "P3D_Filt/p3dTime.h"   
%include "P3D_Filt/p3dFilt.h"
%include "P3D_Filt/Common/p3dCoordsQueue.h"
%include "P3D_Filt/Common/p3dCoordsT.h"
%include "P3D_Filt/Common/p3dRingRemoverCommon.h"




