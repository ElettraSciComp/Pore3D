/* File: p3dMedianFilter.i */
%module p3dSkel

%include typemaps.i

%{
   #define SWIG_FILE_WITH_INIT
   
#include "P3D_Skel/p3dSkel.h"
#include "P3D_Skel/p3dTime.h"
#include "P3D_Skel/Common/p3dBoundingBoxT.h"
#include "P3D_Skel/Common/p3dConnectedComponentsLabeling.h"
#include "P3D_Skel/Common/p3dCoordsList.h"
#include "P3D_Skel/Common/p3dCoordsQueue.h"
#include "P3D_Skel/Common/p3dCoordsT.h"
#include "P3D_Skel/Common/p3dFCoordsList.h"
#include "P3D_Skel/Common/p3dSquaredEuclideanDT.h"
#include "P3D_Skel/Common/p3dThinning.h"
#include "P3D_Skel/Common/p3dUIntList.h"
#include "P3D_Skel/GVFSkeletonization/p3dComputeCoreSkeleton.h"
#include "P3D_Skel/GVFSkeletonization/p3dComputeEigenVal.h"
#include "P3D_Skel/GVFSkeletonization/p3dComputeHierarchicalSkeleton.h"
#include "P3D_Skel/GVFSkeletonization/p3dCriticalPoints.h"
#include "P3D_Skel/GVFSkeletonization/p3dCritPointList.h"
#include "P3D_Skel/GVFSkeletonization/p3dCritPointT.h"
#include "P3D_Skel/GVFSkeletonization/p3dGetHighDivPoints.h"
#include "P3D_Skel/GVFSkeletonization/p3dGVF.h"
#include "P3D_Skel/GVFSkeletonization/p3dHighDivPointList.h"
#include "P3D_Skel/GVFSkeletonization/p3dHighDivPointT.h"
 
 
%}

%include <cmalloc.i>
%include <cdata.i>

%allocators(unsigned char, uchar)
%allocators(unsigned short, ushort)
%allocators(struct SkeletonStats, PSkeletonStats)

%inline %{

    #include <float.h>
    void PrintSkelStruct(struct SkeletonStats* out_stats, char* filename)
    {
        // creating file pointer to work with files
        FILE *fptr;

        // opening file in writing mode
        fptr = fopen(filename, "w");

        // exiting program 
        if (fptr == NULL) {
            printf("Error!");
            exit(1);
        }

        int Digs = FLT_DECIMAL_DIG;

        fprintf(fptr, "End_Width\n"); 
        for (int i = 0; i<out_stats->End_Counter; i++)
           {
                //print to file
                fprintf(fptr, "%.*e\n", Digs, out_stats->End_Width[i]);
           }


        fprintf(fptr, "EndToEnd_Length,EndToEnd_MaxWidth,EndToEnd_MeanWidth,EndToEnd_MinWidth\n"); 
        for (int i = 0; i<out_stats->EndToEnd_Counter; i++)
           {
                //print to file
                fprintf(fptr, "%.*e,", Digs, out_stats->EndToEnd_Length[i]);
                fprintf(fptr, "%.*e,", Digs, out_stats->EndToEnd_MaxWidth[i]);
                fprintf(fptr, "%.*e,", Digs, out_stats->EndToEnd_MeanWidth[i]);
                fprintf(fptr, "%.*e\n", Digs, out_stats->EndToEnd_MinWidth[i]);
           }
        
        fprintf(fptr, "NodeToNode_Length,NodeToNode_MaxWidth,NodeToNode_MeanWidth,NodeToNode_MinWidth\n"); 
        for (int i = 0; i<out_stats->NodeToNode_Counter; i++)
           {
                //print to file
                fprintf(fptr, "%.*e,", Digs, out_stats->NodeToNode_Length[i]);
                fprintf(fptr, "%.*e,", Digs, out_stats->NodeToNode_MaxWidth[i]);
                fprintf(fptr, "%.*e,", Digs, out_stats->NodeToNode_MeanWidth[i]);
                fprintf(fptr, "%.*e\n", Digs, out_stats->NodeToNode_MinWidth[i]);
           }
        
        
        fprintf(fptr, "Node_Width\n"); 
        for (int i = 0; i<out_stats->Node_Counter; i++)
           {
                //print to file
                fprintf(fptr, "%.*e\n", Digs, out_stats->Node_Width[i]);
           }


        fclose(fptr);
        }
%}



%include "P3D_Skel/p3dSkel.h"
%include "P3D_Skel/p3dTime.h"
%include "P3D_Skel/Common/p3dBoundingBoxT.h"
%include "P3D_Skel/Common/p3dConnectedComponentsLabeling.h"
%include "P3D_Skel/Common/p3dCoordsList.h"
%include "P3D_Skel/Common/p3dCoordsQueue.h"
%include "P3D_Skel/Common/p3dCoordsT.h"
%include "P3D_Skel/Common/p3dFCoordsList.h"
%include "P3D_Skel/Common/p3dSquaredEuclideanDT.h"
%include "P3D_Skel/Common/p3dThinning.h"
%include "P3D_Skel/Common/p3dUIntList.h"
%include "P3D_Skel/GVFSkeletonization/p3dComputeCoreSkeleton.h"
%include "P3D_Skel/GVFSkeletonization/p3dComputeEigenVal.h"
%include "P3D_Skel/GVFSkeletonization/p3dComputeHierarchicalSkeleton.h"
%include "P3D_Skel/GVFSkeletonization/p3dCriticalPoints.h"
%include "P3D_Skel/GVFSkeletonization/p3dCritPointList.h"
%include "P3D_Skel/GVFSkeletonization/p3dCritPointT.h"
%include "P3D_Skel/GVFSkeletonization/p3dGetHighDivPoints.h"
%include "P3D_Skel/GVFSkeletonization/p3dGVF.h"
%include "P3D_Skel/GVFSkeletonization/p3dHighDivPointList.h"
%include "P3D_Skel/GVFSkeletonization/p3dHighDivPointT.h"
