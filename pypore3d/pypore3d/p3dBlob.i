
%module p3dBlob

%include typemaps.i
%include "exception.i"


%{
   #define SWIG_FILE_WITH_INIT
   
#include <omp.h>    

#include "P3D_Blob/p3dBlob.h"
#include "P3D_Blob/p3dTime.h"
#include "P3D_Blob/Common/p3dBoundingBoxList.h"
#include "P3D_Blob/Common/p3dBoundingBoxT.h"
#include "P3D_Blob/Common/p3dConnectedComponentsLabeling.h"
#include "P3D_Blob/Common/p3dCoordsList.h"
#include "P3D_Blob/Common/p3dCoordsQueue.h"
#include "P3D_Blob/Common/p3dCoordsT.h"
#include "P3D_Blob/Common/p3dDoubleList.h"
#include "P3D_Blob/Common/p3dFCoordsList.h"
//#include "P3D_Blob/Common/p3dSquaredEuclideanDT.h"
#include "P3D_Blob/Common/p3dUIntList.h"
#include "P3D_Blob/Common/p3dUtils.h"

#include <assert.h>
int myErr = 0; // flag to save error state

%}


  

////////////////
%include <cmalloc.i>
%include <cdata.i>
%allocators(unsigned char, uchar)
%allocators(unsigned short, ushort)
%allocators(double, doub)

%allocators(BasicStats, BasicStat)
%allocators(AnisotropyStats, AnisotropyStat)
%allocators(MorphometricStats, MorphometricStat)
%allocators(TextureStats, TextureStat)
%allocators(BlobStats, BlobStat)



////////////////
%inline %{

//#include "P3D_Blob/p3dBlob.h"
#include <float.h>

void PrintBlobStruct(BlobStats* out_stats, char* filename)
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
    
    //fprintf(fptr, "Count:%d\n", out_stats->blobCount);
    fprintf(fptr, "VOLUME,MAX_SPHERE,EQ_SPHERE,MIN_AXIS,MAX_AXIS,SPHERICITY,ASPECT_RATIO,Extent\n"); 
    for (int i = 0; i<out_stats->blobCount; i++)
       {
       
            // adjust sphericity and aspect ration
            //float sphericity = 0.0;
            //if (out_stats->eq_sph[i]>=0.00000001)
            //    sphericity = out_stats->max_sph[i]/out_stats->eq_sph[i];
            //out_stats->sphericity[i] = sphericity;
                
            //float aspect_ratio = 0.0;
            //if (out_stats->l_max[i]>=0.00000001)
            //    aspect_ratio  = out_stats->l_min[i]/out_stats->l_max[i];
            //out_stats->aspect_ratio[i] = aspect_ratio ;
                
            //print to file
            fprintf(fptr, "%.*e,", Digs, out_stats->volume[i]);
            fprintf(fptr, "%.*e,", Digs, out_stats->max_sph[i]);
            fprintf(fptr, "%.*e,", Digs, out_stats->eq_sph[i]);
            fprintf(fptr, "%.*e,", Digs, out_stats->l_min[i]);
            fprintf(fptr, "%.*e,", Digs,out_stats->l_max[i]);
            fprintf(fptr, "%.*e,", Digs, out_stats->sphericity[i]);
            fprintf(fptr, "%.*e,", Digs, out_stats->aspect_ratio[i]);
            fprintf(fptr, "%.*e\n", Digs, out_stats->extent[i]);
            
                
       }
   
       
    fclose(fptr);
}


void invert_vol(unsigned char* in_im, int dimx, int dimy, int dimz)
{
       for (int k = 0; k < dimz; k++)
               for (int j = 0; j < dimy; j++)
                       for (int i = 0; i < dimx; i++)
                               {
                               
                                       int ind = (i*dimz*dimy) + (j*dimz) + k;
                                       in_im[ind] = 255 -in_im[ind];
                               /*        if (in_im[ind] >0) 
                                               in_im[ind] = 0; 
                                       else 
                                                in_im[ind] = 255 ;
                                                */
                               }               
}


void invert_vol_16(unsigned short* in_im, int dimx, int dimy, int dimz)
{
       for (int k = 0; k < dimz; k++)
               for (int j = 0; j < dimy; j++)
                       for (int i = 0; i < dimx; i++)
                               {
                               	int ind = (i*dimz*dimy) + (j*dimz) + k;
                                       in_im[ind] = (unsigned short) (255 -in_im[ind]);
                                       /*
                                       int ind = (i*dimz*dimy) + (j*dimz) + k;
                                       if (in_im[ind] >0.0) 
                                               in_im[ind] = 0.0; 
                                       else 
                                                in_im[ind] = 1.0 ;
                                                */
                               }               
}


%}



////////////////



////////////////
%include "P3D_Blob/p3dBlob.h"
%include "P3D_Blob/p3dTime.h"
%include "P3D_Blob/Common/p3dBoundingBoxList.h"
%include "P3D_Blob/Common/p3dBoundingBoxT.h"
%include "P3D_Blob/Common/p3dConnectedComponentsLabeling.h"
%include "P3D_Blob/Common/p3dCoordsList.h"
%include "P3D_Blob/Common/p3dCoordsQueue.h"
%include "P3D_Blob/Common/p3dCoordsT.h"
%include "P3D_Blob/Common/p3dDoubleList.h"
%include "P3D_Blob/Common/p3dFCoordsList.h"
/*%include "P3D_Blob/Common/p3dSquaredEuclideanDT.h"*/
%include "P3D_Blob/Common/p3dUIntList.h"
%include "P3D_Blob/Common/p3dUtils.h"

 


