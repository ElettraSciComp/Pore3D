#ifdef __cplusplus
extern "C" {
#endif


    /*
            Constants:
     */
#ifndef P3D_ITKWRAPPED_DEFINED
#define P3D_ITKWRAPPED_DEFINED

#define P3D_FALSE				-1 
#define P3D_TRUE				1 

#define P3D_AUTH_ERROR                          -1
#define P3D_ERROR                               NULL	/* Leave it NULL for simplify tests */
#define P3D_SUCCESS				2	/* Any number */

#define BACKGROUND				0
#define OBJECT					UCHAR_MAX	

    
#define CONN4   611
#define CONN8   612

    
    // Constants for 3D connectivity:
#define CONN6   711
#define CONN18  712
#define CONN26  713

#endif
    
    // Filters:
    int p3dGradientAnisotropicDiffusionFilter2D_8( unsigned char*, unsigned char*, const int,  const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dGradientAnisotropicDiffusionFilter2D_16( unsigned short*, unsigned short*, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dGradientAnisotropicDiffusionFilter3D_8( unsigned char*, unsigned char*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dGradientAnisotropicDiffusionFilter3D_16( unsigned short*, unsigned short*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    
    int p3dCurvatureAnisotropicDiffusionFilter2D_8( unsigned char*, unsigned char*, const int,  const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dCurvatureAnisotropicDiffusionFilter2D_16( unsigned short*, unsigned short*, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dCurvatureAnisotropicDiffusionFilter3D_8( unsigned char*, unsigned char*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    int p3dCurvatureAnisotropicDiffusionFilter3D_16( unsigned short*, unsigned short*, const int, const int, const int, const int, const double, int (*wr_log)(const char*, ...));
    
    int p3dCurvatureFlowFilter2D_8( unsigned char*, unsigned char*, const int,  const int, const int, int (*wr_log)(const char*, ...));
    int p3dCurvatureFlowFilter2D_16( unsigned short*, unsigned short*, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dCurvatureFlowFilter3D_8( unsigned char*, unsigned char*, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dCurvatureFlowFilter3D_16( unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    
    int p3dMinMaxCurvatureFlowFilter2D_8( unsigned char*, unsigned char*, const int,  const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dMinMaxCurvatureFlowFilter2D_16( unsigned short*, unsigned short*, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dMinMaxCurvatureFlowFilter3D_8( unsigned char*, unsigned char*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
    int p3dMinMaxCurvatureFlowFilter3D_16( unsigned short*, unsigned short*, const int, const int, const int, const int, const int, int (*wr_log)(const char*, ...));
  
    
    // Segmentation:
    int p3dConfidenceConnectedRegionGrowing2D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const double,	const unsigned int, const unsigned int,	const unsigned int, const unsigned int,	int (*wr_log)( const char*, ... ));
    int p3dConfidenceConnectedRegionGrowing2D_16 ( unsigned short*, unsigned char*, const unsigned int,const unsigned int, const double, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    int p3dConfidenceConnectedRegionGrowing3D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const double, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    int p3dConfidenceConnectedRegionGrowing3D_16 ( unsigned short*, unsigned char*, const unsigned int,	const unsigned int, const unsigned int,	const double, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    
    int p3dMultipleOtsuThresholding2D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    int p3dMultipleOtsuThresholding2D_16 ( unsigned short*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    int p3dMultipleOtsuThresholding3D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    int p3dMultipleOtsuThresholding3D_16 ( unsigned short*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)( const char*, ... ));
    
    // Morphology:
    int p3dBinaryDilateFilter2D( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)(const char*, ...));
    int p3dBinaryDilateFilter3D( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)(const char*, ...));
        
    int p3dBinaryErodeFilter2D( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)(const char*, ...));
    int p3dBinaryErodeFilter3D( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, int (*wr_log)(const char*, ...));
    
    int p3dHMinimaFilter2D_8 (unsigned char*,  unsigned char*, const unsigned int, const unsigned int, const int, const unsigned char, int (*wr_log)(const char*, ...));
    int p3dHMinimaFilter2D_16 (unsigned short*, unsigned short*, const unsigned int, const unsigned int, const int, const unsigned short, int (*wr_log)(const char*, ...));
    int p3dHMinimaFilter3D_8 (unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const int, const unsigned char, int (*wr_log)(const char*, ...));
    int p3dHMinimaFilter3D_16 (unsigned short*, unsigned short*, const unsigned int, const unsigned int, const unsigned int, const int, const unsigned short, int (*wr_log)(const char*, ...));
    
    int p3dWatershedSegmentation2D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const int, int (*wr_log)(const char*, ...));
    int p3dWatershedSegmentation2D_16 (	unsigned short*, unsigned char*, const unsigned int, const unsigned int, const int, int (*wr_log)(const char*, ...));    
    int p3dWatershedSegmentation3D_8 ( unsigned char*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const int, int (*wr_log)(const char*, ...));
    int p3dWatershedSegmentation3D_16 (	unsigned short*, unsigned char*, const unsigned int, const unsigned int, const unsigned int, const int,	int (*wr_log)(const char*, ...));
	    
#ifdef __cplusplus
}
#endif
