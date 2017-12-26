#ifdef __cplusplus
extern "C" {
#endif


	/*
	Constants:
	*/
#ifndef P3D_SKEL_DEFINED  
	#define P3D_SKEL_DEFINED

	#define P3D_FALSE				-1 
	#define P3D_TRUE				1 

	#define P3D_ERROR				0
	#define P3D_MEM_ERROR           NULL	/* Leave it NULL for simplify tests */
	#define P3D_SUCCESS				2		/* Any number */

	#define BACKGROUND				0
	#define OBJECT					UCHAR_MAX	

		// Skeleton labels on OUT skeleton:
	#define NODE_LABEL        1
	#define END_LABEL         2
	#define ISOLATED_LABEL    3
	#define NODETONODE_LABEL  4
	#define NODETOEND_LABEL   5
	#define ENDTOEND_LABEL    6

#endif

	/*
	Macros:
	*/

#ifndef P3D_SKEL_MACROS
	#define P3D_SKEL_MACROS

	#define I(i,j,k,N,M)    ( (j)*(N) + (i) + (k)*(N)*(M) ) 
	#define I2(i,j,N)       ( (j)*(N) + (i) ) 
	#define MIN(x,y)        (((x) < (y))?(x):(y))
	#define MAX(x,y)        (((x) > (y))?(x):(y))

		/* A sort of TRY-CATCH constructor: */
	#define P3D_MEM_TRY( function ) if ( ((function) == P3D_MEM_ERROR) ) { goto MEM_ERROR; }
	#define P3D_TRY( function ) if ( ((function) == P3D_ERROR) ) { goto MEM_ERROR; }

#endif

#ifndef P3D_SKEL_STRUCTS_DEFINED
#define P3D_SKEL_STRUCTS_DEFINED  

	// Struct for results of skeleton analysis:
	struct SkeletonStats {

		// Connectivity indexes:
		double   ConnectivityDensity;
		int*     CoordinationNumber;	             

		// END points:
		unsigned short  End_Counter;	
		double*  End_Width; 

		// END-TO-END branches:
		unsigned short EndToEnd_Counter;
		double*  EndToEnd_Length;
		double*  EndToEnd_MaxWidth;    // Max thickness
		double*  EndToEnd_MeanWidth;   // Mean thickness
		double*  EndToEnd_MinWidth;    // Throat size

		// NODE-TO-END branches:
		unsigned short NodeToEnd_Counter;
		double*  NodeToEnd_Length;
		double*  NodeToEnd_MaxWidth;    // Max thickness
		double*  NodeToEnd_MeanWidth;   // Mean thickness
		double*  NodeToEnd_MinWidth;    // Throat size

		// NODE-TO-NODE branches:
		unsigned short NodeToNode_Counter;
		double*  NodeToNode_Length; 
		double*  NodeToNode_MaxWidth;    // Max thickness
		double*  NodeToNode_MeanWidth;   // Mean thickness
		double*  NodeToNode_MinWidth;    // Throat size

		// NODE points:
		unsigned short  Node_Counter;	
		double*  Node_Width; 	

		// Tortuosity indexes:
		/*unsigned short  Tort_Counter;
		double*   Tort_X;		
		double*   Tort_Y;		
		double*   Tort_Z;*/
	};

#endif

	int p3dGVFSkeletonization (	
		unsigned char* in_im,
		unsigned char* skl_im,
		const int dimx,	
		const int dimy,	
		const int dimz,	
		const double mu,
		const double eps,
		const double thresh,
		const double scale,
		int (*wr_log)(const char*, ...)
		);

	int p3dThinningSkeletonization(   
		unsigned char* in_im, 
		unsigned char* out_im, 
		const int dimx,
		const int dimy, 
		const int dimz,
		int (*wr_log)(const char*, ...)
		);

	int p3dLKCSkeletonization(   
		unsigned char* in_im, 
		unsigned char* out_im, 
		const int dimx,
		const int dimy, 
		const int dimz,
		int (*wr_log)(const char*, ...)
		);

	int p3dSimpleSkeletonPruning(   
		unsigned char* in_im,		// IN: Input (binary) skeleton
		unsigned char* out_im,		// OUT: Labeled skeleton
		const int dimx,
		const int dimy, 
		const int dimz,
		const int thresh,
		int (*wr_log)(const char*, ...)
		);

	int p3dIterativeSkeletonPruning(   
		unsigned char* in_im,		// IN: Input (binary) skeleton
		unsigned char* out_im,		// OUT: Labeled skeleton
		const int dimx,
		const int dimy, 
		const int dimz,
		const int thresh,
		int (*wr_log)(const char*, ...)
		);

	int p3dUltimateSkeletonPruning(   
		unsigned char* in_im,		// IN: Input (binary) skeleton
		unsigned char* out_im,		// OUT: Labeled skeleton
		const int dimx,
		const int dimy, 
		const int dimz,
		const int iterative,
		int (*wr_log)(const char*, ...)
		);

	int p3dSkeletonLabeling (   
		unsigned char* in_im,		// IN: Input (binary) skeleton
		unsigned char* out_im,		// OUT: Labeled skeleton
		const int dimx,
		const int dimy, 
		const int dimz,
		int (*wr_log)(const char*, ...)
		);

	int p3dSkeletonAnalysis( 
		unsigned char* vol_im,			    // IN: Input segmented (binary) volume
		unsigned char* skl_im,				// IN: Input (binary) skeleton of the segmented volume
		struct SkeletonStats* out_stats,	// OUT: Skeleton statistics
		unsigned char* nodes_im,			// OUT: Image with maximal balls on skeleton node points
		unsigned char* pores_im,			// OUT: Image with maximal balls on pores
		unsigned char* ends_im,				// OUT: Image with maximal balls on skeleton end points
		unsigned char* throats_im,			// OUT: Image with maximal balls on skeleton throats
		const int dimx,
		const int dimy, 
		const int dimz,
		const double merging_factor,
		const int tortuosity_depth,
		const double voxelsize,			// IN: voxel resolution
		int (*wr_log)(const char*, ...)
		);

	int p3dSkeletonAnalysis_2( 
		unsigned char* vol_im,			    // IN: Input segmented (binary) volume
		unsigned char* skl_im,				// IN: Input (binary) skeleton of the segmented volume
		struct SkeletonStats* out_stats,	// OUT: Skeleton statistics
		unsigned char* nodes_im,			// OUT: Image with maximal balls on skeleton node points
		unsigned char* pores_im,			// OUT: Image with maximal balls on pores
		unsigned char* ends_im,				// OUT: Image with maximal balls on skeleton end points
		unsigned char* throats_im,			// OUT: Image with maximal balls on skeleton throats
		const int dimx,
		const int dimy, 
		const int dimz,
		const double merging_factor,
		const int tortuosity_depth,
		const double voxelsize,			// IN: voxel resolution
		int (*wr_log)(const char*, ...)
		);

	int p3dSkeletonAnalysisFeasibility (
		unsigned char*  in_im,		// IN: binary volume
		unsigned char*	 sk_im,		// IN: skeleton mask
		double* ratio,				// OUT: percentage of total pore space occupied by the skeleton.
		const int dimx,
		const int dimy, 
		const int dimz,	
		int (*wr_log)(const char*, ...)
		);

#ifdef __cplusplus
}
#endif




