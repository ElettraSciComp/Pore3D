#include "../Common/p3dFCoordsList.h"
#include "p3dHighDivPointList.h"

int p3dComputeHierarchicalSkeleton (    
	highDiv_point_list_t*   highDiv_point_list,
	fcoords_list_t*         skel_point_list,
    float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,
	const double step,
	const double close_dist
	);