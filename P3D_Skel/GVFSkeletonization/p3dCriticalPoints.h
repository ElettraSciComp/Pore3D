#include "p3dCritPointList.h"

int p3dGetCriticalPoints (
	unsigned char* in_im,
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,	
	crit_point_list_t* crit_point_list
	);

int p3dClassifyCriticalPoints (
    crit_point_list_t crit_point_list,
    float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz	
	);