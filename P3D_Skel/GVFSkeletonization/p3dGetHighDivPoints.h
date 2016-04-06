#include "p3dHighDivPointList.h"

int p3dGetHighDivPoints (
	unsigned char* in_im,
	float* gvf_x,
	float* gvf_y,
	float* gvf_z,
	const int dimx,	
	const int dimy,	
	const int dimz,	
	const double thresh,
	highDiv_point_list_t* highDiv_point_list	
	);