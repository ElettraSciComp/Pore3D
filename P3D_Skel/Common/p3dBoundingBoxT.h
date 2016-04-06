/********************************************************************* 
 * 
 * bb_t type definitions (bounding box).
 *
 *********************************************************************/

#ifndef BB_T_DEFINED

typedef struct {
	int min_x;
	int max_x;
	int min_y;
	int max_y;
	int min_z;
	int max_z;
} bb_t;

#endif

#define BB_T_DEFINED