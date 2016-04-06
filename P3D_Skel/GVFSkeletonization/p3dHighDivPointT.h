/*
	Type definitions:
*/
 
#ifndef HIGHDIV_POINT_T_DEFINED

typedef struct {

	/* Coordinates: */
	double x;
	double y;
	double z;

	/* Divergence value */
	double div;

} highDiv_point_t;

#endif

#define HIGHDIV_POINT_T_DEFINED