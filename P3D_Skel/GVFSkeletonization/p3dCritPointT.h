/*
	Type definitions:
*/
 
#ifndef CRIT_POINT_T_DEFINED

typedef struct {

	/* Coordinates: */
	double x;
	double y;
	double z;

	/* Type { CPT_SADDLE | CPT_ATTRACTING_NODE | CPT_REPELLING_NODE } */
	int type;

	/* Eigenvalues: */
	double eval0;
	double eval1;
	double eval2;

	/* Eigenvectors: */
	double evect0_x;
	double evect0_y;
	double evect0_z;

	double evect1_x;
	double evect1_y;
	double evect1_z;

	double evect2_x;
	double evect2_y;
	double evect2_z;

} crit_point_t;

#endif

#define CRIT_POINT_T_DEFINED