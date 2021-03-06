/***************************************************************************/
/* (C) 2016 Elettra - Sincrotrone Trieste S.C.p.A.. All rights reserved.   */
/*                                                                         */
/*                                                                         */
/* This file is part of Pore3D, a software library for quantitative        */
/* analysis of 3D (volume) images.                                         */
/*                                                                         */
/* Pore3D is free software: you can redistribute it and/or modify it       */
/* under the terms of the GNU General Public License as published by the   */
/* Free Software Foundation, either version 3 of the License, or (at your  */
/* option) any later version.                                              */
/*                                                                         */
/* Pore3D is distributed in the hope that it will be useful, but WITHOUT   */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   */
/* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    */
/* for more details.                                                       */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with Pore3D. If not, see <http://www.gnu.org/licenses/>.          */
/*                                                                         */
/***************************************************************************/

//
// Author: Francesco Brun
// Last modified: Sept, 28th 2016
//

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