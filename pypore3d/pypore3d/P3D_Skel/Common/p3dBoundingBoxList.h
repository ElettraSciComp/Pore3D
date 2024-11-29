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

#include "p3dBoundingBoxT.h"

/********************************************************************* 
 * 
 * uint_list_t type definitions. 
 *
 *********************************************************************/
#ifndef BB_L_DEFINED
	#define BB_L_DEFINED  

	struct bb_lelem_t {
		bb_t elem;
		struct bb_lelem_t	*next;
	};

	typedef struct bb_lelem_t bb_list_elem_t;
		
	typedef bb_list_elem_t* bb_list_t;

#endif
/********************************************************************* 
 * 
 * Interface for the queue (FIFO). 
 *
 *********************************************************************/

void bb_list_init (bb_list_t *list);

int bb_list_add (bb_list_t *list, bb_t item);

int bb_list_isempty(bb_list_t list);

bb_t* bb_list_toarray (bb_list_t *list, int numel );

void bb_list_clear (bb_list_t *list);


