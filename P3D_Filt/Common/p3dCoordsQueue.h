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

#include "p3dCoordsT.h"

/********************************************************************* 
 * 
 * coords_queue_t type definitions. Do NOT modify. coords_queue_t implements a FIFO 
 * policy: elements are pushed into the tail and popped from the head.
 *
 *********************************************************************/

#ifndef COORDS_Q_DEFINED
	#define COORDS_Q_DEFINED  

	struct coords_qelem_t {
		coords_t	item;
		struct coords_qelem_t	*next;
	};

	typedef struct coords_qelem_t coords_queue_elem_t;
		
	struct coords_q_t{
		coords_queue_elem_t *tail;
		coords_queue_elem_t *head;
	};
        
        typedef struct coords_q_t coords_queue_t;

#endif
/********************************************************************* 
 * 
 * Interface for the queue (FIFO). 
 *
 *********************************************************************/

void coords_queue_init (coords_queue_t *queue);

void coords_queue_push(coords_queue_t *queue, coords_t elem);

coords_t coords_queue_pop(coords_queue_t *queue);

int coords_queue_isempty(coords_queue_t queue);

