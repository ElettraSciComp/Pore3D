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

/********************************************************************* 
 * 
 * uint_list_t type definitions. 
 *
 *********************************************************************/

struct uint_lelem_t {
	unsigned int	ct;
	struct uint_lelem_t	*next;
};

typedef struct uint_lelem_t uint_list_elem_t;
	
typedef uint_list_elem_t* uint_list_t;


/********************************************************************* 
 * 
 * Interface for the queue (FIFO). 
 *
 *********************************************************************/

void uint_list_init (uint_list_t *list);

void uint_list_add(uint_list_t *list, unsigned int ct);

int uint_list_isempty(uint_list_t list);

unsigned int* uint_list_toarray (uint_list_t *list, int numel );

void uint_list_clear (uint_list_t *list);

