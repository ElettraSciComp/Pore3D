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

#include <stdlib.h>
#include <stdio.h>

#include "p3dFCoordsList.h"
#include "p3dUtils.h"

void fcoords_list_init (fcoords_list_t *list)
{
    *list = NULL;    
}

int fcoords_list_push (fcoords_list_t *list, fcoords_t item)
{
	fcoords_list_elem_t* new_l;
		
	/* Allocate memory for the new item: */
	new_l = (fcoords_list_elem_t*) malloc (sizeof(fcoords_list_elem_t));
	if ( new_l == NULL ) return P3D_ERROR;

	/* Push item into list: */
	new_l->elem = item;
	new_l->next = *list;

	*list = new_l;

	return P3D_SUCCESS;
}

fcoords_t fcoords_list_pop (fcoords_list_t *list)
{
	fcoords_list_elem_t* tmp_l;
	fcoords_t tmp_coords;



	/* Perform deletion: */
	tmp_coords = (*list)->elem;

	/* Save temporary pointer to the element to delete: */
	tmp_l = *list;		

	/* The list will point on the next element: */
	(*list) = (*list)->next;

	/* Delete first element using the previously set 
	   temporary pointer: */
	free(tmp_l);

	/* Return coordinates: */
	return tmp_coords;
}


int fcoords_list_isempty (fcoords_list_t list)
{
	return ( (list == NULL) ? P3D_TRUE : P3D_FALSE );
}


fcoords_t* fcoords_list_toarray (fcoords_list_t *list, int numel )
{
	fcoords_list_elem_t* tmp_l;


	fcoords_t* v;
	int i;

	/* Allocate memory for output array: */
	v = (fcoords_t*) malloc (numel*sizeof(fcoords_t));
	
	if ( v != NULL ) 
	{
		/* Convert list to array:	*/
		for (i = (numel - 1); i >= 0; i-- )
		{
			v[i] = (*list)->elem;

			/* Perform deletion: */
			tmp_l = *list;	
		
			(*list) = (*list)->next;

			free(tmp_l);	
		}
	}
	else
	{
		/* Delete list in any case: */
		while ( list != NULL )
			fcoords_list_pop ( list );
	}

	return v;
}