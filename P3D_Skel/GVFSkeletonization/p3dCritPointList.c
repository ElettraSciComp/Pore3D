#include <stdlib.h>
#include <stdio.h>

#include "p3dCritPointList.h"

#include "../Common/p3dUtils.h"

void crit_point_list_init (crit_point_list_t *list)
{
    *list = NULL;    
}

int crit_point_list_push (crit_point_list_t *list, crit_point_t item)
{
	crit_point_list_elem_t* new_l;
	

	// Alloc memory for the new item:
	new_l = (crit_point_list_elem_t*) malloc (sizeof(crit_point_list_elem_t));
	if ( new_l == NULL ) return P3D_MEM_ERROR;

	// Push item into queue:
	new_l->elem = item;
	new_l->next = *list;

	*list = new_l;

	return P3D_SUCCESS;
}

crit_point_t crit_point_list_pop (crit_point_list_t *list)
{
	crit_point_list_elem_t* tmp_l;
	crit_point_t tmp_crit_point;



	// Perform deletion:
	tmp_crit_point = (*list)->elem;

	// Save temporary pointer to the element to delete:
	tmp_l = *list;		

	// The list will point on the next element:
	(*list) = (*list)->next;

	// Delete first element using the previously set 
	// temporary pointer:
	free(tmp_l);

	// Return coordinates:
	return tmp_crit_point;
}


int crit_point_list_isempty (crit_point_list_t list)
{
	return ( (list == NULL) ? P3D_TRUE : P3D_FALSE );
}

