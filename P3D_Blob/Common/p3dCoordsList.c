#include <stdlib.h>
#include <stdio.h>

#include "p3dCoordsList.h"
#include "p3dUtils.h"


void coords_list_init (coords_list_t *list)
{
    *list = NULL;    
}

int coords_list_push (coords_list_t *list, coords_t item)
{
	coords_list_elem_t* new_l;
		
	/* Allocate memory for the new item: */
	new_l = (coords_list_elem_t*) malloc (sizeof(coords_list_elem_t));
	if ( new_l == NULL ) return P3D_MEM_ERROR;

	/* Push item into list: */
	new_l->elem = item;
	new_l->next = *list;

	*list = new_l;

	return P3D_SUCCESS;
}

coords_t coords_list_pop (coords_list_t *list)
{
	coords_list_elem_t* tmp_l;
	coords_t tmp_coords;



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


int coords_list_isempty (coords_list_t list)
{
	return ( (list == NULL) ? P3D_TRUE : P3D_FALSE );
}


coords_t* coords_list_toarray (coords_list_t *list, int numel )
{
	coords_list_elem_t* tmp_l;


	coords_t* v;
	int i;

	/* Allocate memory for output array: */
	v = (coords_t*) malloc (numel*sizeof(coords_t));
	
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
			coords_list_pop ( list );
	}

	return v;
}