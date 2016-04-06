#include <stdlib.h>
#include <stdio.h>

#include "p3dBoundingBoxList.h"
#include "p3dUtils.h"

void bb_list_init (bb_list_t *list)
{
    *list = NULL;    
}

int bb_list_add (bb_list_t *list, bb_t item)
{
	bb_list_elem_t* new_l;

	// Alloc memory for the new item:
	new_l = (bb_list_elem_t*) malloc (sizeof(bb_list_elem_t));

	// Push item into queue:
	new_l->elem = item;
	new_l->next = *list;

	*list = new_l;

	return P3D_TRUE;
}


int bb_list_isempty (bb_list_t list)
{
    return ( (list == NULL) ? P3D_TRUE : P3D_FALSE );
}

// List is deleted after conversion:
bb_t* bb_list_toarray (bb_list_t *list, unsigned int numel )
{
	bb_list_elem_t* tmp_l;

	bb_t* v;
	unsigned int i;
	FILE *fp;

	v = (bb_t*) calloc ( (unsigned int) numel, sizeof(bb_t));	

	// Convert list to array:	
	for (i = (unsigned int) (numel - 1); i > 0; i = (unsigned int) (i - 1) )
	{
		if (list != NULL)
		{
			v[i] = (*list)->elem;

			// Perform deletion:
			tmp_l = *list;	
	
			(*list) = (*list)->next;

			free(tmp_l);			
		}
	}

	// Last element:
	if (list != NULL)
	{
		v[0] = (*list)->elem;	
		
		tmp_l = *list;	
	
		(*list) = (*list)->next;
		
		free(tmp_l);
	}

	return v;
}

void bb_list_clear (bb_list_t *list)
{
	bb_list_elem_t* tmp_l;

	// Scan whole list:
	while( *list != NULL )
	{
		// Perform deletion:
		tmp_l = *list;	
	
		(*list) = (*list)->next;

		free(tmp_l);	
	}

	// Final assignment for safety:
	*list = NULL;
}