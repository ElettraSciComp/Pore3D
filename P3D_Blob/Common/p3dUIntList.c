#include <stdlib.h>
#include <stdio.h>

#include "p3dUIntList.h"
#include "p3dUtils.h"

void uint_list_init (uint_list_t *list)
{
    *list = NULL;    
}

void uint_list_add (uint_list_t *list, unsigned int elem)
{
	uint_list_elem_t* new_l;

	// Alloc memory for the new item:
	new_l = (uint_list_elem_t*) malloc (sizeof(uint_list_elem_t));

	// Push item into queue:
	new_l->ct = elem;
	new_l->next = *list;

	*list = new_l;
}


int uint_list_isempty (uint_list_t list)
{
    return ( (list == NULL) ? P3D_TRUE : P3D_FALSE );
}

// List is deleted after conversion:
unsigned int* uint_list_toarray (uint_list_t *list, unsigned int numel )
{
	uint_list_elem_t* tmp_l;

	unsigned int* v;
	unsigned int i;	

	v = (unsigned int*) calloc ( numel, sizeof(unsigned int));

	// Convert list to array:	
	for (i = (unsigned int) (numel - 1); i > 0; i = (unsigned int) (i-1) )
	{
		if (list != NULL)
		{
			v[i] = (*list)->ct;
	
			// Perform deletion:
			tmp_l = *list;	
	
			(*list) = (*list)->next;

			free(tmp_l);			
		}
	}

	// Last element:
	if (list != NULL)
	{
		v[0] = (*list)->ct;	
		
		tmp_l = *list;	
	
		(*list) = (*list)->next;
		
		free(tmp_l);
	}

	return v;
}

void uint_list_clear (uint_list_t *list)
{
	uint_list_elem_t* tmp_l;

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