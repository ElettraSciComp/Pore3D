#include <stdlib.h>
#include <stdio.h>

#include "p3dDoubleList.h"
#include "p3dUtils.h"

void double_list_init (double_list_t *list)
{
    *list = NULL;    
}

void double_list_add (double_list_t *list, double elem)
{
	double_list_elem_t* new_l;

	// Alloc memory for the new item:
	new_l = (double_list_elem_t*) malloc (sizeof(double_list_elem_t));

	// Push item into queue:
	new_l->ct = elem;
	new_l->next = *list;

	*list = new_l;
}


int double_list_isempty (double_list_t list)
{
    return ( ( list == NULL )? P3D_TRUE : P3D_FALSE );
}

// List is deleted after conversion:
double* double_list_toarray (double_list_t *list, int numel )
{
	double_list_elem_t* tmp_l;


	double* v;
	int i;

	v = (double*) malloc (numel*sizeof(double));

	// Convert list to array:	
	for (i = (numel - 1); i >= 0; i-- )
	{
		v[i] = (*list)->ct;

		// Perform deletion:
		tmp_l = *list;	
	
		(*list) = (*list)->next;

		free(tmp_l);	
	}

	return v;
}