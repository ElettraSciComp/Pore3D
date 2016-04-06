#include <stdlib.h>

#include "p3dCoordsQueue.h"
#include "p3dUtils.h"

void coords_queue_init (coords_queue_t *queue)
{
    queue->head = NULL;
    queue->tail = NULL;
}

void coords_queue_push (coords_queue_t *queue, coords_t elem)
{
	coords_queue_elem_t *new_q;
	
	// Alloc memory for the new item:
	new_q = (coords_queue_elem_t *) malloc (sizeof(coords_queue_elem_t));

	// Push item into queue:
	new_q->item = elem;
	new_q->next = NULL;

	// Handle first element pushed:
    if (queue->tail != NULL)
	{
		queue->tail->next = new_q;
		queue->tail = queue->tail->next;
	}
	else
	{
		queue->tail = new_q;
        queue->head = queue->tail;        
    }
}

coords_t coords_queue_pop (coords_queue_t *queue)
{    
	coords_t elem;
	coords_queue_elem_t *temp;

	// Pop item from queue:
	elem = queue->head->item;

	// Free memory:
	temp = queue->head;
	queue->head = queue->head->next;
	free(temp);		

	// Handle last element popped:
	if (queue->head == NULL) 
	{
		queue->tail = NULL;
	}

	// Return the item previously popped:
	return elem;
}

int coords_queue_isempty (coords_queue_t queue)
{
	return ( ((queue.head == NULL) && (queue.tail == NULL)) ? P3D_TRUE : P3D_FALSE );
}