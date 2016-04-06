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
		
	typedef struct {
		coords_queue_elem_t *tail;
		coords_queue_elem_t *head;
	} coords_queue_t;

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

