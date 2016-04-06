#include "p3dBoundingBoxT.h"

/********************************************************************* 
 * 
 * uint_list_t type definitions. 
 *
 *********************************************************************/
#ifndef BB_L_DEFINED
	#define BB_L_DEFINED  

	struct bb_lelem_t {
		bb_t elem;
		struct bb_lelem_t	*next;
	};

	typedef struct bb_lelem_t bb_list_elem_t;
		
	typedef bb_list_elem_t* bb_list_t;

#endif
/********************************************************************* 
 * 
 * Interface for the queue (FIFO). 
 *
 *********************************************************************/

void bb_list_init (bb_list_t *list);

int bb_list_add (bb_list_t *list, bb_t item);

int bb_list_isempty(bb_list_t list);

bb_t* bb_list_toarray (bb_list_t *list, unsigned int numel );

void bb_list_clear (bb_list_t *list);


