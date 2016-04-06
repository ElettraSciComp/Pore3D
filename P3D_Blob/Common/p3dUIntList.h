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

unsigned int* uint_list_toarray (uint_list_t *list, unsigned int numel );

void uint_list_clear (uint_list_t *list);

