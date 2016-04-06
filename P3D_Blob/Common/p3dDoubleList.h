/********************************************************************* 
 * 
 * double_list_t type definitions. 
 *
 *********************************************************************/

struct double_lelem_t {
	double	ct;
	struct double_lelem_t	*next;
};

typedef struct double_lelem_t double_list_elem_t;
	
typedef double_list_elem_t* double_list_t;


/********************************************************************* 
 * 
 * Interface for the queue (FIFO). 
 *
 *********************************************************************/

void double_list_init (double_list_t *list);

void double_list_add(double_list_t *list, double item);

int double_list_isempty(double_list_t list);

double* double_list_toarray (double_list_t *list, int numel );

