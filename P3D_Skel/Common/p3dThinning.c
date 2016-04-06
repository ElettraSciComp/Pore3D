#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "p3dUtils.h"

#define BORDER 100	// Custom label to mark a border voxel
#define SIMPLE 200  // Custom label to mark a simple voxel


#define UP_SOUTH    0
#define NORTH_EAST  1
#define WEST_DOWN   2 

#define EAST_SOUTH  3
#define UP_WEST     4
#define NORTH_DOWN  5

#define SOUTH_WEST  6
#define UP_NORTH    7
#define EAST_DOWN   8

#define NORTH_WEST  9
#define UP_EAST    10
#define SOUTH_DOWN 11

#define UP         12
#define DOWN       13 
#define EAST       14 
#define WEST       15 
#define NORTH      16
#define SOUTH      17


int  checkTemplate ( int n[27] ) 
{
	// T1
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) ) 
		&&
		((n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,2,3,3)] == P3D_FALSE) && (n[I(1,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,2,3,3)] == P3D_FALSE) )
		&&
		((n[I(0,0,0,3,3)] == P3D_TRUE) || (n[I(1,0,0,3,3)] == P3D_TRUE) || 
		 (n[I(2,0,0,3,3)] == P3D_TRUE) || (n[I(0,0,1,3,3)] == P3D_TRUE) || 
		 (n[I(2,0,1,3,3)] == P3D_TRUE) || (n[I(0,0,2,3,3)] == P3D_TRUE) || 
		 (n[I(1,0,2,3,3)] == P3D_TRUE) || (n[I(2,0,2,3,3)] == P3D_TRUE) ||
		 (n[I(0,1,0,3,3)] == P3D_TRUE) || (n[I(1,1,0,3,3)] == P3D_TRUE) || 
		 (n[I(2,1,0,3,3)] == P3D_TRUE) || (n[I(0,1,1,3,3)] == P3D_TRUE) || 
		 (n[I(2,1,1,3,3)] == P3D_TRUE) || (n[I(0,1,2,3,3)] == P3D_TRUE) || 
		 (n[I(1,1,2,3,3)] == P3D_TRUE) || (n[I(2,1,2,3,3)] == P3D_TRUE) )
		)
	{
		return P3D_TRUE;
	}

	// T2
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,1,2,3,3)] == P3D_TRUE)) 
		&&
		((n[I(0,0,0,3,3)] == P3D_FALSE) && (n[I(1,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,0,0,3,3)] == P3D_FALSE) && (n[I(0,1,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(2,1,0,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE))
		&&
		((n[I(0,0,1,3,3)] == P3D_TRUE) || (n[I(1,0,1,3,3)] == P3D_TRUE) || 
		 (n[I(2,0,1,3,3)] == P3D_TRUE) || (n[I(0,1,1,3,3)] == P3D_TRUE) || 
		 (n[I(2,1,1,3,3)] == P3D_TRUE) || (n[I(0,2,1,3,3)] == P3D_TRUE) || 
		 (n[I(1,2,1,3,3)] == P3D_TRUE) || (n[I(2,2,1,3,3)] == P3D_TRUE) ||
		 (n[I(0,0,2,3,3)] == P3D_TRUE) || (n[I(1,0,2,3,3)] == P3D_TRUE) || 
		 (n[I(2,0,2,3,3)] == P3D_TRUE) || (n[I(0,1,2,3,3)] == P3D_TRUE) || 
		 (n[I(2,1,2,3,3)] == P3D_TRUE) || (n[I(0,2,2,3,3)] == P3D_TRUE) || 
		 (n[I(1,2,2,3,3)] == P3D_TRUE) || (n[I(2,2,2,3,3)] == P3D_TRUE))
		)
	{
		return P3D_TRUE;
	}

	// T3
	if(( (n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,2,3,3)] == P3D_TRUE) ) 
		&&
		((n[I(0,0,0,3,3)] == P3D_FALSE) && (n[I(1,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,0,0,3,3)] == P3D_FALSE) && (n[I(0,1,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(2,1,0,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,2,3,3)] == P3D_FALSE) && (n[I(1,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,2,3,3)] == P3D_FALSE) )
		&&
		((n[I(0,0,1,3,3)] == P3D_TRUE) || (n[I(2,0,1,3,3)] == P3D_TRUE) || 
		 (n[I(0,0,2,3,3)] == P3D_TRUE) || (n[I(2,0,2,3,3)] == P3D_TRUE) ||
		 (n[I(0,1,1,3,3)] == P3D_TRUE) || (n[I(2,1,1,3,3)] == P3D_TRUE) || 
		 (n[I(0,1,2,3,3)] == P3D_TRUE) || (n[I(2,1,2,3,3)] == P3D_TRUE) )
		 )
	{
		return P3D_TRUE;
	}

	// T4
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(1,1,2,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(0,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,0,3,3)] == P3D_FALSE) && (n[I(2,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) )
		&&
		/*!(n[0][1][0] && n[0][2][1])
		&&
		!(n[2][1][0] && n[2][2][1])*/
		( (n[I(0,1,0,3,3)] == P3D_FALSE) || (n[I(0,2,1,3,3)] == P3D_FALSE) ) 
		&&
		( (n[I(2,1,0,3,3)] == P3D_FALSE) || (n[I(2,2,1,3,3)] == P3D_FALSE) ) 
		)
	{
		return P3D_TRUE;
	}

	// T5
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(1,1,2,3,3)] == P3D_TRUE) && (n[I(2,2,0,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(0,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,1,3,3)] == P3D_FALSE) )
		&&
		( (n[I(0,1,0,3,3)] == P3D_FALSE) || (n[I(0,2,1,3,3)] == P3D_FALSE) ) 
		&&
		(((n[I(2,1,0,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_TRUE)) 
		|| 
		((n[I(2,1,0,3,3)] == P3D_TRUE) && (n[I(2,2,1,3,3)] == P3D_FALSE)))
		)
	{
		return P3D_TRUE;
	}

	// T6
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(1,1,2,3,3)] == P3D_TRUE) && (n[I(0,2,0,3,3)] == P3D_TRUE)) 
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,1,3,3)] == P3D_FALSE) )
		&& 
		(((n[I(0,1,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_TRUE)) || 
		((n[I(0,1,0,3,3)] == P3D_TRUE) && (n[I(0,2,1,3,3)] == P3D_FALSE)))
		&&
		//!(n[2][1][0] && n[2][2][1])
		((n[I(2,1,0,3,3)] == P3D_FALSE) || (n[I(2,2,1,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T7
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(2,1,1,3,3)] == P3D_TRUE) && (n[I(1,1,2,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(0,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,1,3,3)] == P3D_FALSE))
		&&
		//!(n[0][1][0] && n[0][2][1])
		((n[I(0,1,0,3,3)] == P3D_FALSE) || (n[I(0,2,1,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T8
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(0,1,1,3,3)] == P3D_TRUE) && (n[I(1,1,2,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,1,3,3)] == P3D_FALSE))
		&&
		((n[I(2,1,0,3,3)] == P3D_FALSE) || (n[I(2,2,1,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T9
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(2,1,1,3,3)] == P3D_TRUE) && (n[I(1,1,2,3,3)] == P3D_TRUE) && 
		(n[I(0,2,0,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		(n[I(1,2,1,3,3)] == P3D_FALSE))
		&&
		(((n[I(0,1,0,3,3)] == P3D_TRUE) && (n[I(0,2,1,3,3)] == P3D_FALSE)) || 
		((n[I(0,1,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_TRUE)))
		)
	{
		return P3D_TRUE;
	}

	// T10
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,1,3,3)] == P3D_TRUE) && 
		(n[I(0,1,1,3,3)] == P3D_TRUE) && (n[I(1,1,2,3,3)] == P3D_TRUE) && 
		(n[I(2,2,0,3,3)] == P3D_TRUE))
		&&
		((n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) 
		&& (n[I(1,2,1,3,3)] == P3D_FALSE))
		&&
		(((n[I(2,1,0,3,3)] == P3D_TRUE) && (n[I(2,2,1,3,3)] == P3D_FALSE)) || 
		((n[I(2,1,0,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_TRUE)))
		)
	{
		return P3D_TRUE;
	}

	// T11
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(2,0,1,3,3)] == P3D_TRUE) && 
		(n[I(1,0,2,3,3)] == P3D_TRUE))
		&&
		((n[I(0,0,0,3,3)] == P3D_FALSE) && (n[I(1,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(0,1,0,3,3)] == P3D_FALSE) && (n[I(1,1,0,3,3)] == P3D_FALSE) && 
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,2,3,3)] == P3D_FALSE) && (n[I(1,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,2,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T12
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(0,0,1,3,3)] == P3D_TRUE) && 
		(n[I(1,0,2,3,3)] == P3D_TRUE))
		&&
		((n[I(1,0,0,3,3)] == P3D_FALSE) && (n[I(2,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(2,1,0,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) && (n[I(2,2,1,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,2,3,3)] == P3D_FALSE) && (n[I(1,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,2,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T13
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,2,3,3)] == P3D_TRUE) && 
		(n[I(2,1,2,3,3)] == P3D_TRUE))
		&&
		((n[I(0,0,0,3,3)] == P3D_FALSE) && (n[I(1,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,0,0,3,3)] == P3D_FALSE) && (n[I(0,1,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(2,1,0,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(0,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,1,3,3)] == P3D_FALSE) && (n[I(0,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(1,2,2,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	// T14
	if(((n[I(1,1,1,3,3)] == P3D_TRUE) && (n[I(1,0,2,3,3)] == P3D_TRUE) && 
		(n[I(0,1,2,3,3)] == P3D_TRUE))
		&&
		((n[I(0,0,0,3,3)] == P3D_FALSE) && (n[I(1,0,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,0,0,3,3)] == P3D_FALSE) && (n[I(0,1,0,3,3)] == P3D_FALSE) && 
		 (n[I(1,1,0,3,3)] == P3D_FALSE) && (n[I(2,1,0,3,3)] == P3D_FALSE) &&
		 (n[I(0,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,0,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,0,3,3)] == P3D_FALSE) && (n[I(1,2,1,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,1,3,3)] == P3D_FALSE) && (n[I(1,2,2,3,3)] == P3D_FALSE) && 
		 (n[I(2,2,2,3,3)] == P3D_FALSE))
		)
	{
		return P3D_TRUE;
	}

	return P3D_FALSE;
}



void transformNeigh ( int n[27], int dir, int USn[27] ) 
{
	int i, j, k;
	int tmp[27];

	switch(dir) {
  case 0:  //UP_SOUTH = 0, 
	  // just copy
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,j,k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

	  // The following cases are clearly rotations only
  case 3: //   EAST_SOUTH,

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,j,k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 6: //  SOUTH_WEST,

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(j,2-i,k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 10: //   UP_EAST, 

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(k,j,2-i,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;
  case 4: //  UP_WEST, 

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(2-k,j,i,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;  

  case 11: //    SOUTH_DOWN  
	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,2-j,k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 7: //  UP_NORTH, 

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,j,2-k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 5: //   NORTH_DOWN, 

	  // OR - two reflections
	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  tmp[I(i,j,2-k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  // 2
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,2-j,k,3,3)] = tmp[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

	case 8: //   EAST_DOWN, 

	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  tmp[I(k,j,2-i,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  // 2
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,2-j,k,3,3)] = tmp[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 2: // WEST_DOWN, 
	  // 1
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  tmp[I(2-k,j,i,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  // 2
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(i,2-j,k,3,3)] = tmp[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 1:  //NORTH_EAST,

	  // 1 - reflection
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  tmp[I(i,j,2-k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  // 2 - rotation
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(2-j,i,k,3,3)] = tmp[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

  case 9: //   NORTH_WEST,

	  // 1 - reflection
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  tmp[I(i,j,2-k,3,3)] = n[I(i,j,k,3,3)];
			  }
		  }
	  }
	  // 2 - rotation
	  for(k=0; k < 3; k++) {
		  for(j=0; j < 3; j++) {
			  for(i=0; i < 3; i++) {
				  USn[I(j,2-i,k,3,3)] = tmp[I(i,j,k,3,3)];
			  }
		  }
	  }
	  break;

	}
}


void markBoundary ( unsigned char *vol, int L, int M, int N, int dir )
{
	int slsz = L*M;
	int idx;
	int i, j, k;


	// neighbor index in 18 directions (only last 6 used)
	int nb[18] = 
	{
		+L - slsz,  // UP_SOUTH,   0
		+slsz + 1,	// NORTH_EAST, 1
		-1 - L,     // WEST_DOWN,  2 

		+1 -slsz,	// EAST_SOUTH, 3
		+L - 1,     // UP_WEST,    4
		+slsz - L,  // NORTH_DOWN, 5 

		-slsz - 1,	// SOUTH_WEST, 6
		+L + slsz,  // UP_NORTH,   7
		+1 - L,     // EAST_DOWN,  8

		+slsz - 1,	// NORTH_WEST, 9
		+L + 1,     // UP_EAST,   10
		-slsz - L,  // SOUTH_DOWN,11

		+L,			// UP,        12
		-L,			// DOWN,      13
		+1,			// EAST,      14
		-1,			// WEST,      15
		+slsz,		// NORTH,     16
		-slsz		// SOUTH,     17
	};


	for( k = 1; k < (N-1); k++ ) 
		for( j = 1; j < (M-1); j++ ) 
			for( i = 1; i < (L-1); i++ ) 
			{
				idx = k*slsz + j*L + i;

				if( (vol[idx] == OBJECT) && ( vol[idx + nb[dir]] == 0) ) 
				{
					vol[idx] = BORDER;
				}
			}	
}









void bufferizeNeigh ( unsigned char *vol, int idx, int nb[27], const int volNeigh[27] )
{
	int nidx;
	int i, j, k, ii;

	ii = 0;
	for( k = 0; k < 3; k++ ) 
		for( j = 0; j < 3; j++ ) 
			for( i = 0; i < 3; i++ ) 
			{
				nidx = idx + volNeigh[ii++];

				nb[I(i,j,k,3,3)] = ( vol[nidx] != BACKGROUND ) ? P3D_TRUE : P3D_FALSE;
			}
}


int p3dThinning (   
	unsigned char* in_rev, 
	const int dimx,
	const int dimy, 
	const int dimz	
	) 
{	
	int nrDel;

	int dir;
	int nb[27];
	int USn[27];
	int idx;
	int i, j, k;

	// Initialize neighbors array:
	const int volNeighbors[27] = 
	{		
		// lower plane
		(-(dimx*dimy) -dimx -1),
		(-(dimx*dimy) -dimx +0),
		(-(dimx*dimy) -dimx +1),
		(-(dimx*dimy) +0 -1),
		(-(dimx*dimy) +0 +0),
		(-(dimx*dimy) +0 +1),
		(-(dimx*dimy) +dimx -1),
		(-(dimx*dimy) +dimx +0),
		(-(dimx*dimy) +dimx +1),
		// same plane
		(+0 -dimx -1),
		(+0 -dimx +0),
		(+0 -dimx +1),
		(+0 +0 -1),
		(+0 +0 +0),
		(+0 +0 +1),
		(+0 +dimx -1),
		(+0 +dimx +0),
		(+0 +dimx +1),
		// upper plane
		(+(dimx*dimy) -dimx -1),
		(+(dimx*dimy) -dimx +0),
		(+(dimx*dimy) -dimx +1),
		(+(dimx*dimy) +0 -1),
		(+(dimx*dimy) +0 +0),
		(+(dimx*dimy) +0 +1),
		(+(dimx*dimy) +dimx -1),
		(+(dimx*dimy) +dimx +0),
		(+(dimx*dimy) +dimx +1),
	};

	nrDel = 1;


	// Loop until no more voxels can be deleted:
	while( nrDel > 0 ) 
	{
		nrDel = 0;

		// Perform the 12 sub-iteration thinning (in parallel ???)
		//#pragma omp parallel for private(i, j, k, idx, nb, USn)  reduction (+ : nrDel )
		for(dir = 0; dir < 12; dir++) 
		{
			switch(dir) 
			{
				  case UP_SOUTH:
					  // UP
					  markBoundary(in_rev, dimx, dimy, dimz, UP);
					  // SOUTH
					  markBoundary(in_rev, dimx, dimy, dimz, SOUTH);
					  break;
				  case NORTH_EAST:
					  // NOTH
					  markBoundary(in_rev, dimx, dimy, dimz, NORTH);
					  // EAST
					  markBoundary(in_rev, dimx, dimy, dimz, EAST);
					  break;
				  case WEST_DOWN:
					  // WEST
					  markBoundary(in_rev, dimx, dimy, dimz, WEST);
					  // DOWN
					  markBoundary(in_rev, dimx, dimy, dimz, DOWN);
					  break;
				  case EAST_SOUTH:
					  // EAST
					  markBoundary(in_rev, dimx, dimy, dimz, EAST);
					  // SOUTH
					  markBoundary(in_rev, dimx, dimy, dimz, SOUTH);
					  break;
				  case UP_WEST:
					  // UP
					  markBoundary(in_rev, dimx, dimy, dimz, UP);
					  // WEST
					  markBoundary(in_rev, dimx, dimy, dimz, WEST);
					  break;
				  case NORTH_DOWN:
					  // NORTH
					  markBoundary(in_rev, dimx, dimy, dimz, NORTH);
					  // DOWN
					  markBoundary(in_rev, dimx, dimy, dimz, DOWN);
					  break;
				  case SOUTH_WEST:
					  // SOUTH
					  markBoundary(in_rev, dimx, dimy, dimz, SOUTH);
					  // WEST
					  markBoundary(in_rev, dimx, dimy, dimz, WEST);
					  break;
				  case UP_NORTH:
					  // UP
					  markBoundary(in_rev, dimx, dimy, dimz, UP);
					  // NORTH
					  markBoundary(in_rev, dimx, dimy, dimz, NORTH);
					  break;
				  case EAST_DOWN:
					  // EAST
					  markBoundary(in_rev, dimx, dimy, dimz, EAST);
					  // DOWN
					  markBoundary(in_rev, dimx, dimy, dimz, DOWN);	
					  break;
				  case NORTH_WEST:
					  // NORTH
					  markBoundary(in_rev, dimx, dimy, dimz, NORTH);
					  // WEST
					  markBoundary(in_rev, dimx, dimy, dimz, WEST);
					  break;
				  case UP_EAST:
					  // UP
					  markBoundary(in_rev, dimx, dimy, dimz, UP);
					  // EAST
					  markBoundary(in_rev, dimx, dimy, dimz, EAST);
					  break;
				  case SOUTH_DOWN:
					  // SOUTH
					  markBoundary(in_rev, dimx, dimy, dimz, SOUTH);
					  // DOWN
					  markBoundary(in_rev, dimx, dimy, dimz, DOWN);	
					  break;
			}

			
			// Check each boundary point and remove it if it macthes a template:
			#pragma omp parallel for private(i, j, idx, nb, USn)  reduction (+ : nrDel )
			for( k = 1; k < (dimz-1); k++ ) 			
				for( j = 1; j < (dimy-1); j++ ) 				
					for( i = 1; i < (dimx-1); i++ ) 
					{
						idx = I(i,j,k,dimx,dimy);

						if( in_rev[ idx ] == BORDER ) 
						{
							// Copy neighborhood into buffer:							
							bufferizeNeigh(in_rev, idx, nb, volNeighbors);

							transformNeigh(nb, dir, USn);								

							if( checkTemplate(USn) == P3D_TRUE ) 
							{		  
								// Mark as SIMPLE a point that can be removed:
								in_rev[ idx ] = SIMPLE;
								nrDel++;								
							}
						}
					}
				

			// Reset all object voxels to OBJECT and delete simple points:
			#pragma omp parallel for
			for( idx = 0; idx < (dimx*dimy*dimz); idx++ ) 
			{
				if(in_rev[idx] == SIMPLE) in_rev[idx] = BACKGROUND;
				if(in_rev[idx] != BACKGROUND) in_rev[idx] = OBJECT;
			}
		}
	}	

	// Return OK:
	return P3D_SUCCESS;
}

