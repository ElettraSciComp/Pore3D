#include <string.h>
#include <omp.h>

#include "p3dUtils.h"



int findNeighbor ( 
	unsigned char* im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,
	const int j,
	const int k,
	coords_t* coords
	)
{
	int a,b,c;
	int ct = 0;

	// 6-connection:
	c = k - 1; b = j; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	// Do other 12 tests for 18-connection

	// On k-1 plane:
	c = k - 1; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}


	// On k plane:
	c = k; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	// Do other 8 tests for 26-connectivity
	
	// On k-1 plane:
	c = k - 1; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k - 1; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	c = k + 1; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;
	if (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) { 
		coords->x = a; 
		coords->y = b;
		coords->z = c;
	}

	// Return number of voxels in the neighborhood:
	return ct;
}

// Return the number of voxels in the neighborhood:
int countNeighbors ( 
	unsigned char* im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,
	const int j,
	const int k
	)
{
	int a,b,c;
	int ct = 0;

	// 6-connection:
	c = k - 1; b = j; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	// Do other 12 tests for 18-connection

	// On k-1 plane:
	c = k - 1; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j + 1; a = i;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;


	// On k plane:
	c = k; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	// Do other 8 tests for 26-connectivity
	
	// On k-1 plane:
	c = k - 1; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k - 1; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;


	// On k+1 plane:
	c = k + 1; b = j - 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j - 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j + 1; a = i - 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	c = k + 1; b = j + 1; a = i + 1;
	ct += (im[ I(a,b,c,dimx,dimy) ] != BACKGROUND) ? 1 : 0;

	// Return number of voxels in the neighborhood:
	return ct;
}

/** 
 * Octree_labeling [Lee94]
 * This is a recursive method that calulates the number of connected
 * components in the 3D neighbourhood after the center pixel would
 * have been removed.
 */
void Octree_labeling(int octant, int label, int *cube)
{
  // check if there are points in the octant with value 1
  if( octant==1 )
  {
  	// set points in this octant to current label
  	// and recurseive labeling of adjacent octants
    if( cube[0] == 1 )
      cube[0] = label;
    if( cube[1] == 1 )
    {
      cube[1] = label;        
      Octree_labeling( 2, label, cube);
    }
    if( cube[3] == 1 )
    {
      cube[3] = label;        
      Octree_labeling( 3, label, cube);
    }
    if( cube[4] == 1 )
    {
      cube[4] = label;        
      Octree_labeling( 2, label, cube);
      Octree_labeling( 3, label, cube);
      Octree_labeling( 4, label, cube);
    }
    if( cube[9] == 1 )
    {
      cube[9] = label;        
      Octree_labeling( 5, label, cube);
    }
    if( cube[10] == 1 )
    {
      cube[10] = label;        
      Octree_labeling( 2, label, cube);
      Octree_labeling( 5, label, cube);
      Octree_labeling( 6, label, cube);
    }
    if( cube[12] == 1 )
    {
      cube[12] = label;        
      Octree_labeling( 3, label, cube);
      Octree_labeling( 5, label, cube);
      Octree_labeling( 7, label, cube);
    }
  }
  if( octant==2 )
  {
    if( cube[1] == 1 )
    {
      cube[1] = label;
      Octree_labeling( 1, label, cube);
    }
    if( cube[4] == 1 )
    {
      cube[4] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 3, label, cube);
      Octree_labeling( 4, label, cube);
    }
    if( cube[10] == 1 )
    {
      cube[10] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 5, label, cube);
      Octree_labeling( 6, label, cube);
    }
    if( cube[2] == 1 )
      cube[2] = label;        
    if( cube[5] == 1 )
    {
      cube[5] = label;        
      Octree_labeling( 4, label, cube);
    }
    if( cube[11] == 1 )
    {
      cube[11] = label;        
      Octree_labeling( 6, label, cube);
    }
    if( cube[13] == 1 )
    {
      cube[13] = label;        
      Octree_labeling( 4, label, cube);
      Octree_labeling( 6, label, cube);
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==3 )
  {
    if( cube[3] == 1 )
    {
      cube[3] = label;        
      Octree_labeling( 1, label, cube);
    }
    if( cube[4] == 1 )
    {
      cube[4] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 2, label, cube);
      Octree_labeling( 4, label, cube);
    }
    if( cube[12] == 1 )
    {
      cube[12] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 5, label, cube);
      Octree_labeling( 7, label, cube);
    }
    if( cube[6] == 1 )
      cube[6] = label;        
    if( cube[7] == 1 )
    {
      cube[7] = label;        
      Octree_labeling( 4, label, cube);
    }
    if( cube[14] == 1 )
    {
      cube[14] = label;        
      Octree_labeling( 7, label, cube);
    }
    if( cube[15] == 1 )
    {
      cube[15] = label;        
      Octree_labeling( 4, label, cube);
      Octree_labeling( 7, label, cube);
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==4 )
  {
  	if( cube[4] == 1 )
    {
      cube[4] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 2, label, cube);
      Octree_labeling( 3, label, cube);
    }
  	if( cube[5] == 1 )
    {
      cube[5] = label;        
      Octree_labeling( 2, label, cube);
    }
    if( cube[13] == 1 )
    {
      cube[13] = label;        
      Octree_labeling( 2, label, cube);
      Octree_labeling( 6, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[7] == 1 )
    {
      cube[7] = label;        
      Octree_labeling( 3, label, cube);
    }
    if( cube[15] == 1 )
    {
      cube[15] = label;        
      Octree_labeling( 3, label, cube);
      Octree_labeling( 7, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[8] == 1 )
      cube[8] = label;        
    if( cube[16] == 1 )
    {
      cube[16] = label;        
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==5 )
  {
  	if( cube[9] == 1 )
    {
      cube[9] = label;        
      Octree_labeling( 1, label, cube);
    }
    if( cube[10] == 1 )
    {
      cube[10] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 2, label, cube);
      Octree_labeling( 6, label, cube);
    }
    if( cube[12] == 1 )
    {
      cube[12] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 3, label, cube);
      Octree_labeling( 7, label, cube);
    }
    if( cube[17] == 1 )
      cube[17] = label;        
    if( cube[18] == 1 )
    {
      cube[18] = label;        
      Octree_labeling( 6, label, cube);
    }
    if( cube[20] == 1 )
    {
      cube[20] = label;        
      Octree_labeling( 7, label, cube);
    }
    if( cube[21] == 1 )
    {
      cube[21] = label;        
      Octree_labeling( 6, label, cube);
      Octree_labeling( 7, label, cube);
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==6 )
  {
  	if( cube[10] == 1 )
    {
      cube[10] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 2, label, cube);
      Octree_labeling( 5, label, cube);
    }
    if( cube[11] == 1 )
    {
      cube[11] = label;        
      Octree_labeling( 2, label, cube);
    }
    if( cube[13] == 1 )
    {
      cube[13] = label;        
      Octree_labeling( 2, label, cube);
      Octree_labeling( 4, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[18] == 1 )
    {
      cube[18] = label;        
      Octree_labeling( 5, label, cube);
    }
    if( cube[21] == 1 )
    {
      cube[21] = label;        
      Octree_labeling( 5, label, cube);
      Octree_labeling( 7, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[19] == 1 )
      cube[19] = label;        
    if( cube[22] == 1 )
    {
      cube[22] = label;        
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==7 )
  {
  	if( cube[12] == 1 )
    {
      cube[12] = label;        
      Octree_labeling( 1, label, cube);
      Octree_labeling( 3, label, cube);
      Octree_labeling( 5, label, cube);
    }
  	if( cube[14] == 1 )
    {
      cube[14] = label;        
      Octree_labeling( 3, label, cube);
    }
    if( cube[15] == 1 )
    {
      cube[15] = label;        
      Octree_labeling( 3, label, cube);
      Octree_labeling( 4, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[20] == 1 )
    {
      cube[20] = label;        
      Octree_labeling( 5, label, cube);
    }
    if( cube[21] == 1 )
    {
      cube[21] = label;        
      Octree_labeling( 5, label, cube);
      Octree_labeling( 6, label, cube);
      Octree_labeling( 8, label, cube);
    }
    if( cube[23] == 1 )
      cube[23] = label;        
    if( cube[24] == 1 )
    {
      cube[24] = label;        
      Octree_labeling( 8, label, cube);
    }
  }
  if( octant==8 )
  {
  	if( cube[13] == 1 )
    {
      cube[13] = label;        
      Octree_labeling( 2, label, cube);
      Octree_labeling( 4, label, cube);
      Octree_labeling( 6, label, cube);
    }
  	if( cube[15] == 1 )
    {
      cube[15] = label;        
      Octree_labeling( 3, label, cube);
      Octree_labeling( 4, label, cube);
      Octree_labeling( 7, label, cube);
    }
  	if( cube[16] == 1 )
    {
      cube[16] = label;        
      Octree_labeling( 4, label, cube);
    }
  	if( cube[21] == 1 )
    {
      cube[21] = label;        
      Octree_labeling( 5, label, cube);
      Octree_labeling( 6, label, cube);
      Octree_labeling( 7, label, cube);
    }
  	if( cube[22] == 1 )
    {
      cube[22] = label;        
      Octree_labeling( 6, label, cube);
    }
  	if( cube[24] == 1 )
    {
      cube[24] = label;        
      Octree_labeling( 7, label, cube);
    }
  	if( cube[25] == 1 )
      cube[25] = label;        
  } 
}




int isSimplePoint(   
	unsigned char* in_rev,		// IN: Input (binary) original volume
	const int dimx,
	const int dimy, 
	const int dimz,
	const int x,
	const int y,
	const int z
	)
{
	// copy neighbors for labeling
	int cube[26];
	int i = 0;	
	int label = 2; // set initial label

	int a,b,c;

	/*for( i = 0; i < 13; i++ )  // i =  0..12 -> cube[0..12]
    cube[i] = neighbors[i];*/

	// From 0 to 8:
	c = z - 1;
	for ( b = (y - 1); b <= (y + 1); b++ )
		for ( a = (x - 1); a <= (x + 1); a++ )
		{			
			cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					
		}

	// From 9 to 11:
	c = z;
	b = y - 1;
	for (a = (x - 1); a <= (x + 1); a++)
	{		
		cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					
	}	

	// The 12th:
	c = z;
	b = y;
	a = x - 1;

	cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					

	// i != 13 : ignore center pixel when counting (see [Lee94])
	/*for( i = 14; i < 27; i++ ) // i = 14..26 -> cube[13..25]
		cube[i-1] = neighbors[i];*/

	// The 14th:
	c = z;
	b = y;
	a = x + 1;

	cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					

	// From 15 to 18:
	c = z;
	b = y + 1;
	for (a = (x - 1); a <= (x + 1); a++)
	{		
		cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					
	}	

	// From 19 to 27:
	c = z + 1;
	for ( b = (y - 1); b <= (y + 1); b++ )
		for ( a = (x - 1); a <= (x + 1); a++ )
		{			
			cube[i++] = (in_rev[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					
		}
  
	

  // for all points in the neighborhood
  for( i = 0; i < 26; i++ )
  {
    if( cube[i]==1 )     // voxel has not been labelled yet
    {
      // start recursion with any octant that contains the point i
      switch( i )
      {
      case 0:
      case 1:
      case 3:
      case 4:
      case 9:
      case 10:
      case 12:
        Octree_labeling(1, label, cube );
        break;
      case 2:
      case 5:
      case 11:
      case 13:
        Octree_labeling(2, label, cube );
        break;
      case 6:
      case 7:
      case 14:
      case 15:
        Octree_labeling(3, label, cube );
        break;
      case 8:
      case 16:
        Octree_labeling(4, label, cube );
        break;
      case 17:
      case 18:
      case 20:
      case 21:
        Octree_labeling(5, label, cube );
        break;
      case 19:
      case 22:
        Octree_labeling(6, label, cube );
        break;
      case 23:
      case 24:
        Octree_labeling(7, label, cube );
        break;
      case 25:
        Octree_labeling(8, label, cube );
        break;
      }
      label++;
      if( label-2 >= 2 )
      {
        return P3D_FALSE;
      }
    }
  }
  //return label-2; in [Lee94] if the number of connected components would be needed
  return P3D_TRUE;
}


int p3dCrop3D_uchar2uchar (	
    unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;


	// Compute dimensions of padded REV:
	a_dimx = dimx - size*2;
	a_dimy = dimy - size*2;
	a_dimz = dimz - size*2;

	// Copy original (internal) values:
	for (k = 0; k < a_dimz; k++)
		for (j = 0; j < a_dimy; j++)
			for (i = 0; i < a_dimx; i++)			
				out_rev[ I( i, j, k, a_dimx, a_dimy ) ] = 
					in_rev[ I( i + size, j + size, k + size, dimx, dimy) ];

	// Return OK:
	return P3D_SUCCESS;
}

int p3dCrop3D_ushort2ushort (	
    unsigned short* in_rev,
	unsigned short* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;


	// Compute dimensions of padded REV:
	a_dimx = dimx - size*2;
	a_dimy = dimy - size*2;
	a_dimz = dimz - size*2;

	// Copy original (internal) values:
	for (k = 0; k < a_dimz; k++)
		for (j = 0; j < a_dimy; j++)
			for (i = 0; i < a_dimx; i++)			
				out_rev[ I( i, j, k, a_dimx, a_dimy ) ] = 
					in_rev[ I( i + size, j + size, k + size, dimx, dimy) ];

	// Return OK:
	return P3D_SUCCESS;
}

int p3dCrop3D_uint2uint (	
    unsigned int* in_rev,
	unsigned int* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;


	// Compute dimensions of padded REV:
	a_dimx = dimx - size*2;
	a_dimy = dimy - size*2;
	a_dimz = dimz - size*2;

	// Copy original (internal) values:
	for (k = 0; k < a_dimz; k++)
		for (j = 0; j < a_dimy; j++)
			for (i = 0; i < a_dimx; i++)			
				out_rev[ I( i, j, k, a_dimx, a_dimy ) ] = 
					in_rev[ I( i + size, j + size, k + size, dimx, dimy) ];

	// Return OK:
	return P3D_SUCCESS;
}

int p3dCrop3D_float2float (	
    float* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;


	// Compute dimensions of padded REV:
	a_dimx = dimx - size*2;
	a_dimy = dimy - size*2;
	a_dimz = dimz - size*2;

	// Copy original (internal) values:
	for (k = 0; k < a_dimz; k++)
		for (j = 0; j < a_dimy; j++)
			for (i = 0; i < a_dimx; i++)			
				out_rev[ I( i, j, k, a_dimx, a_dimy ) ] = 
					in_rev[ I( i + size, j + size, k + size, dimx, dimy) ];

	// Return OK:
	return P3D_SUCCESS;

}

int p3dZeroPadding3D_uchar2float (	
    unsigned char* in_rev,
	float* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size*2;
	a_dimy = dimy + size*2;
	a_dimz = dimz + size*2;

	// Set to zero all values:
	memset( out_rev, 0, a_dimx*a_dimy*a_dimz*sizeof(float) );


	// Copy original (internal) values:
	for (k = 0; k < dimz; k++)
		for (j = 0; j < dimy; j++)
			for (i = 0; i < dimx; i++)			
				out_rev[ I( i + size, j + size, k + size, a_dimx, a_dimy ) ] = 
					(float) (in_rev[ I( i, j, k, dimx, dimy ) ]);

	// Return OK:
	return P3D_SUCCESS;
}


int p3dZeroPadding3D_uchar2uchar (	
    unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;

	// Compute dimensions of padded REV:
	a_dimx = dimx + size*2;
	a_dimy = dimy + size*2;
	a_dimz = dimz + size*2;

	// Set to zero all values:
	memset( out_rev, 0, a_dimx*a_dimy*a_dimz*sizeof(unsigned char) );


	// Copy original (internal) values:
	for (k = 0; k < dimz; k++)
		for (j = 0; j < dimy; j++)
			for (i = 0; i < dimx; i++)			
				out_rev[ I( i + size, j + size, k + size, a_dimx, a_dimy ) ] = 
					in_rev[ I( i, j, k, dimx, dimy ) ];

	// Return OK:
	return P3D_SUCCESS;
}


int p3dReplicatePadding3D_uchar2uchar (	
	unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k,ct;

	
	// Compute dimensions of padded REV:
	a_dimx = dimx + size*2;
	a_dimy = dimy + size*2;
	a_dimz = dimz + size*2;


	// Perform first zero padding:
	p3dZeroPadding3D_uchar2uchar (in_rev, out_rev, dimx, dimy, dimz, size);


	// Replicate border values:
	for (ct = size; ct > 0; ct--)
	{
		// Faces:
		
		for (i = ct; i < (a_dimx - ct); i++)
			for (j = ct; j < (a_dimy - ct); j++)
			{
				out_rev[ I( i, j, ct - 1, a_dimx, a_dimy) ] = 
					out_rev[ I( i, j, ct, a_dimx, a_dimy ) ];
				out_rev[ I( i, j, a_dimz - ct, a_dimx, a_dimy)] = 
					out_rev[ I( i, j, a_dimz - 1 - ct, a_dimx, a_dimy ) ];
			}

		for (j = ct; j < (a_dimy - ct); j++)
			for (k = ct; k < (a_dimz - ct); k++)
			{
				out_rev[ I( ct - 1, j, k, a_dimx, a_dimy) ] = 
					out_rev[ I( ct, j, k, a_dimx, a_dimy ) ];
				out_rev[ I( a_dimx - ct, j, k, a_dimx, a_dimy)] = 
					out_rev[ I( a_dimx - 1 - ct, j, k, a_dimx, a_dimy ) ];
			}
		
		for (i = ct; i < (a_dimx - ct); i++)
			for (k = ct; k < (a_dimz - ct); k++)
			{
				out_rev[ I( i, ct - 1, k, a_dimx, a_dimy) ] = 
					out_rev[ I( i, ct, k, a_dimx, a_dimy ) ];
				out_rev[ I( i, a_dimy - ct, k, a_dimx, a_dimy)] = 
					out_rev[ I( i, a_dimy - 1 - ct, k, a_dimx, a_dimy ) ];
			}

		// Edges:
		
		for (i = ct; i < (a_dimx - ct); i++)
		{
			out_rev[ I( i, ct - 1, ct - 1, a_dimx, a_dimy) ] = 
				out_rev[ I( i, ct, ct, a_dimx, a_dimy ) ];
			out_rev[ I( i, ct - 1, a_dimz - ct, a_dimx, a_dimy)] = 
				out_rev[ I( i, ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];
			out_rev[ I( i, a_dimy - ct, ct - 1, a_dimx, a_dimy) ] = 
				out_rev[ I( i, a_dimy - 1 - ct, ct, a_dimx, a_dimy ) ];
			out_rev[ I( i, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] = 
				out_rev[ I( i, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];
		}

		for (j = ct; j < (a_dimy - ct); j++)
		{
			out_rev[ I( ct - 1, j, ct - 1, a_dimx, a_dimy) ] = 
				out_rev[ I( ct, j, ct, a_dimx, a_dimy ) ];
			out_rev[ I( ct - 1, j, a_dimz - ct, a_dimx, a_dimy)] = 
				out_rev[ I( ct, j, a_dimz - 1 - ct, a_dimx, a_dimy ) ];
			out_rev[ I( a_dimx - ct, j, ct - 1, a_dimx, a_dimy) ] = 
				out_rev[ I( a_dimx - 1 - ct, j, ct, a_dimx, a_dimy ) ];
			out_rev[ I( a_dimx - ct, j, a_dimz - ct, a_dimx, a_dimy)] = 
				out_rev[ I( a_dimx - 1 - ct, j, a_dimz - 1 - ct, a_dimx, a_dimy ) ];
		}

		for (k = ct; k < (a_dimz - ct); k++)
		{
			out_rev[ I( ct - 1, ct - 1, k, a_dimx, a_dimy) ] = 
				out_rev[ I( ct, ct, k, a_dimx, a_dimy ) ];
			out_rev[ I( a_dimx - ct, ct - 1, k, a_dimx, a_dimy)] = 
				out_rev[ I( a_dimx - 1 - ct, ct, k, a_dimx, a_dimy ) ];
			out_rev[ I( ct - 1, a_dimy - ct, k, a_dimx, a_dimy) ] = 
				out_rev[ I( ct, a_dimy - 1 - ct, k, a_dimx, a_dimy ) ];
			out_rev[ I( a_dimx - ct, a_dimy - ct, k, a_dimx, a_dimy)] = 
				out_rev[ I( a_dimx - 1 - ct, a_dimy - 1 - ct, k, a_dimx, a_dimy ) ];
		}

		// Corners:

		out_rev[ I( ct - 1, ct - 1, ct - 1, a_dimx, a_dimy)] = 
			out_rev[ I( ct, ct, ct, a_dimx, a_dimy ) ];

		out_rev[ I( a_dimx - ct, ct - 1, ct - 1, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - ct, ct, ct, a_dimx, a_dimy ) ];

		out_rev[ I( ct - 1, a_dimy - ct, ct - 1, a_dimx, a_dimy)] = 
			out_rev[ I( ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy ) ];

		out_rev[ I( ct - 1, ct - 1, a_dimz - ct, a_dimx, a_dimy)] = 
			out_rev[ I( ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];


		out_rev[ I( a_dimx - ct, a_dimy - ct, ct - 1, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - ct, a_dimy - 1 - ct, ct, a_dimx, a_dimy ) ];

		out_rev[ I( ct - 1, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] = 
			out_rev[ I( ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];

		out_rev[ I( a_dimx - ct, ct - 1, a_dimz - ct, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - ct, ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];

		out_rev[ I( a_dimx - ct, a_dimy - ct, a_dimz - ct, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - ct, a_dimy - 1 - ct, a_dimz - 1 - ct, a_dimx, a_dimy ) ];


	}

	// Return OK:
	return P3D_SUCCESS;
}

double interpolation ( 
	float* gvf,
	int dimx, 
	int dimy, 
	int dimz, 
	const double x, 
	const double y, 
	const double z	
	)
{
	double alpha, beta, gamma;
     
	alpha = x - (int) (x);
	beta  = y - (int) (y);
	gamma = z - (int) (z);
  
  
	return (gvf [ I( (int) (x), (int) (y), (int) (z), dimx, dimy) ]*(1-alpha)*(1-beta)*(1-gamma)
		  + gvf [ I( (int) (x), (int) (y), (int) (z+1), dimx, dimy) ]*(1-alpha)*(1-beta)*gamma
		  + gvf [ I( (int) (x), (int) (y+1), (int) (z), dimx, dimy) ]*(1-alpha)*beta*(1-gamma)
		  + gvf [ I( (int) (x+1), (int) (y), (int) (z), dimx, dimy) ]*alpha*(1-beta)*(1-gamma)
		  + gvf [ I( (int) (x+1), (int) (y), (int) (z+1), dimx, dimy) ]*alpha*(1-beta)*gamma
		  + gvf [ I( (int) (x+1), (int) (y+1), (int) (z), dimx, dimy) ]*alpha*beta*(1-gamma)
		  + gvf [ I( (int) (x), (int) (y+1), (int) (z+1), dimx, dimy) ]*(1-alpha)*beta*gamma
		  + gvf [ I( (int) (x+1), (int) (y+1), (int) (z+1), dimx, dimy) ]*(alpha*beta*gamma));	
}










int isBoundary ( 
	unsigned char* in_im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,    
	const int j,    
	const int k
	)
{
	int a,b,c;  

	c = k - 1; b = j; a = i;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	c = k + 1; b = j; a = i;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	c = k; b = j - 1; a = i;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	c = k; b = j + 1; a = i;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	c = k; b = j; a = i - 1;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	c = k; b = j; a = i + 1;
	if ( in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND ) return P3D_TRUE;

	return P3D_FALSE;
}

int isFullNeighborhood ( 
	unsigned char* in_im, 
	const int dimx, 
	const int dimy, 
	const int dimz,
	const int i,    
	const int j,    
	const int k
	)
{
	int a,b,c;   

	// Neighbors touched by interpolation:
	c = k + 1; b = j; a = i;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;

	c = k; b = j + 1; a = i;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;

	c = k; b = j; a = i + 1;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;
	
	c = k + 1; b = j + 1; a = i;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;

	c = k; b = j + 1; a = i + 1;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;

	c = k + 1; b = j; a = i + 1;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;

	c = k + 1; b = j + 1; a = i + 1;
	if ( (in_im[ I(a,b,c,dimx,dimy) ] == BACKGROUND) ||
		(isBoundary( in_im, dimx, dimy, dimz, a, b, c ) == P3D_TRUE) ) return P3D_FALSE;


	// No BACKGROUND voxel found:
	return P3D_TRUE;
}


int p3dSpecialPadding3D_uchar2uchar (	
	unsigned char* in_rev,
	unsigned char* out_rev,
	const int dimx, // ncols
	const int dimy, // nrows
	const int dimz, // nplanes
	const int size
	)
{
	int a_dimx, a_dimy, a_dimz;
	int i,j,k;

	
	// Compute dimensions of padded REV:
	a_dimx = dimx + size*2;
	a_dimy = dimy + size*2;
	a_dimz = dimz + size*2;


	// Set to zero all values:
	memset( out_rev, 0, a_dimx*a_dimy*a_dimz*sizeof(unsigned char) );


	// Copy original (internal) values converting 255 to 1:
	for (k = 0; k < dimz; k++)
		for (j = 0; j < dimy; j++)
			for (i = 0; i < dimx; i++)			
				if ( in_rev[ I( i, j, k, dimx, dimy ) ] == OBJECT ) 
					out_rev[ I( i + size, j + size, k + size, a_dimx, a_dimy ) ] = 1;



	// Replicate border values:

	// Faces:		
	/*for (i = size; i < (a_dimx - size); i++)
		for (j = size; j < (a_dimy - size); j++)
		{
			out_rev[ I( i, j, size - 1, a_dimx, a_dimy) ] = 
				out_rev[ I( i, j, size, a_dimx, a_dimy ) ];
			out_rev[ I( i, j, a_dimz - size, a_dimx, a_dimy)] = 
				out_rev[ I( i, j, a_dimz - 1 - size, a_dimx, a_dimy ) ];
		}

	for (j = size; j < (a_dimy - size); j++)
		for (k = size; k < (a_dimz - size); k++)
		{
			out_rev[ I( size - 1, j, k, a_dimx, a_dimy) ] = 
				out_rev[ I( size, j, k, a_dimx, a_dimy ) ];
			out_rev[ I( a_dimx - size, j, k, a_dimx, a_dimy)] = 
				out_rev[ I( a_dimx - 1 - size, j, k, a_dimx, a_dimy ) ];
		}
	
	for (i = size; i < (a_dimx - size); i++)
		for (k = size; k < (a_dimz - size); k++)
		{
			out_rev[ I( i, size - 1, k, a_dimx, a_dimy) ] = 
				out_rev[ I( i, size, k, a_dimx, a_dimy ) ];
			out_rev[ I( i, a_dimy - size, k, a_dimx, a_dimy)] = 
				out_rev[ I( i, a_dimy - 1 - size, k, a_dimx, a_dimy ) ];
		}

	// Edges:
	
	for (i = size; i < (a_dimx - size); i++)
	{
		out_rev[ I( i, size - 1, size - 1, a_dimx, a_dimy) ] = 
			out_rev[ I( i, size, size, a_dimx, a_dimy ) ];
		out_rev[ I( i, size - 1, a_dimz - size, a_dimx, a_dimy)] = 
			out_rev[ I( i, size, a_dimz - 1 - size, a_dimx, a_dimy ) ];
		out_rev[ I( i, a_dimy - size, size - 1, a_dimx, a_dimy) ] = 
			out_rev[ I( i, a_dimy - 1 - size, size, a_dimx, a_dimy ) ];
		out_rev[ I( i, a_dimy - size, a_dimz - size, a_dimx, a_dimy)] = 
			out_rev[ I( i, a_dimy - 1 - size, a_dimz - 1 - size, a_dimx, a_dimy ) ];
	}

	for (j = size; j < (a_dimy - size); j++)
	{
		out_rev[ I( size - 1, j, size - 1, a_dimx, a_dimy) ] = 
			out_rev[ I( size, j, size, a_dimx, a_dimy ) ];
		out_rev[ I( size - 1, j, a_dimz - size, a_dimx, a_dimy)] = 
			out_rev[ I( size, j, a_dimz - 1 - size, a_dimx, a_dimy ) ];
		out_rev[ I( a_dimx - size, j, size - 1, a_dimx, a_dimy) ] = 
			out_rev[ I( a_dimx - 1 - size, j, size, a_dimx, a_dimy ) ];
		out_rev[ I( a_dimx - size, j, a_dimz - size, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - size, j, a_dimz - 1 - size, a_dimx, a_dimy ) ];
	}

	for (k = size; k < (a_dimz - size); k++)
	{
		out_rev[ I( size - 1, size - 1, k, a_dimx, a_dimy) ] = 
			out_rev[ I( size, size, k, a_dimx, a_dimy ) ];
		out_rev[ I( a_dimx - size, size - 1, k, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - size, size, k, a_dimx, a_dimy ) ];
		out_rev[ I( size - 1, a_dimy - size, k, a_dimx, a_dimy) ] = 
			out_rev[ I( size, a_dimy - 1 - size, k, a_dimx, a_dimy ) ];
		out_rev[ I( a_dimx - size, a_dimy - size, k, a_dimx, a_dimy)] = 
			out_rev[ I( a_dimx - 1 - size, a_dimy - 1 - size, k, a_dimx, a_dimy ) ];
	}

	// Corners:

	out_rev[ I( size - 1, size - 1, size - 1, a_dimx, a_dimy)] = 
		out_rev[ I( size, size, size, a_dimx, a_dimy ) ];

	out_rev[ I( a_dimx - size, size - 1, size - 1, a_dimx, a_dimy)] = 
		out_rev[ I( a_dimx - 1 - size, size, size, a_dimx, a_dimy ) ];

	out_rev[ I( size - 1, a_dimy - size, size - 1, a_dimx, a_dimy)] = 
		out_rev[ I( size, a_dimy - 1 - size, size, a_dimx, a_dimy ) ];

	out_rev[ I( size - 1, size - 1, a_dimz - size, a_dimx, a_dimy)] = 
		out_rev[ I( size, size, a_dimz - 1 - size, a_dimx, a_dimy ) ];


	out_rev[ I( a_dimx - size, a_dimy - size, size - 1, a_dimx, a_dimy)] = 
		out_rev[ I( a_dimx - 1 - size, a_dimy - 1 - size, size, a_dimx, a_dimy ) ];

	out_rev[ I( size - 1, a_dimy - size, a_dimz - size, a_dimx, a_dimy)] = 
		out_rev[ I( size, a_dimy - 1 - size, a_dimz - 1 - size, a_dimx, a_dimy ) ];

	out_rev[ I( a_dimx - size, size - 1, a_dimz - size, a_dimx, a_dimy)] = 
		out_rev[ I( a_dimx - 1 - size, size, a_dimz - 1 - size, a_dimx, a_dimy ) ];

	out_rev[ I( a_dimx - size, a_dimy - size, a_dimz - size, a_dimx, a_dimy)] = 
		out_rev[ I( a_dimx - 1 - size, a_dimy - 1 - size, a_dimz - 1 - size, a_dimx, a_dimy ) ];*/
	

	// Return OK:
	return P3D_SUCCESS;
}
