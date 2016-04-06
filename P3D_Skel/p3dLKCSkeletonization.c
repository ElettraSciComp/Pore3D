#include <stdlib.h>
#include <omp.h>

#include "p3dSkel.h"
#include "p3dTime.h"

#include "Common/p3dCoordsList.h"
#include "Common/p3dUtils.h"


void fillEulerLUT ( int *LUT )
{
	LUT[1]  =  1;
	LUT[3]  = -1;
	LUT[5]  = -1;
	LUT[7]  =  1;
	LUT[9]  = -3;
	LUT[11] = -1;
	LUT[13] = -1;
	LUT[15] =  1;
	LUT[17] = -1;
	LUT[19] =  1;
	LUT[21] =  1;
	LUT[23] = -1;
	LUT[25] =  3;
	LUT[27] =  1;
	LUT[29] =  1;
	LUT[31] = -1;
	LUT[33] = -3;
	LUT[35] = -1;
	LUT[37] =  3;
	LUT[39] =  1;
	LUT[41] =  1;
	LUT[43] = -1;
	LUT[45] =  3;
	LUT[47] =  1;
	LUT[49] = -1;
	LUT[51] =  1;

	LUT[53] =  1;
	LUT[55] = -1;
	LUT[57] =  3;
	LUT[59] =  1;
	LUT[61] =  1;
	LUT[63] = -1;
	LUT[65] = -3;
	LUT[67] =  3;
	LUT[69] = -1;
	LUT[71] =  1;
	LUT[73] =  1;
	LUT[75] =  3;
	LUT[77] = -1;
	LUT[79] =  1;
	LUT[81] = -1;
	LUT[83] =  1;
	LUT[85] =  1;
	LUT[87] = -1;
	LUT[89] =  3;
	LUT[91] =  1;
	LUT[93] =  1;
	LUT[95] = -1;
	LUT[97] =  1;
	LUT[99] =  3;
	LUT[101] =  3;
	LUT[103] =  1;

	LUT[105] =  5;
	LUT[107] =  3;
	LUT[109] =  3;
	LUT[111] =  1;
	LUT[113] = -1;
	LUT[115] =  1;
	LUT[117] =  1;
	LUT[119] = -1;
	LUT[121] =  3;
	LUT[123] =  1;
	LUT[125] =  1;
	LUT[127] = -1;
	LUT[129] = -7;
	LUT[131] = -1;
	LUT[133] = -1;
	LUT[135] =  1;
	LUT[137] = -3;
	LUT[139] = -1;
	LUT[141] = -1;
	LUT[143] =  1;
	LUT[145] = -1;
	LUT[147] =  1;
	LUT[149] =  1;
	LUT[151] = -1;
	LUT[153] =  3;
	LUT[155] =  1;

	LUT[157] =  1;
	LUT[159] = -1;
	LUT[161] = -3;
	LUT[163] = -1;
	LUT[165] =  3;
	LUT[167] =  1;
	LUT[169] =  1;
	LUT[171] = -1;
	LUT[173] =  3;
	LUT[175] =  1;
	LUT[177] = -1;
	LUT[179] =  1;
	LUT[181] =  1;
	LUT[183] = -1;
	LUT[185] =  3;
	LUT[187] =  1;
	LUT[189] =  1;
	LUT[191] = -1;
	LUT[193] = -3;
	LUT[195] =  3;
	LUT[197] = -1;
	LUT[199] =  1;
	LUT[201] =  1;
	LUT[203] =  3;
	LUT[205] = -1;
	LUT[207] =  1;

	LUT[209] = -1;
	LUT[211] =  1;
	LUT[213] =  1;
	LUT[215] = -1;
	LUT[217] =  3;
	LUT[219] =  1;
	LUT[221] =  1;
	LUT[223] = -1;
	LUT[225] =  1;
	LUT[227] =  3;
	LUT[229] =  3;
	LUT[231] =  1;
	LUT[233] =  5;
	LUT[235] =  3;
	LUT[237] =  3;
	LUT[239] =  1;
	LUT[241] = -1;
	LUT[243] =  1;
	LUT[245] =  1;
	LUT[247] = -1;
	LUT[249] =  3;
	LUT[251] =  1;
	LUT[253] =  1;
	LUT[255] = -1;
}

/** 
 * Check for Euler invariance. (see [Lee94])
 */
int isEulerInvariant(
	unsigned char* in_im,		// IN: Input (binary) original volume
	const int dimx,
	const int dimy, 
	const int dimz,
	const int x,
	const int y,
	const int z,
	const int *LUT
	)
{
	// calculate Euler characteristic for each octant and sum up
	int EulerChar = 0;
	unsigned char n;

	// copy neighbors for labeling
	int neighbors[27];
	int a,b,c;
	int i = 0;

	// Fill the neieghbors:
	for ( c = (z - 1); c <= (z + 1); c++ )
		for ( b = (y - 1); b <= (y + 1); b++ )
			for ( a = (x - 1); a <= (x + 1); a++ )
			{			
				neighbors[i++] = (in_im[ I( a, b, c, dimx, dimy ) ] == OBJECT) ? 1 : 0;					
			}


	// Octant SWU
	n = 1;
	if( neighbors[24]==1 )
		n |= 128;
	if( neighbors[25]==1 )
		n |=  64;
	if( neighbors[15]==1 )
		n |=  32;
	if( neighbors[16]==1 )
		n |=  16;
	if( neighbors[21]==1 )
		n |=   8;
	if( neighbors[22]==1 )
		n |=   4;
	if( neighbors[12]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant SEU
	n = 1;
	if( neighbors[26]==1 )
		n |= 128;
	if( neighbors[23]==1 )
		n |=  64;
	if( neighbors[17]==1 )
		n |=  32;
	if( neighbors[14]==1 )
		n |=  16;
	if( neighbors[25]==1 )
		n |=   8;
	if( neighbors[22]==1 )
		n |=   4;
	if( neighbors[16]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant NWU
	n = 1;
	if( neighbors[18]==1 )
		n |= 128;
	if( neighbors[21]==1 )
		n |=  64;
	if( neighbors[9]==1 )
		n |=  32;
	if( neighbors[12]==1 )
		n |=  16;
	if( neighbors[19]==1 )
		n |=   8;
	if( neighbors[22]==1 )
		n |=   4;
	if( neighbors[10]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant NEU
	n = 1;
	if( neighbors[20]==1 )
		n |= 128;
	if( neighbors[23]==1 )
		n |=  64;
	if( neighbors[19]==1 )
		n |=  32;
	if( neighbors[22]==1 )
		n |=  16;
	if( neighbors[11]==1 )
		n |=   8;
	if( neighbors[14]==1 )
		n |=   4;
	if( neighbors[10]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant SWB
	n = 1;
	if( neighbors[6]==1 )
		n |= 128;
	if( neighbors[15]==1 )
		n |=  64;
	if( neighbors[7]==1 )
		n |=  32;
	if( neighbors[16]==1 )
		n |=  16;
	if( neighbors[3]==1 )
		n |=   8;
	if( neighbors[12]==1 )
		n |=   4;
	if( neighbors[4]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant SEB
	n = 1;
	if( neighbors[8]==1 )
		n |= 128;
	if( neighbors[7]==1 )
		n |=  64;
	if( neighbors[17]==1 )
		n |=  32;
	if( neighbors[16]==1 )
		n |=  16;
	if( neighbors[5]==1 )
		n |=   8;
	if( neighbors[4]==1 )
		n |=   4;
	if( neighbors[14]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant NWB
	n = 1;
	if( neighbors[0]==1 )
		n |= 128;
	if( neighbors[9]==1 )
		n |=  64;
	if( neighbors[3]==1 )
		n |=  32;
	if( neighbors[12]==1 )
		n |=  16;
	if( neighbors[1]==1 )
		n |=   8;
	if( neighbors[10]==1 )
		n |=   4;
	if( neighbors[4]==1 )
		n |=   2;
	EulerChar += LUT[n];
	// Octant NEB
	n = 1;
	if( neighbors[2]==1 )
		n |= 128;
	if( neighbors[1]==1 )
		n |=  64;
	if( neighbors[11]==1 )
		n |=  32;
	if( neighbors[10] == 1 )
		n |=  16;
	if( neighbors[5]  == 1 )
		n |=   8;
	if( neighbors[4]  == 1 )
		n |=   4;
	if( neighbors[14] == 1 )
		n |=   2;
	EulerChar += LUT[n];

	if( EulerChar == 0 )
		return P3D_TRUE;
	else
		return P3D_FALSE;
}

int p3dLKCSkeletonization(   
	unsigned char* in_im, 
	unsigned char* out_im, 
	const int dimx,
	const int dimy, 
	const int dimz,
	int (*wr_log)(const char*, ...)
	)
{

	unsigned char* tmp_im;	
	
	// Define Euler LUT:
	int* eulerLUT; 


	// Padding/cropping size and dimensions
	const int a_rad = 1;					
	const int a_dimx = dimx + a_rad*2;
	const int a_dimy = dimy + a_rad*2;
	const int a_dimz = dimz + a_rad*2;	

	int unchangedBorders;
	int isBorderPoint;
	int currentBorder;
	int noChange;

	int i,j,k;
	int a,b,c;

	coords_t       simple_point;
	coords_list_t  simple_point_list;

    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dLKCSkeletonization");
    if (auth_code == '0') goto AUTH_ERROR;*/
        


	// Start tracking computational time:
	if (wr_log != NULL)
	{
		p3dResetStartTime(); 
		wr_log ("Pore3D - Performing LKC skeletonization..." );
	}



	// Init output voume with input volume values zero padded:
	P3D_TRY( eulerLUT = (int*) calloc(256,sizeof(int)) );
	P3D_TRY( tmp_im = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) ) );
	P3D_TRY( p3dZeroPadding3D_uchar2uchar ( in_im, tmp_im, dimx, dimy, dimz, a_rad ) );

	// Prepare Euler LUT:
	fillEulerLUT( eulerLUT );
  
	
    // Loop through the image several times until there is no change.

	unchangedBorders = 0;

	// Loop until no change for all the six border types:
	while( unchangedBorders < 6 )  
	{
		unchangedBorders = 0;

		#pragma omp parallel for private(i, j, k, simple_point_list, isBorderPoint, a, b, c, simple_point, noChange) reduction ( + :  unchangedBorders)
		for( currentBorder = 1; currentBorder <= 6; currentBorder++)
		{
			coords_list_init ( &simple_point_list );

			// Loop through the image:
			for (k = a_rad; k < (a_dimz - a_rad); k++)
				for (j = a_rad; j < (a_dimy - a_rad); j++)
					for (i = a_rad; i < (a_dimx - a_rad); i++)		
					{ 
						// Check if point is foreground:
						if ( tmp_im[ I(i,j,k,a_dimx,a_dimy)] == BACKGROUND)
						{
							continue;
						}

						// Check 6-neighbors if point is a border point of type currentBorder
						isBorderPoint = P3D_FALSE;
						if( (currentBorder == 1) && ( tmp_im[ I(i - 1,j,k,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;
						if( (currentBorder == 2) && ( tmp_im[ I(i + 1,j,k,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;
						if( (currentBorder == 3) && ( tmp_im[ I(i,j - 1,k,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;
						if( (currentBorder == 4) && ( tmp_im[ I(i,j + 1,k,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;
						if( (currentBorder == 5) && ( tmp_im[ I(i,j,k - 1,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;
						if( (currentBorder == 6) && ( tmp_im[ I(i,j,k + 1,a_dimx,a_dimy)] == BACKGROUND ) )
						  isBorderPoint = P3D_TRUE;

						if( isBorderPoint == P3D_FALSE)
						{
							continue;         // current point is not deletable
						}        
		        
						// Check if point is the end of an arc:
						if ( countNeighbors( tmp_im, a_dimx, a_dimy, a_dimz, i, j, k ) == 1 )
						{
							continue;
						}

						// Check if point is Euler invariant:
						if( isEulerInvariant( tmp_im, a_dimx, a_dimy, a_dimz, i, j, k, eulerLUT ) == P3D_FALSE)
						{
							continue;         // current point is not deletable
						}

						// Check if point is simple:
						if( isSimplePoint( tmp_im, a_dimx, a_dimy, a_dimz, i, j, k ) == P3D_FALSE )
						{
							continue;         // current point is not deletable
						}

						// Add current simple border point to a list for 
						// further sequential re-checking:
						simple_point.x = i;
						simple_point.y = j;
						simple_point.z = k;

						coords_list_push ( &simple_point_list, simple_point );

					} // end image iteration loop

			// sequential re-checking to preserve connectivity when
			// deleting in a parallel way
			noChange = P3D_TRUE;
			
			while ( coords_list_isempty( simple_point_list ) == P3D_FALSE )			
			{
      			// 1. Set simple border point to BACKGROUND:
				simple_point = coords_list_pop ( &simple_point_list);

				a = simple_point.x;
				b = simple_point.y;
				c = simple_point.z;

				tmp_im[ I(a,b,c,a_dimx,a_dimy) ] = BACKGROUND;
				
				// 2. Check if neighborhood is still connected:
				if ( isSimplePoint ( tmp_im, a_dimx, a_dimy, a_dimz, a, b, c ) == P3D_FALSE )
				{
					// We cannot delete current point, so reset:
					tmp_im[ I(a,b,c,a_dimx,a_dimy) ] = OBJECT;
				}
				else
				{
					noChange = P3D_FALSE;
				}
			}

			if( noChange == P3D_TRUE )
				unchangedBorders++;

		} // end currentBorder for loop
	} // end unchangedBorders while loop


	// Crop output:
	P3D_TRY( p3dCrop3D_uchar2uchar ( tmp_im, out_im, a_dimx, a_dimy, a_dimz, a_rad ) );



	// Print elapsed time (if required):
	if (wr_log != NULL)
	{			
		wr_log ("Pore3D - LKC skeletonization performed successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
	}	
	
    // Release resources:	
	if ( eulerLUT != NULL ) free ( eulerLUT );
	if ( tmp_im != NULL ) free( tmp_im );	

	// Return OK:
	return P3D_SUCCESS;
	

MEM_ERROR:

	if (wr_log != NULL)
	{
		wr_log ("Pore3D - Not enough (contiguous) memory. Program will exit.");
	}

    // Release resources:	
	if ( tmp_im != NULL ) free(tmp_im);	

	
	return P3D_MEM_ERROR;
        
        /*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
}

