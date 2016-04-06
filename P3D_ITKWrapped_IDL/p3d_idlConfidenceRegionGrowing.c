
// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dITKWrapped.h"


IDL_VPTR p3d_idlConfidenceRegionGrowing (int argc, IDL_VPTR argv[], char* argk )  
{  
	typedef struct {  
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		int iterations;  
		int it_there;
		double multiplier;  
		int mu_there;  
		IDL_LONG seed_data[3];    // Maximum size
		int seed_there;  
		IDL_MEMINT seed_n;  
		int winsize;  
		int ws_there;
	} KW_RESULT;  

	static IDL_KW_ARR_DESC_R seed = { (char*) IDL_KW_OFFSETOF(seed_data), 2, 3, (IDL_LONG*) IDL_KW_OFFSETOF(seed_n) };
   
	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {  
		IDL_KW_FAST_SCAN,  
		{ "ITERATIONS", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(it_there), (char*) IDL_KW_OFFSETOF(iterations) },  				
		{ "K", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(mu_there), (char*) IDL_KW_OFFSETOF(multiplier) },  				
		{ "SEED", IDL_TYP_LONG, 1, IDL_KW_ARRAY, (int*) IDL_KW_OFFSETOF(seed_there), (char*) (&seed) },  						
		{ "WIDTH", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(ws_there), (char*) IDL_KW_OFFSETOF(winsize) },  		
    	{ NULL }  
	 };  

	KW_RESULT kw;
	
	IDL_VPTR idl_out_rev, idl_in_rev;
	unsigned char *in_rev8, *out_rev;  	
	unsigned short *in_rev16;

	// Keywords:
	double multiplier = 2.5;     // default
	unsigned int iterations = 0; // default
	unsigned int winsize = 3; // default
	unsigned int seedX, seedY, seedZ;
	
	// Keywords counter:
	int keywords_ct = 0;

	int err_code;
	


	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 

	
	// Get input data in IDL format:
	idl_in_rev = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_rev);  
	IDL_ENSURE_ARRAY(idl_in_rev);  
	
	
	// Get the ITERATIONS input argument:
	if (kw.it_there)
	{
		// Check values:
		if ( ( kw.iterations < 0 ) || ( kw.multiplier > 50 ) )
			_p3d_idlPrintNamedError("ITERATIONS must be within the range [1,40].");

		// Get values:
		iterations = (unsigned int) kw.iterations;

		keywords_ct++;
	}

	// Get the MULTIPLIER input argument:
	if (kw.mu_there)
	{
		// Check values:
		if ( ( kw.multiplier < 1 ) || ( kw.multiplier > 50 ) )
			_p3d_idlPrintNamedError("MULTIPLIER must be within the range [1,40].");

		// Get values:
		multiplier = (double) kw.multiplier;

		keywords_ct++;
	}

	// Get the WINSIZE input argument:
	if (kw.ws_there)
	{
		// Check values:
		if ( ( kw.winsize < 3 ) || ( kw.multiplier > 51 ) )
			_p3d_idlPrintNamedError("WINSIZE must be within the range [3,51].");

		// Get values:
		winsize = (unsigned int) kw.winsize;

		keywords_ct++;
	}


	// Allocate memory for output:
	if (!(idl_in_rev->flags & IDL_V_TEMP))  
		out_rev = (unsigned char *) IDL_MakeTempArray(
			IDL_TYP_BYTE,
			idl_in_rev->value.arr->n_dim,  
			idl_in_rev->value.arr->dim,  
			IDL_ARR_INI_NOP, 
			&idl_out_rev
			);    

	// Call Pore3D depending on input arguments:
	if ( idl_in_rev->value.arr->n_dim == 2 )
	{	
		// Get SEED input argument:
		if (kw.seed_there)
		{
			// Check values:
			if ( kw.seed_n != 2 )
				_p3d_idlPrintNamedError("Input argument SEED must contain two [ X, Y ] elements.");

			if ( ( kw.seed_data[0] < 0 ) || ( kw.seed_data[0] > idl_in_rev->value.arr->dim[0] ) )
				_p3d_idlPrintNamedError("X value of input argument SEED must be within IMAGE dimensions.");

			if ( ( kw.seed_data[1] < 0 ) || ( kw.seed_data[1] > idl_in_rev->value.arr->dim[1] ) )
				_p3d_idlPrintNamedError("Y value of input argument SEED must be within IMAGE dimensions.");

			// Get values:
			seedX = (unsigned int) kw.seed_data[0];
			seedY = (unsigned int) kw.seed_data[1];

			keywords_ct++;
		}
		else
		{
			// Set default values for centerX and centerY:
			seedX = 0;
			seedY = 0;
		}


		// Extract first input (volume to filter) in C format:
		if (idl_in_rev->type == IDL_TYP_BYTE)			
		{
			in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;  
 

			// Call Pore3D without the mask:
			err_code = p3dConfidenceConnectedRegionGrowing2D_8 ( 
				in_rev8,
				out_rev,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				multiplier,
				winsize,
				iterations,
				seedX,
				seedY,
				_p3d_idlPrintInfo
			);	

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else if (idl_in_rev->type == IDL_TYP_UINT)		
		{
			in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;  
  

			// Call Pore3D without the mask:
			err_code = p3dConfidenceConnectedRegionGrowing2D_16 ( 
				in_rev16,
				out_rev,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				multiplier,
				winsize,
				iterations,
				seedX,
				seedY,
				_p3d_idlPrintInfo
			);	

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else
		{
			_p3d_idlPrintNamedError("Input argument VOLUME must be of type BYTE or UINT.");
		}		
	}
	else if ( idl_in_rev->value.arr->n_dim == 3 )
	{	
		// Get SEED input argument:
		if (kw.seed_there)
		{
			// Check values:
			if ( kw.seed_n != 3)
				_p3d_idlPrintNamedError("Input argument SEED must contain three [ X, Y, Z ] elements.");

			if ( ( kw.seed_data[0] < 0 ) || ( kw.seed_data[0] > idl_in_rev->value.arr->dim[0] ) )
				_p3d_idlPrintNamedError("X value of input argument SEED must be within IMAGE dimensions.");

			if ( ( kw.seed_data[1] < 0 ) || ( kw.seed_data[1] > idl_in_rev->value.arr->dim[1] ) )
				_p3d_idlPrintNamedError("Y value of input argument SEED must be within IMAGE dimensions.");

			if ( ( kw.seed_data[2] < 0 ) || ( kw.seed_data[2] > idl_in_rev->value.arr->dim[2] ) )
				_p3d_idlPrintNamedError("Z value of input argument SEED must be within IMAGE dimensions.");

			// Get values:
			seedX = (unsigned int) kw.seed_data[0];
			seedY = (unsigned int) kw.seed_data[1];
			seedZ = (unsigned int) kw.seed_data[2];

			keywords_ct++;
		}
		else
		{
			// Set default values for centerX and centerY:
			seedX = 0;
			seedY = 0;
			seedZ = 0;
		}


		// Extract first input (volume to filter) in C format:
		if (idl_in_rev->type == IDL_TYP_BYTE)			
		{
			in_rev8 = (unsigned char *) idl_in_rev->value.arr->data; 
 

			// Call Pore3D:
			err_code = p3dConfidenceConnectedRegionGrowing3D_8 ( 
				in_rev8,
				out_rev,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				multiplier,
				winsize,
				iterations,
				seedX,
				seedY,
				seedZ,
				_p3d_idlPrintInfo
			);	

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else if (idl_in_rev->type == IDL_TYP_UINT)		
		{
			in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;  
 

			// Call Pore3D without the mask:
			err_code = p3dConfidenceConnectedRegionGrowing3D_16 ( 
				in_rev16,
				out_rev,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				multiplier,
				winsize,
				iterations,
				seedX,
				seedY,
				seedZ,
				_p3d_idlPrintInfo
			);	

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else
		{
			_p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
		}		
	}
	else
	{
		_p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
	}

	// Free resources:
	IDL_KW_FREE;

	// Return output in IDL Format
	return(idl_out_rev);
}  


