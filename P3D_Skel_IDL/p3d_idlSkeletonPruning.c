// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dSkel.h"


IDL_VPTR p3d_idlSkeletonPruning (int argc, IDL_VPTR argv[], char* argk)  
{  
	typedef struct {  
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		IDL_LONG iterative;
		IDL_LONG thresh;
		int thresh_there;
		IDL_LONG ultimate;
	} KW_RESULT;  

   	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {  
		IDL_KW_FAST_SCAN,   
		{ "ITERATIVE", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(iterative) },
		{ "THRESH", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(thresh_there), (char*) IDL_KW_OFFSETOF(thresh) }, 
		{ "ULTIMATE", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(ultimate) },
		{ NULL }  
	 };  

	KW_RESULT kw;
	
	IDL_VPTR idl_out_im, idl_in_im;
	unsigned char *in_im8, *out_im8;  		
	int keywords_ct = 0;

	int thresh  = 3; // default value
	int ultimate = P3D_FALSE;
	int iterative = P3D_FALSE;
		
	int err_code;
	
	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 

	// Get input data in IDL format:
	idl_in_im = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_im);  
	IDL_ENSURE_ARRAY(idl_in_im);  


	// Get the MU keyword:
	if (kw.thresh_there)
	{
		// Check values:
		if ( kw.thresh <= 1 )
			_p3d_idlPrintNamedError("THRESH must be greater than 1.");

		// Get values:
		thresh = (int) kw.thresh;
                

		keywords_ct++;
	}
	// Get the ITERATIVE keyword:
	iterative  = kw.iterative ? P3D_TRUE : P3D_FALSE;
	// Get the ULTIMATE keyword:
	ultimate   = kw.ultimate ? P3D_TRUE : P3D_FALSE;
	

	// Call Pore3D depending on input arguments:
	if ( idl_in_im->value.arr->n_dim == 3 )
	{	
		// Extract first input (volume to filter) in C format:
		if (idl_in_im->type == IDL_TYP_BYTE)			
		{
			in_im8 = (unsigned char *) idl_in_im->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_im->flags & IDL_V_TEMP))  
				out_im8 = (unsigned char *) IDL_MakeTempArray(
					IDL_TYP_BYTE,
					idl_in_im->value.arr->n_dim,  
					idl_in_im->value.arr->dim,  
					IDL_ARR_INI_ZERO, 
					&idl_out_im
					);  

			if ( ( iterative == P3D_FALSE ) && ( ultimate == P3D_FALSE ) )
			{
				// Call Pore3D:
				err_code = p3dSimpleSkeletonPruning ( 
					in_im8,
					out_im8,
					(int) idl_in_im->value.arr->dim[0],
					(int) idl_in_im->value.arr->dim[1],
					(int) idl_in_im->value.arr->dim[2],
					thresh,
					_p3d_idlPrintInfo
				);
			} 
			if ( ( iterative == P3D_TRUE ) && ( ultimate == P3D_FALSE ) )
			{
				// Call Pore3D:
				err_code = p3dIterativeSkeletonPruning ( 
					in_im8,
					out_im8,
					(int) idl_in_im->value.arr->dim[0],
					(int) idl_in_im->value.arr->dim[1],
					(int) idl_in_im->value.arr->dim[2],
					thresh,
					_p3d_idlPrintInfo
				);
			}
			if ( ultimate == P3D_TRUE )
			{
				// Call Pore3D:
				err_code = p3dUltimateSkeletonPruning ( 
					in_im8,
					out_im8,
					(int) idl_in_im->value.arr->dim[0],
					(int) idl_in_im->value.arr->dim[1],
					(int) idl_in_im->value.arr->dim[2],
					iterative,
					_p3d_idlPrintInfo
				);
			}

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}		
		else
		{
			_p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");
		}		
	}
	else
	{
		_p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
	}

		
	// Free keywords resources:
	IDL_KW_FREE;

	// Return output in IDL format:
	return idl_out_im;  
} 

