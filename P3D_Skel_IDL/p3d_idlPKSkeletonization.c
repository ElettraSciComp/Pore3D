// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dSkel.h"


IDL_VPTR p3d_idlPKSkeletonization (int argc, IDL_VPTR argv[], char* argk)  
{  
	IDL_VPTR idl_out_im, idl_in_im;
	unsigned char *in_im8, *out_im8;  		
	int keywords_ct = 0;	
		
	int err_code;
	

	// Get input data in IDL format:
	idl_in_im = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_im);  
	IDL_ENSURE_ARRAY(idl_in_im);  
	

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
 

			// Call Pore3D without the mask:
			err_code = p3dThinningSkeletonization ( 
				in_im8,
				out_im8,
				(int) idl_in_im->value.arr->dim[0],
				(int) idl_in_im->value.arr->dim[1],
				(int) idl_in_im->value.arr->dim[2],				
				_p3d_idlPrintInfo
			);

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

	// Return output in IDL format:
	return idl_out_im;  
} 

