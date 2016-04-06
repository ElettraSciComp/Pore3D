
// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dITKWrapped.h"


IDL_VPTR p3d_idlMultipleOtsuThresholding(int argc, IDL_VPTR argv[], char* argk)
{ 
	IDL_VPTR idl_out_rev, idl_in_rev;
	unsigned char *in_rev8, *out_rev8;  	
	unsigned short *in_rev16;

	unsigned int regions;

	int err_code;
		
	// Get input data in IDL format:
	idl_in_rev = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_rev);  
	IDL_ENSURE_ARRAY(idl_in_rev);  

	// 
	//Extract inputs in C format:
	//
	if ( argv[1]->type == IDL_TYP_INT)
		regions = (unsigned int) argv[1]->value.i;
	else 
		_p3d_idlPrintNamedError("Input argument REGIONS must be of type INT.");

	
	// Call Pore3D depending on input arguments:
	if ( idl_in_rev->value.arr->n_dim == 2 )
	{	
		// Extract first input (volume to filter) in C format:
		if (idl_in_rev->type == IDL_TYP_BYTE)			
		{
			in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_rev->flags & IDL_V_TEMP))  
				out_rev8 = (unsigned char *) IDL_MakeTempArray(
					IDL_TYP_BYTE,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					); 
 

			// Call Pore3D:
			err_code = p3dMultipleOtsuThresholding2D_8 ( 
				in_rev8,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				regions,
				_p3d_idlPrintInfo
			);

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else if (idl_in_rev->type == IDL_TYP_UINT)		
		{
			in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_rev->flags & IDL_V_TEMP))  
				out_rev8 = (unsigned char *) IDL_MakeTempArray(
					IDL_TYP_BYTE,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					);   
 
			
			// Call Pore3D:
			err_code = p3dMultipleOtsuThresholding2D_16 ( 
				in_rev16,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				regions,
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
	else if ( idl_in_rev->value.arr->n_dim == 3 )
	{
		// Extract first input (volume to filter) in C format:
		if (idl_in_rev->type == IDL_TYP_BYTE)			
		{
			in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_rev->flags & IDL_V_TEMP))  
				out_rev8 = (unsigned char *) IDL_MakeTempArray(
					IDL_TYP_BYTE,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					); 
 
			
			// Call Pore3D without the mask:
			err_code = p3dMultipleOtsuThresholding3D_8 ( 
				in_rev8,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				regions,
				_p3d_idlPrintInfo
			);	

			// On exception print error:
			if (err_code == P3D_ERROR)
				_p3d_idlPrintNamedError("Error on code execution.");
 
		}
		else if (idl_in_rev->type == IDL_TYP_UINT)		
		{
			in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_rev->flags & IDL_V_TEMP))  
				out_rev8 = (unsigned char *) IDL_MakeTempArray(
					IDL_TYP_BYTE,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					); 
 
			
			// Call Pore3D without the mask:
			err_code = p3dMultipleOtsuThresholding3D_16 ( 
				in_rev16,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				regions,
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

	// Return output in IDL Format
	return(idl_out_rev);  
}  


