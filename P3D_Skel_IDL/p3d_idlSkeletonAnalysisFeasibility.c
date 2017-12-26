// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dSkel.h"


IDL_VPTR p3d_idlSkeletonAnalysisFeasibility(int argc, IDL_VPTR argv[], char* argk)  
{  
		
	IDL_VPTR idl_in_rev, idl_skl_rev, idl_out_ratio;  
	unsigned char *in_rev, *skl_rev;
	double* ratio;

	int err_code;

	// Get input data in IDL format:
	idl_in_rev = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_rev);  
	IDL_ENSURE_ARRAY(idl_in_rev);  // Get input data in IDL format:		

	// Get input data in IDL format:
	idl_skl_rev = argv[1]; 

	IDL_ENSURE_SIMPLE(idl_skl_rev);  
	IDL_ENSURE_ARRAY(idl_skl_rev);  // Get input data in IDL format:	

	// Check input format:
	if (idl_in_rev->value.arr->n_dim != idl_skl_rev->value.arr->n_dim)
		_p3d_idlPrintNamedError("IMAGE argument and MASK argument must be 3D matrices with identical dimensions.");

	if (idl_in_rev->value.arr->dim[0] != idl_skl_rev->value.arr->dim[0])
		_p3d_idlPrintNamedError("IMAGE argument and MASK argument must be 3D matrices with identical dimensions.");

	if (idl_in_rev->value.arr->dim[1] != idl_skl_rev->value.arr->dim[1])
		_p3d_idlPrintNamedError("IMAGE argument and MASK argument must be 3D matrices with identical dimensions.");

	if (idl_in_rev->value.arr->dim[2] != idl_skl_rev->value.arr->dim[2])
		_p3d_idlPrintNamedError("IMAGE argument and MASK argument must be 3D matrices with identical dimensions.");


	if ( idl_in_rev->value.arr->n_dim == 3 )
	{
		// Extract input in C format
		if (idl_in_rev->type == IDL_TYP_BYTE)			
			in_rev = (unsigned char *) idl_in_rev->value.arr->data;  
		else 
			_p3d_idlPrintNamedError("Input argument VOLUME must be of type BYTE.");	

		// Extract input in C format
		if (idl_skl_rev->type == IDL_TYP_BYTE)			
			skl_rev = (unsigned char *) idl_skl_rev->value.arr->data;  
		else 
			_p3d_idlPrintNamedError("Input argument VOLUME must be of type BYTE.");	
	}
	else
	{
		_p3d_idlPrintNamedError("Input argument VOLUME must be a 3D matrix.");
	}

	ratio = (double *) IDL_MakeTempVector(	IDL_TYP_DOUBLE, 1,  
						IDL_ARR_INI_NOP, &idl_out_ratio	); 
		
	// Call Pore3D:
	err_code = p3dSkeletonAnalysisFeasibility ( 
		in_rev,
		skl_rev,
		ratio,
		(int) idl_in_rev->value.arr->dim[0],
		(int) idl_in_rev->value.arr->dim[1],
		(int) idl_in_rev->value.arr->dim[2],
		_p3d_idlPrintInfo
	);	

	if (err_code == P3D_ERROR)
	{
		// Print error:
		_p3d_idlPrintNamedError("Error on internal code execution.");
	}				

	// Return output in IDL Format:	
	return idl_out_ratio;  
} 
