
// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dFilt.h"


IDL_VPTR p3d_idlBilateralFilter(int argc, IDL_VPTR argv[], char* argk)
{ 
	typedef struct {  
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		IDL_LONG nsize;  
		int ns_there;  
		IDL_LONG iter;  
		int it_there;  
		double sigma_d;
		int sd_there;
		double sigma_r;  
		int sr_there; 		
	} KW_RESULT;  

	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {  
		IDL_KW_FAST_SCAN,  
		{ "ITERATIONS", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(it_there), (char*) IDL_KW_OFFSETOF(iter) },  		
		{ "SIGMA_D", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(sd_there), (char*) IDL_KW_OFFSETOF(sigma_d) },  		
		{ "SIGMA_R", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(sr_there), (char*) IDL_KW_OFFSETOF(sigma_r) },  		
		{ "WIDTH", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(ns_there), (char*) IDL_KW_OFFSETOF(nsize) },  		
		{ NULL }  
	};  

	KW_RESULT kw;

	IDL_VPTR idl_out_rev, idl_in_rev;
	unsigned char *in_rev8, *out_rev8;  	
	unsigned short *in_rev16, *out_rev16;
	int keywords_ct = 0;

	int iter = 10;
	int size = 3; // default
	double sigma_d = 1.0;  // default
	double sigma_r = 3.0;  // default


	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 


	// Get input data in IDL format:
	idl_in_rev = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_rev);  
	IDL_ENSURE_ARRAY(idl_in_rev);  


	// Get the SIZE input argument:
	if (kw.ns_there)
	{
		// Check values:
		if ( ( kw.nsize < 3 ) || ( kw.nsize > 51 ) )
			_p3d_idlPrintNamedError("WIDTH must be an odd integer value within the range [3,51].");

		if ( ( kw.nsize % 2) == 0 )
			_p3d_idlPrintNamedError("WIDTH must be an odd integer value within the range [3,51].");

		// Get values:
		size = (int) kw.nsize;

		keywords_ct++;
	}

	// Get the SIZE input argument:
	if (kw.it_there)
	{
		// Check values:
		if ( kw.iter <= 0 ) 
			_p3d_idlPrintNamedError("ITERATIONS must be greater than 0.");

		// Get values:
		iter = (int) kw.iter;

		keywords_ct++;
	}


	// Get the SIGMA_D input argument:
	if (kw.sd_there)
	{
		// Check values:
		if ( kw.sigma_d <= 0 )
			_p3d_idlPrintNamedError("SIGMA_D must be greater than 0.");

		// Get values:
		sigma_d = (double) kw.sigma_d;

		keywords_ct++;
	}

	// Get the SIGMA_R input argument:
	if (kw.sr_there)
	{
		// Check values:
		if ( kw.sigma_r <= 0 )
			_p3d_idlPrintNamedError("SIGMA_R must be greater than 0.");

		// Get values:
		sigma_r = (double) kw.sigma_r;

		keywords_ct++;
	}


	// Call Pore3D depending on input arguments:
	if ( idl_in_rev->value.arr->n_dim == 3 )
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
			p3dBilateralFilter3D_8 ( 
				in_rev8,
				out_rev8,
				(int) idl_in_rev->value.arr->dim[0],
				(int) idl_in_rev->value.arr->dim[1],
				(int) idl_in_rev->value.arr->dim[2],
				size,
				sigma_d,
				sigma_r,
				iter,
				_p3d_idlPrintInfo,
				NULL
				);
		}
		else if (idl_in_rev->type == IDL_TYP_UINT)		
		{
			in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;  

			// Allocate memory for output:
			if (!(idl_in_rev->flags & IDL_V_TEMP))  
				out_rev16 = (unsigned short *) IDL_MakeTempArray(
				IDL_TYP_UINT,
				idl_in_rev->value.arr->n_dim,  
				idl_in_rev->value.arr->dim,  
				IDL_ARR_INI_NOP, 
				&idl_out_rev
				);   

			// Call Pore3D:
			p3dBilateralFilter3D_16 ( 
				in_rev16,
				out_rev16,
				(int) idl_in_rev->value.arr->dim[0],
				(int) idl_in_rev->value.arr->dim[1],
				(int) idl_in_rev->value.arr->dim[2],
				size,
				sigma_d,
				sigma_r,
				iter,
				_p3d_idlPrintInfo,
				NULL
				);	
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