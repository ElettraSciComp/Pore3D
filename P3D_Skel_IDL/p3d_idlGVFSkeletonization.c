// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dSkel.h"


IDL_VPTR p3d_idlGVFSkeletonization (int argc, IDL_VPTR argv[], char* argk)  
{  
	typedef struct {  
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		double eps;
		int eps_there;
		double hierarc;
		int hierarc_there;
		double mu;
		int mu_there;
		double scale;
		int scale_there;		
	} KW_RESULT;  

   	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {  
		IDL_KW_FAST_SCAN,  
		{ "EPS", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(eps_there), (char*) IDL_KW_OFFSETOF(eps) },  		
		{ "HIERARC", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(hierarc_there), (char*) IDL_KW_OFFSETOF(hierarc) },  		
		{ "MU", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(mu_there), (char*) IDL_KW_OFFSETOF(mu) },  		
		{ "SCALE", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(scale_there), (char*) IDL_KW_OFFSETOF(scale) }, 
		{ NULL }  
	 };  

	KW_RESULT kw;
	
	IDL_VPTR idl_out_im, idl_in_im;
	unsigned char *in_im8, *out_im8;  		
	int keywords_ct = 0;

	double mu  = 0.15;
	double eps = 1E-4;
	double scale = 1.0;
	double hierarc = 0.0;
		
	int err_code;
	
	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 

	// Get input data in IDL format:
	idl_in_im = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_im);  
	IDL_ENSURE_ARRAY(idl_in_im);  


	// Get the MU keyword:
	if (kw.mu_there)
	{
		// Check values:
		if ( kw.mu <= 0 )
			_p3d_idlPrintNamedError("MU must be greater than 0.");

		// Get values:
		mu = (double) kw.mu;

		keywords_ct++;
	}
	// Get the EPS keyword:
	if (kw.eps_there)
	{
		// Check values:
		if ( kw.eps <= 0 )
			_p3d_idlPrintNamedError("EPS must be greater than 0.");

		// Get values:
		eps = (double) kw.eps;

		keywords_ct++;
	}
	// Get the SCALE keyword:
	if (kw.scale_there)
	{
		// Check values:
		if ( kw.scale <= 0 )
			_p3d_idlPrintNamedError("SCALE must be greater than 0.");

		// Get values:
		scale = (double) kw.scale;

		keywords_ct++;
	}
	// Get the HIERARC keyword:
	if (kw.hierarc_there)
	{
		// Check values:
		if ( kw.hierarc <= 0 )
			_p3d_idlPrintNamedError("HIERARC must be greater than 0.");

		// Get values:
		hierarc = (double) kw.hierarc;

		keywords_ct++;
	}

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
			err_code = p3dGVFSkeletonization ( 
				in_im8,
				out_im8,
				idl_in_im->value.arr->dim[0],
				idl_in_im->value.arr->dim[1],
				idl_in_im->value.arr->dim[2],
				mu,
				eps,
				hierarc,
				scale,
				_p3d_idlPrintInfo
			);

			// On exception print error:
			if (err_code == P3D_MEM_ERROR)
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

