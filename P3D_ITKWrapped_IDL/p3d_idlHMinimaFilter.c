
// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dITKWrapped.h"


IDL_VPTR p3d_idlHMinimaFilter(int argc, IDL_VPTR argv[], char* argk)
{ 
	typedef struct {  
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		IDL_LONG conn;
		int cn_there;  		
	} KW_RESULT;  
  
	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {  
		IDL_KW_FAST_SCAN,  
		{ "CONN", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(cn_there), (char*) IDL_KW_OFFSETOF(conn) },  		
    	{ NULL }  
	 };  

	KW_RESULT kw;

	IDL_VPTR idl_out_rev, idl_in_rev;
	unsigned char *in_rev8, *out_rev8;  	
	unsigned short *in_rev16, *out_rev16;

	int err_code;
	int thresh;

	int keywords_ct = 0;

	int conn2D = CONN8;
	int conn3D = CONN26;

	
	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 

		
	// Get input data in IDL format:
	idl_in_rev = argv[0]; 

	IDL_ENSURE_SIMPLE(idl_in_rev);  
	IDL_ENSURE_ARRAY(idl_in_rev);  

	if ( argv[1]->type == IDL_TYP_BYTE ) 
		thresh = (int) argv[1]->value.c;	
	else if ( argv[1]->type == IDL_TYP_INT ) 
		thresh = (int) argv[1]->value.i;	
	else if ( argv[1]->type == IDL_TYP_UINT ) 
		thresh = (int) argv[1]->value.ui;
	else if ( argv[1]->type == IDL_TYP_LONG )
		thresh = (int) argv[1]->value.l;	
	else if ( argv[1]->type == IDL_TYP_ULONG )
		thresh = (int) argv[1]->value.ul;	
	else if ( argv[1]->type == IDL_TYP_LONG64 )
		thresh = (int) argv[1]->value.l64;	
	else if ( argv[1]->type == IDL_TYP_ULONG64 )
		thresh = (int) argv[1]->value.ul64;				
	else if ( argv[1]->type == IDL_TYP_FLOAT )
		thresh = (int) argv[1]->value.f;			
	else if ( argv[1]->type == IDL_TYP_DOUBLE )
		thresh = (int) argv[1]->value.d;

	// Get the CONN input argument:
	if (kw.cn_there)
	{
		if ( idl_in_rev->value.arr->n_dim == 2 )
		{
			// Check values:
			if ( ( kw.conn != 4 ) && ( kw.conn != 8) )
				_p3d_idlPrintNamedError("CONN must be a value of the set {4,8}.");
			
			// Get values:
			if ( kw.conn == 4 )
				conn2D = CONN4;
			else if ( kw.conn == 8 )
				conn2D = CONN8;
		}
		else if ( idl_in_rev->value.arr->n_dim == 3 )
		{
			// Check values:
			if ( ( kw.conn != 6 ) && ( kw.conn != 18 ) && ( kw.conn != 26 ) )
				_p3d_idlPrintNamedError("CONN must be a value of the set {6,18,26}.");
			
			// Get values:
			if ( kw.conn == 6 )
				conn3D = CONN6;
			else if ( kw.conn == 18 )
				conn3D = CONN18;
			else if ( kw.conn == 26 )
				conn3D = CONN26;
		}
		// else: error on input arguments with further handling.

		keywords_ct++;
	}

	
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
 

			// Call Pore3D without the mask:
			err_code = p3dHMinimaFilter2D_8 ( 
				in_rev8,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				conn2D,
				(unsigned char) thresh,
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
				out_rev16 = (unsigned short *) IDL_MakeTempArray(
					IDL_TYP_UINT,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					);   
			
			// Call Pore3D without the mask:
			err_code = p3dHMinimaFilter2D_16 ( 
				in_rev16,
				out_rev16,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				conn2D,
				(unsigned short) thresh,
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
			err_code = p3dHMinimaFilter3D_8 ( 
				in_rev8,
				out_rev8,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				conn3D,
				(unsigned char) thresh,
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
				out_rev16 = (unsigned short *) IDL_MakeTempArray(
					IDL_TYP_UINT,
					idl_in_rev->value.arr->n_dim,  
					idl_in_rev->value.arr->dim,  
					IDL_ARR_INI_NOP, 
					&idl_out_rev
					); 
 
			
			// Call Pore3D without the mask:
			err_code = p3dHMinimaFilter3D_16 ( 
				in_rev16,
				out_rev16,
				(unsigned int) idl_in_rev->value.arr->dim[0],
				(unsigned int) idl_in_rev->value.arr->dim[1],
				(unsigned int) idl_in_rev->value.arr->dim[2],
				conn3D,
				(unsigned short) thresh,
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

