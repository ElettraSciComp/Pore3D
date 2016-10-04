/***************************************************************************/
/* (C) 2016 Elettra - Sincrotrone Trieste S.C.p.A.. All rights reserved.   */
/*                                                                         */
/*                                                                         */
/* This file is part of Pore3D, a software library for quantitative        */
/* analysis of 3D (volume) images.                                         */
/*                                                                         */
/* Pore3D is free software: you can redistribute it and/or modify it       */
/* under the terms of the GNU General Public License as published by the   */
/* Free Software Foundation, either version 3 of the License, or (at your  */
/* option) any later version.                                              */
/*                                                                         */
/* Pore3D is distributed in the hope that it will be useful, but WITHOUT   */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   */
/* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    */
/* for more details.                                                       */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with Pore3D. If not, see <http://www.gnu.org/licenses/>.          */
/*                                                                         */
/***************************************************************************/

//
// Author: Francesco Brun
// Last modified: Sept, 28th 2016
//

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
				idl_in_im->value.arr->dim[0],
				idl_in_im->value.arr->dim[1],
				idl_in_im->value.arr->dim[2],				
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

	// Return output in IDL format:
	return idl_out_im;  
} 

