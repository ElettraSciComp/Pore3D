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

#include "p3dFilt.h"

IDL_VPTR p3d_idlFrom16To8(int argc, IDL_VPTR argv[], char* argk) {

    IDL_VPTR idl_out_rev, idl_in_rev, idl_dims;
    unsigned char *in_rev8, *out_rev8;
    unsigned short *in_rev16, *out_rev16;
    IDL_MEMINT tmp_dims[2];
    IDL_INT* dims;
    
    unsigned short min, max;
    int err_code;

    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);

    // Get DIMX and DIMY input arguments:		
    idl_dims = argv[1];

    IDL_ENSURE_SIMPLE(idl_dims);
    IDL_ENSURE_ARRAY(idl_dims);


    if ((idl_dims->type != IDL_TYP_BYTE) && (idl_dims->type != IDL_TYP_INT) &&
            (idl_dims->type != IDL_TYP_UINT) && (idl_dims->type != IDL_TYP_LONG) &&
            (idl_dims->type != IDL_TYP_ULONG) && (idl_dims->type != IDL_TYP_LONG64) &&
            (idl_dims->type != IDL_TYP_ULONG64))
        _p3d_idlPrintNamedError("Input argument RANGE must be an array of integer type.");

        
    // User wants to read 2D RAW data:
    if (idl_dims->value.arr->n_elts == 2) {


        if (idl_dims->type == IDL_TYP_BYTE)
            dims = (UCHAR*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_INT)
            dims = (IDL_INT*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_UINT)
            dims = (IDL_UINT*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_LONG)
            dims = (IDL_LONG*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_ULONG)
            dims = (IDL_ULONG*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_LONG64)
            dims = (IDL_LONG64*) idl_dims->value.arr->data;
        else if (idl_dims->type == IDL_TYP_ULONG64)
            dims = (IDL_ULONG64*) idl_dims->value.arr->data;

        tmp_dims[0] = dims[0];
        tmp_dims[1] = dims[1];
    } else {
        _p3d_idlPrintNamedError("Input argument RANGE must contain two [ MIN, MAX ] elements.");
    }

    if (idl_in_rev->value.arr->n_dim == 3) {
        if (idl_in_rev->type == IDL_TYP_UINT) {
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
            err_code = p3dFrom16To8(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    (unsigned short) tmp_dims[0],
                    (unsigned short) tmp_dims[1],
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");

        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type UINT.");
        }
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
    }


    // Return output in IDL Format
    return (idl_out_rev);
}
