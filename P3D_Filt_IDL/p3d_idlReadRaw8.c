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

IDL_VPTR p3d_idlReadRaw8(int argc, IDL_VPTR argv[], char* argk) {
    IDL_VPTR idl_out_im;
    IDL_VPTR idl_dims;

    char* filename;
    int err_code;

    unsigned char *out_im;


    //
    // Extract inputs in C format:
    //
    if (argv[0]->type == IDL_TYP_STRING)
        filename = IDL_VarGetString(argv[0]);
    else
        _p3d_idlPrintNamedError("Input argument FILENAME must be a string.");

    // Get input data in IDL format:
    idl_dims = argv[1];

    IDL_ENSURE_SIMPLE(idl_dims);
    IDL_ENSURE_ARRAY(idl_dims);

    if (idl_dims->type != IDL_TYP_INT)
        _p3d_idlPrintNamedError("Input argument DIMS must be an array of INT.");

    if ((idl_dims->value.arr->n_elts < 2) || (idl_dims->value.arr->n_elts > 3))
        _p3d_idlPrintNamedError("Input argument DIMS must contain two [ X, Y ] or three [ X, Y, Z ] elements.");


    // User wants to read 2D RAW data:
    if (idl_dims->value.arr->n_elts == 2) {
        // Get input data
        IDL_INT* dims = (IDL_INT*) idl_dims->value.arr->data;
        IDL_MEMINT tmp_dims[2];

        tmp_dims[0] = dims[0];
        tmp_dims[1] = dims[1];

        //
        // Allocate memory for output:
        //
        out_im = (unsigned char *) IDL_MakeTempArray(
                IDL_TYP_BYTE,
                2,
                tmp_dims,
                IDL_ARR_INI_NOP,
                &idl_out_im
                );

        //
        // Call Pore3D for 2D reading:
        //
        err_code = p3dReadRaw8(
                filename,
                out_im,
                tmp_dims[0],
                tmp_dims[1],
                1,
                _p3d_idlPrintInfo,
                NULL
                );

        // On exception print error:
        if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
            _p3d_idlPrintNamedError("Error on code execution.");
    }        // User wants to read 3D data:
    else {
        // Get input data
        IDL_INT* dims = (IDL_INT*) idl_dims->value.arr->data;
        IDL_MEMINT tmp_dims[3];

        tmp_dims[0] = dims[0];
        tmp_dims[1] = dims[1];
        tmp_dims[2] = dims[2];


        //
        // Allocate memory for output:
        //
        out_im = (unsigned char *) IDL_MakeTempArray(
                IDL_TYP_BYTE,
                3,
                tmp_dims,
                IDL_ARR_INI_NOP,
                &idl_out_im
                );

        //
        // Call Pore3D for 3D reading:
        //
        err_code = p3dReadRaw8(
                filename,
                out_im,
                tmp_dims[0],
                tmp_dims[1],
                tmp_dims[2],
                _p3d_idlPrintInfo,
                NULL
                );

        // On exception print error:
        if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
            _p3d_idlPrintNamedError("Error on code execution.");
    }


    // Return output in IDL format:
    return idl_out_im;
}
