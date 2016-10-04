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

#include "p3dBlob.h"

IDL_VPTR p3d_idlChamferDT(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG w_data[3];
        int w_there;
        IDL_MEMINT w_n;
    } KW_RESULT;

    static IDL_KW_ARR_DESC_R weights = {(char*) IDL_KW_OFFSETOF(w_data), 3, 3, (IDL_LONG*) IDL_KW_OFFSETOF(w_n)};


    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "WEIGHTS", IDL_TYP_LONG, 1, IDL_KW_ARRAY, (int*) IDL_KW_OFFSETOF(w_there), (char*) (&weights)},
        { NULL}
    };

    KW_RESULT kw;


    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev;
    unsigned short *out_rev;
    int keywords_ct = 0;

    unsigned int w1, w2, w3;

    int err_code;


    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get WEIGHTS input keyword:
    if (kw.w_there) {
        // Check values:
        if (kw.w_n == 3) {
            if (kw.w_data[0] < 0)
                _p3d_idlPrintNamedError("W1 value of input keyword WEIGHTS must be an integer value greater than 0.");

            if (kw.w_data[1] < 0)
                _p3d_idlPrintNamedError("W2 value of input keyword WEIGHTS must be an integer value greater than 0.");

            if (kw.w_data[2] < 0)
                _p3d_idlPrintNamedError("W3 value of input keyword WEIGHTS must be an integer value greater than 0.");

            // Get values:
            w1 = (unsigned int) kw.w_data[0];
            w2 = (unsigned int) kw.w_data[1];
            w3 = (unsigned int) kw.w_data[2];
        }
        else {
            _p3d_idlPrintNamedError("Input argument WEIGHTS must contain three [ W1, W2, W3 ] elements.");
        }

        keywords_ct++;
    } else {
        // Set default values for w1, w2 and w3:
        w1 = 3;
        w2 = 4;
        w3 = 5;
    }


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev); // Get input data in IDL format:

    // Allocate memory for output:
    if (!(idl_in_rev->flags & IDL_V_TEMP))
        out_rev = (unsigned short *) IDL_MakeTempArray(
            IDL_TYP_UINT,
            idl_in_rev->value.arr->n_dim,
            idl_in_rev->value.arr->dim,
            IDL_ARR_INI_NOP,
            &idl_out_rev
            );

    if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract input in C format
        if (idl_in_rev->type == IDL_TYP_BYTE)
            in_rev = (unsigned char *) idl_in_rev->value.arr->data;
        else
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");

        // Call Pore3D:
        err_code = p3dChamferDT(
                in_rev,
                out_rev,
                (unsigned int) idl_in_rev->value.arr->dim[0],
                (unsigned int) idl_in_rev->value.arr->dim[1],
                (unsigned int) idl_in_rev->value.arr->dim[2],
                w1,
                w2,
                w3,
                _p3d_idlPrintInfo
                );

        if (err_code == P3D_MEM_ERROR) {
            // Print error:
            _p3d_idlPrintNamedError("Error on internal code execution.");
        }

        //
        // Following code is not necessary 'cause number of input arguments
        // is checked at compile-time and not at run time (see further IDL_Load).
        //
        //else
        //{
        //	// Print error:
        //	_p3d_idlPrintNamedError("Incorrect number of arguments.");
        //}
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
    }

    // Release keyword resources:
    IDL_KW_FREE;

    // Return output in IDL Format:
    return idl_out_rev;
}
