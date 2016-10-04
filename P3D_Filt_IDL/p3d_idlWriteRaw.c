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

void p3d_idlWriteRaw(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG endian;
        IDL_LONG sign;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "BIG_ENDIAN", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(endian)},
        { "UNSIGNED", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(sign)},
        { NULL}
    };

    KW_RESULT kw;
    IDL_VPTR idl_in_rev;

    char* filename;
    int err_code;

    int little_endian;
    int is_signed;

    unsigned char* in_rev8;
    unsigned short* in_rev16;
    unsigned int* in_rev32;


    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, (IDL_VPTR *) 0, 1, &kw);

    little_endian = (kw.endian == 0) ? P3D_TRUE : P3D_FALSE;
    is_signed = (kw.sign == 0) ? P3D_TRUE : P3D_FALSE;

    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    if (argv[1]->type == IDL_TYP_STRING)
        filename = IDL_VarGetString(argv[1]);
    else
        _p3d_idlPrintNamedError("Input argument FILENAME must be a string.");


    // Check if user wants to write a 8-bit or 16-bit format image:
    if (idl_in_rev->type == IDL_TYP_BYTE) {
        // Extract input in C format
        in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;


        // Check if user wants to write a 2D image or a 3D volume:
        if (idl_in_rev->value.arr->n_dim == 2) {
            // Call Pore3D:
            err_code = p3dWriteRaw8(
                    in_rev8,
                    filename,
                    idl_in_rev->value.arr->dim[0],
                    idl_in_rev->value.arr->dim[1],
                    1,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");
        } else if (idl_in_rev->value.arr->n_dim == 3) {


            // Call Pore3D:
            err_code = p3dWriteRaw8(
                    in_rev8,
                    filename,
                    idl_in_rev->value.arr->dim[0],
                    idl_in_rev->value.arr->dim[1],
                    idl_in_rev->value.arr->dim[2],
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");
        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
        }
    } else if (idl_in_rev->type == IDL_TYP_UINT) {
        // Extract input in C format
        in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;

        // Check if user wants to write a 2D image or a 3D volume:
        if (idl_in_rev->value.arr->n_dim == 2) {


            // Call Pore3D:
            err_code = p3dWriteRaw16(
                    in_rev16,
                    filename,
                    idl_in_rev->value.arr->dim[0],
                    idl_in_rev->value.arr->dim[1],
                    1,
                    little_endian,
                    is_signed,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");
        } else if (idl_in_rev->value.arr->n_dim == 3) {


            // Call Pore3D:
            err_code = p3dWriteRaw16(
                    in_rev16,
                    filename,
                    idl_in_rev->value.arr->dim[0],
                    idl_in_rev->value.arr->dim[1],
                    idl_in_rev->value.arr->dim[2],
                    little_endian,
                    is_signed,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");
             
        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix."); 
        }
    } else if (idl_in_rev->type == IDL_TYP_ULONG) {
            // Extract input in C format
            in_rev32 = (unsigned int *) idl_in_rev->value.arr->data;

            // Check if user wants to write a 2D image or a 3D volume:
            if (idl_in_rev->value.arr->n_dim == 2) {


                // Call Pore3D:
                err_code = p3dWriteRaw32(
                        in_rev32,
                        filename,
                        idl_in_rev->value.arr->dim[0],
                        idl_in_rev->value.arr->dim[1],
                        1,
                        little_endian,
                        is_signed,
                        _p3d_idlPrintInfo,
                        NULL
                        );

                // On exception print error:
                if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                    _p3d_idlPrintNamedError("Error on code execution.");
            } else if (idl_in_rev->value.arr->n_dim == 3) {


                // Call Pore3D:
                err_code = p3dWriteRaw32(
                        in_rev32,
                        filename,
                        idl_in_rev->value.arr->dim[0],
                        idl_in_rev->value.arr->dim[1],
                        idl_in_rev->value.arr->dim[2],
                        little_endian,
                        is_signed,
                        _p3d_idlPrintInfo,
                        NULL
                        );

                // On exception print error:
                if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                    _p3d_idlPrintNamedError("Error on code execution.");
            }
            else {
                _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
            }
        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be a BYTE or UINT matrix.");
        }


        // Free keywords resources:
        IDL_KW_FREE;
    }
