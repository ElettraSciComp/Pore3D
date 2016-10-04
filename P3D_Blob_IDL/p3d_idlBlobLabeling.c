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

IDL_VPTR p3d_idlBlobLabeling(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure        
        IDL_LONG rand;
        IDL_LONG skip_borders;
        IDL_LONG uint_output;
        IDL_LONG conn;
        int cn_there;
    } KW_RESULT;


    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "CONN", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(cn_there), (char*) IDL_KW_OFFSETOF(conn)},
        { "RAND", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(rand)},
        { "SKIPBORDERS", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(skip_borders)},
        { "ULONG", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(uint_output)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev8;
    unsigned short *out_rev16;
    unsigned int *out_rev32;

    int keywords_ct = 0;

    int rand;
    int skip_borders;
    int uint_output;
    int conn2D = CONN8;
    int conn3D = CONN6;

    int err_code;

    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);

    rand = (kw.rand == 0) ? P3D_FALSE : P3D_TRUE;
    skip_borders = (kw.skip_borders == 0) ? P3D_FALSE : P3D_TRUE;
    uint_output = (kw.uint_output == 0) ? P3D_FALSE : P3D_TRUE;


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    // Get the CONN input argument:
    if (kw.cn_there) {
        if (idl_in_rev->value.arr->n_dim == 2) {
            // Check values:
            if ((kw.conn != 4) && (kw.conn != 8))
                _p3d_idlPrintNamedError("CONN must be a value of the set {4,8}.");

            // Get values:
            if (kw.conn == 4)
                conn2D = CONN4;
            else if (kw.conn == 8)
                conn2D = CONN8;
        } else if (idl_in_rev->value.arr->n_dim == 3) {
            // Check values:
            if ((kw.conn != 6) && (kw.conn != 18) && (kw.conn != 26))
                _p3d_idlPrintNamedError("CONN must be a value of the set {6,18,26}.");

            // Get values:
            if (kw.conn == 6)
                conn3D = CONN6;
            else if (kw.conn == 18)
                conn3D = CONN18;
            else if (kw.conn == 26)
                conn3D = CONN26;
        }
        // else: error on input arguments with further handling.

        keywords_ct++;
    }

    // Call Pore3D depending on input arguments:
    /*if ( idl_in_rev->value.arr->n_dim == 2 )
    {
            // Extract first input (volume to filter) in C format:
            if (idl_in_rev->type == IDL_TYP_BYTE)
            {
                    in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;

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
                    err_code = p3dBlobLabeling2D (
                            in_rev8,
                            out_rev16,
                            NULL,
                            &num_cc,
                            (int) idl_in_rev->value.arr->dim[0],
                            (int) idl_in_rev->value.arr->dim[1],
                            conn2D,
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
else */if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract first input (volume to filter) in C format:
        if (idl_in_rev->type == IDL_TYP_BYTE) {
            in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;

            if (uint_output == P3D_FALSE) {
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
                err_code = p3dBlobLabeling_ushort(
                        in_rev8,
                        out_rev16,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        conn3D,
                        rand,
                        skip_borders,
                        _p3d_idlPrintInfo
                        );
            } else {
                // Allocate memory for output:
                if (!(idl_in_rev->flags & IDL_V_TEMP))
                    out_rev32 = (unsigned int *) IDL_MakeTempArray(
                        IDL_TYP_ULONG,
                        idl_in_rev->value.arr->n_dim,
                        idl_in_rev->value.arr->dim,
                        IDL_ARR_INI_NOP,
                        &idl_out_rev
                        );


                // Call Pore3D:
                err_code = p3dBlobLabeling_uint(
                        in_rev8,
                        out_rev32,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        conn3D,
                        rand,
                        skip_borders,
                        _p3d_idlPrintInfo
                        );
            }

            // On exception print error:
            if (err_code == P3D_MEM_ERROR)
                _p3d_idlPrintNamedError("Error on code execution.");
        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");
        }
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
    }


    // Free resources:
    IDL_KW_FREE;

    // Return output in IDL Format
    return (idl_out_rev);
}


