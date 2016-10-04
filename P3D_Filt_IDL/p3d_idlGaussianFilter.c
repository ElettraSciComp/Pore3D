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

IDL_VPTR p3d_idlGaussianFilter(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        double sigma;
        int si_there;
        IDL_LONG nsize;
        int ns_there;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "SIGMA", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(si_there), (char*) IDL_KW_OFFSETOF(sigma)},
        { "WIDTH", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(ns_there), (char*) IDL_KW_OFFSETOF(nsize)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev8, *out_rev8;
    unsigned short *in_rev16, *out_rev16;
    int keywords_ct = 0;

    int    width = 3;   // default
    double sigma = 1.0; // default

    int err_code;

    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    // Get the WIDTH input argument:
    if (kw.ns_there) {
        // Check values:
        if ((kw.nsize < 3) || (kw.nsize > 51))
            _p3d_idlPrintNamedError("WIDTH must be an odd value within the range [3,51].");

        if ((kw.nsize % 2) == 0)
            _p3d_idlPrintNamedError("WIDTH must be an odd value within the range [3,51].");

        // Get values:
        width = (int) kw.nsize;

        keywords_ct++;
    }

    // Get the SIGMA input argument:
    if (kw.si_there) {
        // Check values:
        if (kw.sigma < 0)
            _p3d_idlPrintNamedError("SIGMA must be greater than zero.");

        // Get values:
        sigma = (double) kw.sigma;

        keywords_ct++;
    }


    // Call Pore3D depending on input arguments:
    /*if (idl_in_rev->value.arr->n_dim == 2) {
        // Extract first input (volume to filter) in C format:
        if (idl_in_rev->type == IDL_TYP_BYTE) {
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
            err_code = p3dGaussianFilter2D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    width,
                    sigma,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");

        } else if (idl_in_rev->type == IDL_TYP_UINT) {
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
            err_code = p3dGaussianFilter2D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    width,
                    sigma,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");

        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } else*/ if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract first input (volume to filter) in C format:
        if (idl_in_rev->type == IDL_TYP_BYTE) {
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
            err_code = p3dGaussianFilter3D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    width,
                    sigma,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");

        } else if (idl_in_rev->type == IDL_TYP_UINT) {
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
            err_code = p3dGaussianFilter3D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    width,
                    sigma,
                    _p3d_idlPrintInfo,
                    NULL
                    );

            // On exception print error:
            if ((err_code == P3D_IO_ERROR) || (err_code == P3D_MEM_ERROR))
                _p3d_idlPrintNamedError("Error on code execution.");

        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
    }


    // Free resources:
    IDL_KW_FREE;

    // Return output in IDL Format
    return (idl_out_rev);
}

