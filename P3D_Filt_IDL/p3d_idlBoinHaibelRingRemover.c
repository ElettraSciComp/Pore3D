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
#include <limits.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dFilt.h"

IDL_VPTR p3d_idlBoinHaibelRingRemover(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG winsize;
        int ws_there;
        IDL_LONG iterations;
        int it_there;
        double precision;
        int pr_there;
        double cen_data[2];
        int cen_there;
        IDL_MEMINT cen_n;
    } KW_RESULT;

    static IDL_KW_ARR_DESC_R center = {(char*) IDL_KW_OFFSETOF(cen_data), 2, 2, (IDL_LONG*) IDL_KW_OFFSETOF(cen_n)};

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "CENTER", IDL_TYP_DOUBLE, 1, IDL_KW_ARRAY, (int*) IDL_KW_OFFSETOF(cen_there), (char*) (&center)},
        { "ITERATIONS", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(it_there), (char*) IDL_KW_OFFSETOF(iterations)},
        { "PRECISION", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(pr_there), (char*) IDL_KW_OFFSETOF(precision)},
        { "WINSIZE", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(ws_there), (char*) IDL_KW_OFFSETOF(winsize)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev, idl_mask;
    unsigned char *in_rev8, *out_rev8, *mask;
    unsigned short *in_rev16, *out_rev16;

    int iterations = 1; // default
    double precision = 1.5; // default
    int winsize = 51; // default
    double centerX, centerY;

    int keywords_ct = 0;



    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);

    // Extract function specific parameters:

    // Get the ITERATIONS input argument:
    if (kw.it_there) {
        // Check values:
        if (kw.iterations < 0)
            _p3d_idlPrintNamedError("ITERATIONS must be an integer value greater than 0.");

        // Get values:
        iterations = (unsigned int) kw.iterations;

        keywords_ct++;
    }

    // Get the PRECISION input argument:
    if (kw.pr_there) {
        // Check values:
        if ((kw.precision < 1.0) && (kw.precision > 5.0))
            _p3d_idlPrintNamedError("PRECISION must be a decimal value in the range [1.0,5.0].");

        // Get values:
        precision = (double) kw.precision;

        keywords_ct++;
    }

    // Get the WINSIZE input argument:
    if (kw.ws_there) {
        // Check values:
        if ((kw.winsize < 3) || (kw.winsize > 201))
            _p3d_idlPrintNamedError("WINSIZE must be an odd value within the range [3,201].");

        if ((kw.winsize % 2) == 0)
            _p3d_idlPrintNamedError("WINSIZE must be an odd value within the range [3,201].");

        // Get values:
        winsize = (unsigned int) kw.winsize;

        keywords_ct++;
    }


    // Get CENTER input argument:
    if (kw.cen_there) {
        // Check values:
        if (kw.cen_n != 2)
            _p3d_idlPrintNamedError("Input argument CENTER must contain two [ X, Y ] elements.");

        if ((kw.cen_data[0] < 0) || (kw.cen_data[0] > idl_in_rev->value.arr->dim[0]))
            _p3d_idlPrintNamedError("X value of input argument CENTER must be within IMAGE dimensions.");

        if ((kw.cen_data[1] < 0) || (kw.cen_data[1] > idl_in_rev->value.arr->dim[1]))
            _p3d_idlPrintNamedError("Y value of input argument CENTER must be within IMAGE dimensions.");

        // Get values:
        centerX = (double) kw.cen_data[0];
        centerY = (double) kw.cen_data[1];

        keywords_ct++;
    } else {
        // Set default values for centerX and centerY:
        centerX = (double) idl_in_rev->value.arr->dim[0] / 2.0;
        centerY = (double) idl_in_rev->value.arr->dim[1] / 2.0;
    }


    // Call Pore3D depending on input arguments:
    if (idl_in_rev->value.arr->n_dim == 2) {
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
            p3dBoinHaibelRingRemover2D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    centerX,
                    centerY,
                    winsize,
                    iterations,
                    precision,
                    _p3d_idlPrintInfo
                    );

            //
            // Following code is not necessary 'cause number of input arguments
            // is checked at compile-time and not at run time (see further IDL_Load).
            //
            //else
            //{
            //	// Print error:
            //	_p3d_idlPrintNamedError("Incorrect number of arguments.");
            //}
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
            p3dBoinHaibelRingRemover2D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    centerX,
                    centerY,
                    winsize,
                    iterations,
                    precision,
                    _p3d_idlPrintInfo
                    );

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
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    }/*else if (idl_in_rev->value.arr->n_dim == 3) {
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


            // Check if MASK is present:
            if ((argc - keywords_ct) == 2) {
                // Call Pore3D without the mask:
                p3dSijbersPostnovRingRemover3D_8(
                        in_rev8,
                        out_rev8,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        centerX,
                        centerY,
                        winsize,
                        thresh,
                        iterations,
                        precision,
                        NULL,
                        _p3d_idlPrintInfo
                        );
            } else if ((argc - keywords_ct) == 3) {
                // Get input data in IDL format:
                idl_mask = argv[2];

                IDL_ENSURE_SIMPLE(idl_mask);
                IDL_ENSURE_ARRAY(idl_mask);

                // Extract input in C format
                mask = (unsigned char *) idl_mask->value.arr->data;


                // Call Pore3D using the mask:
                p3dSijbersPostnovRingRemover3D_8(
                        in_rev8,
                        out_rev8,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        centerX,
                        centerY,
                        winsize,
                        thresh,
                        iterations,
                        precision,
                        mask,
                        _p3d_idlPrintInfo
                        );
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


            // Check if MASK is present:
            if ((argc - keywords_ct) == 2) {
                // Call Pore3D without the mask:
                p3dSijbersPostnovRingRemover3D_16(
                        in_rev16,
                        out_rev16,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        centerX,
                        centerY,
                        winsize,
                        thresh,
                        iterations,
                        precision,
                        NULL,
                        _p3d_idlPrintInfo
                        );
            } else if ((argc - keywords_ct) == 3) {
                // Get input data in IDL format:
                idl_mask = argv[2];

                IDL_ENSURE_SIMPLE(idl_mask);
                IDL_ENSURE_ARRAY(idl_mask);

                // Extract input in C format
                mask = (unsigned char *) idl_mask->value.arr->data;


                // Call Pore3D using the mask:
                p3dSijbersPostnovRingRemover3D_16(
                        in_rev16,
                        out_rev16,
                        (int) idl_in_rev->value.arr->dim[0],
                        (int) idl_in_rev->value.arr->dim[1],
                        (int) idl_in_rev->value.arr->dim[2],
                        centerX,
                        centerY,
                        winsize,
                        thresh,
                        iterations,
                        precision,
                        mask,
                        _p3d_idlPrintInfo
                        );
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
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } */
    else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D matrix.");
    }


    // Free resources:
    IDL_KW_FREE;

    // Return output in IDL Format
    return (idl_out_rev);
}

