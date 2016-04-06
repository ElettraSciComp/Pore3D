
// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dITKWrapped.h"

IDL_VPTR p3d_idlGradientAnisotropicDiffusionFilter(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG iterations;
        int it_there;
        double conductance;
        int cd_there;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "ITERATIONS", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(it_there), (char*) IDL_KW_OFFSETOF(iterations)},
        { "K", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(cd_there), (char*) IDL_KW_OFFSETOF(conductance)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev8, *out_rev8;
    unsigned short *in_rev16, *out_rev16;
    int keywords_ct = 0;

    int iterations = 5; // default
    double conductance = 3.0; // default

    int err_code;

    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    // Get the ITERATIONS input argument:
    if (kw.it_there) {
        // Check values:
        if ((kw.iterations < 1) || (kw.iterations > 120))
            _p3d_idlPrintNamedError("ITERATIONS must be an integer value within the range [1,120].");

        // Get values:
        iterations = (int) kw.iterations;

        keywords_ct++;
    }


    // Get the SIGMA_D input argument:
    if (kw.cd_there) {
        // Check values:
        if (kw.conductance < 0)
            _p3d_idlPrintNamedError("K must be greater than 0.");

        // Get values:
        conductance = (double) kw.conductance;

        keywords_ct++;
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


            // Call Pore3D:
            err_code = p3dGradientAnisotropicDiffusionFilter2D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    iterations,
                    conductance,
                    _p3d_idlPrintInfo
                    );

            // On exception print error:
            if (err_code == P3D_ERROR)
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
            err_code = p3dGradientAnisotropicDiffusionFilter2D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    iterations,
                    conductance,
                    _p3d_idlPrintInfo
                    );

            // On exception print error:
            if (err_code == P3D_ERROR)
                _p3d_idlPrintNamedError("Error on code execution.");

        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } else if (idl_in_rev->value.arr->n_dim == 3) {
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
            err_code = p3dGradientAnisotropicDiffusionFilter3D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    iterations,
                    conductance,
                    _p3d_idlPrintInfo
                    );

            // On exception print error:
            if (err_code == P3D_ERROR)
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
            err_code = p3dGradientAnisotropicDiffusionFilter3D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    iterations,
                    conductance,
                    _p3d_idlPrintInfo
                    );

            // On exception print error:
            if (err_code == P3D_ERROR)
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
