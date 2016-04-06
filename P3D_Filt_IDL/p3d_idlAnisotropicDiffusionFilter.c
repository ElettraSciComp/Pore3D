// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dFilt.h"

IDL_VPTR p3d_idlAnisotropicDiffusionFilter(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG m;
        int m_there;
        IDL_LONG iter;
        int it_there;
        double sigma;
        int s_there;
        double lambda;
        int l_there;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "ITERATIONS", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(it_there), (char*) IDL_KW_OFFSETOF(iter)},
        { "LAMBDA", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(l_there), (char*) IDL_KW_OFFSETOF(lambda)},
        { "MU", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(m_there), (char*) IDL_KW_OFFSETOF(m)},
        { "SIGMA", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(s_there), (char*) IDL_KW_OFFSETOF(sigma)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev8, *out_rev8;
    unsigned short *in_rev16, *out_rev16;
    int keywords_ct = 0;

    int iter      = 20;
    int m         = 1; // default
    double sigma  = 0.01; // default
    double lambda = 0.01; // default


    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    // Get the SIZE input argument:
    if (kw.m_there) {
        // Check values:
        if (kw.m <= 0)
            _p3d_idlPrintNamedError("MU must be an integer value greater than zero.");

        // Get values:
        m = (int) kw.m;

        keywords_ct++;
    }

    // Get the SIZE input argument:
    if (kw.it_there) {
        // Check values:
        if (kw.iter <= 0)
            _p3d_idlPrintNamedError("ITERATIONS must be greater than 0.");

        // Get values:
        iter = (int) kw.iter;

        keywords_ct++;
    }


    // Get the SIGMA_D input argument:
    if (kw.s_there) {
        // Check values:
        if (kw.sigma <= 0)
            _p3d_idlPrintNamedError("SIGMA must be greater than 0.");

        // Get values:
        sigma = (double) kw.sigma;

        keywords_ct++;
    }

    // Get the SIGMA_R input argument:
    if (kw.l_there) {
        // Check values:
        if (kw.lambda <= 0)
            _p3d_idlPrintNamedError("LAMBDA must be greater than 0.");

        // Get values:
        lambda = (double) kw.lambda;

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

            // Call Pore3D:
            p3dBilateralFilter2D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    size,
                    sigma_d,
                    sigma_r,
                    iter,
                    _p3d_idlPrintInfo,
                    NULL
                    );
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
            p3dBilateralFilter2D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    size,
                    sigma_d,
                    sigma_r,
                    iter,
                    _p3d_idlPrintInfo,
                    NULL
                    );
        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } else */
    if (idl_in_rev->value.arr->n_dim == 3) {
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
            p3dAnisotropicDiffusionFilter3D_8(
                    in_rev8,
                    out_rev8,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    m,
                    lambda,
                    sigma,
                    iter,
                    _p3d_idlPrintInfo,
                    NULL
                    );
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
            p3dAnisotropicDiffusionFilter3D_16(
                    in_rev16,
                    out_rev16,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    m,
                    lambda,
                    sigma,
                    iter,
                    _p3d_idlPrintInfo,
                    NULL
                    );
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
