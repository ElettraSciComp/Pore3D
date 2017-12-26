// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dBlob.h"

IDL_VPTR p3d_idlSquaredEuclideanDT(int argc, IDL_VPTR argv[], char* argk) {

    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev;
    unsigned int *out_rev;

    int err_code;


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev); // Get input data in IDL format:

    // Allocate memory for output:
    if (!(idl_in_rev->flags & IDL_V_TEMP))
        out_rev = (unsigned int *) IDL_MakeTempArray(
            IDL_TYP_ULONG,
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
        err_code = p3dSquaredEuclideanDT(
                in_rev,
                out_rev,
                (unsigned int) idl_in_rev->value.arr->dim[0],
                (unsigned int) idl_in_rev->value.arr->dim[1],
                (unsigned int) idl_in_rev->value.arr->dim[2],
                _p3d_idlPrintInfo
                );

        if (err_code == P3D_ERROR) {
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

    // Return output in IDL Format:
    return idl_out_rev;
}
