// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dFilt.h"

IDL_VPTR p3d_idlAutoThresholding(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG method;
        int mt_there;
        IDL_VPTR thresh;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "METHOD", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(mt_there), (char*) IDL_KW_OFFSETOF(method)},
        { "THRESH", IDL_TYP_LONG, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(thresh)},
        { NULL}
    };


    KW_RESULT kw;

    IDL_VPTR idl_out_rev, idl_in_rev;
    IDL_VPTR idl_thresh;
    unsigned char *in_rev8, *out_rev8;
    unsigned short *in_rev16;
    int keywords_ct = 0;
    int method = 1; // default = Otsu's
    unsigned char thresh8;
    unsigned short thresh16;

    int err_code;

    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    // Get the METHOD input argument:
    if (kw.mt_there) {
        
        if ((kw.method < 1) || (kw.method > 7))
            _p3d_idlPrintNamedError("METHOD must be an integer value within the range [1,7].");
          
        // Get values:
        method = kw.method; 
        
        keywords_ct++;
    }
    

    // Call Pore3D depending on input arguments:
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
            
            
            if ( method == 1 )
                err_code = p3dKittlerThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 2 )
                err_code = p3dOtsuThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 3 )
                err_code = p3dPunThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 4 )
                 err_code = p3dRidlerThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 5 )
                err_code = p3dKapurThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 6 )
                err_code = p3dJohannsenThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else if ( method == 7 )
                err_code = p3dHuangYagerThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
            else
                // Default case:
                err_code = p3dOtsuThresholding_8(in_rev8, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh8, _p3d_idlPrintInfo, NULL);
                        
            if (kw.thresh) IDL_StoreScalar(kw.thresh, IDL_TYP_BYTE, &thresh8);            

            // On exception print error:
            if (err_code == P3D_MEM_ERROR)
                _p3d_idlPrintNamedError("Error on code execution.");

        } else if (idl_in_rev->type == IDL_TYP_UINT) {
            in_rev16 = (unsigned short *) idl_in_rev->value.arr->data;

            // Allocate memory for output:
            if (!(idl_in_rev->flags & IDL_V_TEMP))
                out_rev8 = (unsigned short *) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,
                    idl_in_rev->value.arr->dim,
                    IDL_ARR_INI_NOP,
                    &idl_out_rev
                    );

            
            // Call Pore3D:
            switch (method) {
                case 1: err_code = p3dKittlerThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 2: err_code = p3dOtsuThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 3: err_code = p3dPunThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 4: err_code = p3dRidlerThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 5: err_code = p3dKapurThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 6: err_code = p3dJohannsenThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                case 7: err_code = p3dHuangYagerThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    break;
                /*default: 
                    _p3d_idlPrintInfo("\tUnsupported method specified. Default one will be used.");
                    err_code = p3dOtsuThresholding_16(in_rev16, out_rev8, (int) idl_in_rev->value.arr->dim[0], (int) idl_in_rev->value.arr->dim[1], (int) idl_in_rev->value.arr->dim[2], &thresh16, _p3d_idlPrintInfo, NULL);
                    //break;*/
            }

            if (kw.thresh) IDL_StoreScalar(kw.thresh, IDL_TYP_UINT, &thresh16);

            // On exception print error:
            if (err_code == P3D_MEM_ERROR)
                _p3d_idlPrintNamedError("Error on code execution.");

        } else {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE or UINT.");
        }
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
    }


    // Free resources:
    IDL_KW_FREE;

    // Return output in IDL Format
    return (idl_out_rev);
}
