// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dBlob.h"

// Morphometric statistics output struct in IDL format.
// NOTE: Order is crucial in reference to the C structure and 
// names should be in upper case format and alphabetical:
static IDL_STRUCT_TAG_DEF MorphStats_tags[] = {
    //{ "FRI", 0, (void *) IDL_TYP_DOUBLE },     
    { "BSBV", 0, (void *) IDL_TYP_DOUBLE },
    { "BVTV", 0, (void *) IDL_TYP_DOUBLE },
    { "TBN", 0, (void *) IDL_TYP_DOUBLE },
    { "TBSP", 0, (void *) IDL_TYP_DOUBLE },
    { "TBTH", 0, (void *) IDL_TYP_DOUBLE },
    { 0 }
};

IDL_VPTR p3d_idlMorphometricAnalysis(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        double resolution;
        int rs_there;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "VOXELSIZE", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(rs_there), (char*) IDL_KW_OFFSETOF(resolution)},
        { NULL}
    };

    KW_RESULT kw;

    // IDL array input and IDL struct output:
    IDL_VPTR idl_in_rev, idl_out_struct, idl_mask;
    unsigned char *in_rev, *mask;
    MorphometricStats stats;
    //MorphometricStats *stats = (MorphometricStats*) malloc(sizeof (MorphometricStats));
    int keywords_ct = 0;

    IDL_MEMINT tmp_dims[IDL_MAX_ARRAY_DIM];
    IDL_MEMINT offset;
    double* d_tmp_ptr;
    char* s_data;

    double resolution = 1.0; // default value

    int err_code;
    IDL_StructDefPtr s;

    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev); // Get input data in IDL format:		


    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);


    // Get the RESOLUTION input keyword:
    if (kw.rs_there) {
        // Check values:
        if (kw.resolution <= 0)
            _p3d_idlPrintNamedError("VOXELSIZE must be greater than 0.");

        // Get value:
        resolution = (double) kw.resolution;

        keywords_ct++;
    }


    // Create the structure definition  
    s = IDL_MakeStruct(NULL, MorphStats_tags);

    // Import the data area into an IDL structure. This code performs an 1-to-1 
    // mapping of IDL struct with the C struct:
    //idl_out_struct = IDL_ImportArray(1, &one, IDL_TYP_STRUCT, (UCHAR *) stats, NULL, s);
    //stats = (MorphometricStats *) idl_out_struct->value.s.arr->data;

    if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract input in C format
        if (idl_in_rev->type == IDL_TYP_BYTE)
            in_rev = (unsigned char *) idl_in_rev->value.arr->data;
        else
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");


        // Check if MASK is present:
        if ((argc - keywords_ct) == 1) {

            // No MASK:
            err_code = p3dMorphometricAnalysis(
                    in_rev,
                    NULL,
                    &stats,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    resolution,
                    _p3d_idlPrintInfo
                    );

        } else if ((argc - keywords_ct) == 2) {

            // Get input data in IDL format:
            idl_mask = argv[1];

            IDL_ENSURE_SIMPLE(idl_mask);
            IDL_ENSURE_ARRAY(idl_mask);

            // Extract input in C format
            mask = (unsigned char *) idl_mask->value.arr->data;


            // Call Pore3D:
            err_code = p3dMorphometricAnalysis(
                    in_rev,
                    mask,
                    &stats,
                    (int) idl_in_rev->value.arr->dim[0],
                    (int) idl_in_rev->value.arr->dim[1],
                    (int) idl_in_rev->value.arr->dim[2],
                    resolution,
                    _p3d_idlPrintInfo
                    );

        }


        if (err_code == P3D_ERROR) {
            // Print error:
            _p3d_idlPrintNamedError("Error on internal code execution.");
        }

    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
    }

    // -------------------------------------------------------------------------
    s = IDL_MakeStruct(NULL, MorphStats_tags);

    // Use temporary variable: doing so IDL will handle memory:
    tmp_dims[0] = 1;
    s_data = (char *) IDL_MakeTempStruct(s, 1, tmp_dims, &idl_out_struct, 0);


    // Get the field of the structure:
    offset = IDL_StructTagInfoByName(s, "BVTV", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr) = stats.BvTv;

    offset = IDL_StructTagInfoByName(s, "BSBV", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr) = stats.BsBv;

    offset = IDL_StructTagInfoByName(s, "TBN", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr) = stats.TbN;

    offset = IDL_StructTagInfoByName(s, "TBTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr) = stats.TbTh;

    offset = IDL_StructTagInfoByName(s, "TBSP", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr) = stats.TbSp;


    // Release keyword resources:
    IDL_KW_FREE;

    // Return output in IDL Format:
    return idl_out_struct;
}


