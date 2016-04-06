// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dBlob.h"

// REV output struct in IDL format:
static IDL_MEMINT REV_tags_dims[IDL_MAX_ARRAY_DIM];

// NOTE: Order is crucial in reference to the C structure and 
// names should be in upper case format and alphabetical:
static IDL_STRUCT_TAG_DEF REV_tags[] = {
    { "POROSITY", REV_tags_dims, (void *) IDL_TYP_DOUBLE},
    { "SIZES", REV_tags_dims, (void *) IDL_TYP_ULONG},
    { 0}
};

IDL_VPTR p3d_idlREVEstimation(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG cen_data[3];
        int cen_there;
        IDL_MEMINT cen_n;
        IDL_LONG stepsize;
        int ss_there;
    } KW_RESULT;

    static IDL_KW_ARR_DESC_R center = {(char*) IDL_KW_OFFSETOF(cen_data), 3, 3, (IDL_LONG*) IDL_KW_OFFSETOF(cen_n)};


    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "CENTER", IDL_TYP_LONG, 1, IDL_KW_ARRAY, (int*) IDL_KW_OFFSETOF(cen_there), (char*) (&center)},
        { "STEP", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(ss_there), (char*) IDL_KW_OFFSETOF(stepsize)},
        { NULL}
    };

    KW_RESULT kw;

    IDL_VPTR idl_in_rev, idl_out_struct;
    unsigned char *in_rev;

    double* rev_porosity;
    unsigned int* rev_sizes;

    unsigned int num_el;

    unsigned int stepsize = 10; // default
    unsigned int centerX = -1; // default
    unsigned int centerY = -1; // default
    unsigned int centerZ = -1; // default

    int keywords_ct = 0;
    int err_code;
    void* s;

    IDL_MEMINT tmp_dims[IDL_MAX_ARRAY_DIM];
    IDL_MEMINT offset;
    double* porosity_ptr;
    unsigned int* sizes_ptr;
    char* s_data;
    int i;

    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);

    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev);


    if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract input in C format
        if (idl_in_rev->type == IDL_TYP_BYTE)
            in_rev = (unsigned char *) idl_in_rev->value.arr->data;
        else
            _p3d_idlPrintNamedError("Input argument VOLUME must be of type BYTE.");
    } else {
        _p3d_idlPrintNamedError("Input argument VOLUME must be a 3D matrix.");
    }


    // Get CENTER keyword:
    if (kw.cen_there) {
        // Check values:
        if (kw.cen_n != 3)
            _p3d_idlPrintNamedError("Input argument CENTER must contain three [ X, Y, Z ] elements.");

        if ((kw.cen_data[0] < 0) || (kw.cen_data[0] > idl_in_rev->value.arr->dim[0]))
            _p3d_idlPrintNamedError("X value of input argument CENTER must be within IMAGE dimension.");

        if ((kw.cen_data[1] < 0) || (kw.cen_data[1] > idl_in_rev->value.arr->dim[1]))
            _p3d_idlPrintNamedError("Y value of input argument CENTER must be within IMAGE dimension.");

        if ((kw.cen_data[2] < 0) || (kw.cen_data[2] > idl_in_rev->value.arr->dim[2]))
            _p3d_idlPrintNamedError("Z value of input argument CENTER must be within IMAGE dimensions.");

        // Get values:
        centerX = (unsigned int) kw.cen_data[0];
        centerY = (unsigned int) kw.cen_data[1];
        centerZ = (unsigned int) kw.cen_data[2];
    }


    // Get the STEPSIZE input keyword:
    if (kw.ss_there) {
        // Check values:
        if (kw.stepsize < 0)
            _p3d_idlPrintNamedError("STEP must be greater than 0.");

        // Get value:
        stepsize = (unsigned int) kw.stepsize;

        keywords_ct++;
    }


    // Call Pore3D:
    err_code = p3dREVEstimation(
            in_rev,
            &rev_porosity,
            &rev_sizes,
            &num_el,
            (int) idl_in_rev->value.arr->dim[0],
            (int) idl_in_rev->value.arr->dim[1],
            (int) idl_in_rev->value.arr->dim[2],
            stepsize,
            centerX,
            centerY,
            centerZ,
            _p3d_idlPrintInfo
            );

    if (err_code == P3D_MEM_ERROR) {
        // Print error:
        _p3d_idlPrintNamedError("Error on internal code execution.");
    }




    // ---- Format struct for IDL output:

    // Now that we know dimension we can create output struct:
    REV_tags_dims[0] = 1;
    REV_tags_dims[1] = num_el;

    s = IDL_MakeStruct(NULL, REV_tags);

    // Use temporary variable: doing so IDL will handle memory:
    tmp_dims[0] = 1;

    s_data = (char *) IDL_MakeTempStruct(s, 1, tmp_dims, &idl_out_struct, 0);

    // Get the field of the structure:
    offset = IDL_StructTagInfoByName(s, "POROSITY", IDL_MSG_LONGJMP, NULL);
    // Get a pointer to that location:
    porosity_ptr = (double *) (s_data + offset);
    // Store the double values into array: 
    for (i = 0; i < num_el; i++)
        *(porosity_ptr++) = rev_porosity[i];

    // Get the field of the structure:
    offset = IDL_StructTagInfoByName(s, "SIZES", IDL_MSG_LONGJMP, NULL);
    // Get a pointer to that location:
    sizes_ptr = (unsigned int*) (s_data + offset);
    // Store the double values into array: 
    for (i = 0; i < num_el; i++)
        *(sizes_ptr++) = rev_sizes[i];

    // Free C memory: data has been copied and will be handled by IDL:
    free(rev_porosity);
    free(rev_sizes);
    
     // Release keyword resources:
    IDL_KW_FREE;

    // Return output in IDL Format:
    return idl_out_struct;
}

