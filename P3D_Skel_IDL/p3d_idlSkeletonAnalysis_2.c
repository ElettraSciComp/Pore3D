// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dSkel.h"

// REV output struct in IDL format:
static IDL_MEMINT Skel_tags_node_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_end_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_node2node_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_node2end_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_end2end_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_coordinationnumber_dims[IDL_MAX_ARRAY_DIM];
static IDL_MEMINT Skel_tags_tortuosity_dims[IDL_MAX_ARRAY_DIM];

// NOTE: Order is crucial in reference to the C structure and 
// names should be in upper case format and alphabetical:
static IDL_STRUCT_TAG_DEF SkelStats_tags[] = {
    { "CONNECTIVITY_DENSITY", 0, (void *) IDL_TYP_DOUBLE},
    { "COORDINATION_NUMBER", Skel_tags_coordinationnumber_dims, (void *) IDL_TYP_LONG},
    { "ENDPOINTS_COUNT", 0, (void *) IDL_TYP_UINT},
    { "ENDPOINTS_WIDTH", Skel_tags_end_dims, (void *) IDL_TYP_DOUBLE},    
    { "ENDTOEND_COUNT", 0, (void *) IDL_TYP_UINT},
    { "ENDTOEND_LENGTH", Skel_tags_end2end_dims, (void *) IDL_TYP_DOUBLE},
    { "ENDTOEND_MAXWIDTH", Skel_tags_end2end_dims, (void *) IDL_TYP_DOUBLE},
    { "ENDTOEND_MEANWIDTH", Skel_tags_end2end_dims, (void *) IDL_TYP_DOUBLE},
    { "ENDTOEND_MINWIDTH", Skel_tags_end2end_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETOEND_COUNT", 0, (void *) IDL_TYP_UINT},
    { "NODETOEND_LENGTH", Skel_tags_node2end_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETOEND_MAXWIDTH", Skel_tags_node2end_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETOEND_MEANWIDTH", Skel_tags_node2end_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETOEND_MINWIDTH", Skel_tags_node2end_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETONODE_COUNT", 0, (void *) IDL_TYP_UINT},
    { "NODETONODE_LENGTH", Skel_tags_node2node_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETONODE_MAXWIDTH", Skel_tags_node2node_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETONODE_MEANWIDTH", Skel_tags_node2node_dims, (void *) IDL_TYP_DOUBLE},
    { "NODETONODE_MINWIDTH", Skel_tags_node2node_dims, (void *) IDL_TYP_DOUBLE},
    { "PORES_COUNT", 0, (void *) IDL_TYP_UINT},
    { "PORES_WIDTH", Skel_tags_node_dims, (void *) IDL_TYP_DOUBLE},
    { "TORTUOSITY_COUNT", 0, (void *) IDL_TYP_UINT},
    { "TORTUOSITY_X", Skel_tags_tortuosity_dims, (void *) IDL_TYP_DOUBLE},
    { "TORTUOSITY_Y", Skel_tags_tortuosity_dims, (void *) IDL_TYP_DOUBLE},     
    { "TORTUOSITY_Z", Skel_tags_tortuosity_dims, (void *) IDL_TYP_DOUBLE},
    { 0 }
};

// Function:		P3DSKELETONSTATISTICS
// Description:	    Computes basic statistics of the input binary volume
// IDL Syntax:		Result = P3DSKELETONSTATISTICS ( VOLUME [, STEPSIZE ] [, CENTER ] )
// Return value:	A struct containing the following values:
// Arguments:		VOLUME:				the 3D matrix of the binary volume
//					RESOLUTION:			the strength of the potential field (default: 6)
// Example:			SKEL = P3DSKELETONSTATISTICS ( VOL, 5, [ 100, 100, 100 ] )

IDL_VPTR p3d_idlSkeletonAnalysis_2(int argc, IDL_VPTR argv[], char* argk) {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_VPTR idl_ends_im;
        double merging_factor;
        int mf_there;
        IDL_VPTR idl_nodes_im;
        IDL_VPTR idl_pores_im;
        IDL_VPTR idl_throats_im;
        int tortuosity_depth;
        int tort_there;
        double resolution;
        int rs_there;
    } KW_RESULT;

    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { "ENDS_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_ends_im)},
        { "MERGING_FACTOR", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(mf_there), (char*) IDL_KW_OFFSETOF(merging_factor)},
        { "NODES_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_nodes_im)},
        { "PORES_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_pores_im)},        
        { "THROATS_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_throats_im)},
        { "TORTUOSITY_DEPTH", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(tort_there), (char*) IDL_KW_OFFSETOF(tortuosity_depth)},
        { "VOXELSIZE", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(rs_there), (char*) IDL_KW_OFFSETOF(resolution)},
        { NULL}
    };

    KW_RESULT kw;


    IDL_VPTR idl_in_rev, idl_skl_rev, idl_out_struct;
    unsigned char *in_rev, *skl_rev;

    unsigned char *nodes_im = NULL;
    unsigned char *pores_im = NULL;
    unsigned char *ends_im = NULL;
    unsigned char *throats_im = NULL;
    IDL_VPTR idl_nodes_im, idl_pores_im, idl_ends_im, idl_throats_im;

    //SkeletonStats *stats = (SkeletonStats*) malloc(sizeof(SkeletonStats));
    struct SkeletonStats stats;
    int keywords_ct = 0;

    double merging_factor = 0.85; // default value
    double resolution = 1.0; // default value
    int tortuosity_depth = 3; //defaul value
    
    int err_code;
    void* s;


    IDL_MEMINT tmp_dims[IDL_MAX_ARRAY_DIM];
    IDL_MEMINT offset;
    double* d_tmp_ptr;
    int* i_tmp_ptr;
    unsigned short* us_tmp_ptr;
    char* s_data;

    unsigned short i;

    // Init skeleton stats:
    stats.Node_Width = NULL;
    stats.End_Width = NULL;

    stats.EndToEnd_Length = NULL;
    stats.EndToEnd_MinWidth = NULL;
    stats.EndToEnd_MeanWidth = NULL;
    stats.EndToEnd_MaxWidth = NULL;

    stats.NodeToEnd_Length = NULL;
    stats.NodeToEnd_MinWidth = NULL;
    stats.NodeToEnd_MeanWidth = NULL;
    stats.NodeToEnd_MaxWidth = NULL;

    stats.NodeToNode_Length = NULL;
    stats.NodeToNode_MinWidth = NULL;
    stats.NodeToNode_MeanWidth = NULL;
    stats.NodeToNode_MaxWidth = NULL;

    stats.CoordinationNumber = NULL;
    stats.Tort_X = NULL;
    stats.Tort_Y = NULL;
    stats.Tort_Z = NULL;

    // Get input data in IDL format:
    idl_in_rev = argv[0];

    IDL_ENSURE_SIMPLE(idl_in_rev);
    IDL_ENSURE_ARRAY(idl_in_rev); // Get input data in IDL format:

    // Get input data in IDL format:
    idl_skl_rev = argv[1];

    IDL_ENSURE_SIMPLE(idl_skl_rev);
    IDL_ENSURE_ARRAY(idl_skl_rev); // Get input data in IDL format:

    // Check input format:
    if (idl_in_rev->value.arr->n_dim != idl_skl_rev->value.arr->n_dim)
        _p3d_idlPrintNamedError("IMAGE argument and SKELETON argument must be 3D matrices with identical dimensions.");

    if (idl_in_rev->value.arr->dim[0] != idl_skl_rev->value.arr->dim[0])
        _p3d_idlPrintNamedError("IMAGE argument and SKELETON argument must be 3D matrices with identical dimensions.");

    if (idl_in_rev->value.arr->dim[1] != idl_skl_rev->value.arr->dim[1])
        _p3d_idlPrintNamedError("IMAGE argument and SKELETON argument must be 3D matrices with identical dimensions.");

    if (idl_in_rev->value.arr->dim[2] != idl_skl_rev->value.arr->dim[2])
        _p3d_idlPrintNamedError("IMAGE argument and SKELETON argument must be 3D matrices with identical dimensions.");


    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);

    // Get the skip_borders keyword:
    //skip_borders = (kw.skip_borders == 0) ? P3D_FALSE : P3D_TRUE;

    // Get the MERGING_FACTOR input keyword:
    if (kw.mf_there) {
        // Check values:
        if ( (kw.merging_factor < 0.0) || (kw.merging_factor > 1.0) )
            _p3d_idlPrintNamedError("MERGING_FACTOR must be in the range [0.0,1.0]");

        // Get value:
        merging_factor = (double) kw.merging_factor;

        keywords_ct++;
    }
    
    // Get the TORTUOSITY DEPTH input keyword:
    if (kw.tort_there) {
        // Check values:
        if ( (kw.tortuosity_depth < 1) || (kw.tortuosity_depth > 10) )
            _p3d_idlPrintNamedError("TORTUOSITY_DEPTH must be an integer value in the range [1,10]");

        // Get value:
        tortuosity_depth = (int) kw.tortuosity_depth;

        keywords_ct++;
    }
    
    
    // Get the RESOLUTION input keyword:
    if (kw.rs_there) {
        // Check values:
        if (kw.resolution <= 0)
            _p3d_idlPrintNamedError("VOXELSIZE must be greater than 0.");

        // Get value:
        resolution = (double) kw.resolution;

        keywords_ct++;
    }

    // Get the NODES_IM pointer from keyword:
    if (kw.idl_nodes_im) {
        // Simple way to delete current memory:
        IDL_StoreScalarZero(kw.idl_nodes_im, IDL_TYP_LONG);

        // Create a new array:
        if (!(idl_in_rev->flags & IDL_V_TEMP)) {
            nodes_im = (unsigned char*) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,
                    idl_in_rev->value.arr->dim,
                    IDL_ARR_INI_ZERO,
                    &idl_nodes_im
                    );
        }

        keywords_ct++;
    }

    // Get the ENDS_IM pointer from keyword:
    if (kw.idl_ends_im) {
        // Simple way to delete current memory:
        IDL_StoreScalarZero(kw.idl_ends_im, IDL_TYP_LONG);

        // Create a new array:
        if (!(idl_in_rev->flags & IDL_V_TEMP)) {
            ends_im = (unsigned char*) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,
                    idl_in_rev->value.arr->dim,
                    IDL_ARR_INI_ZERO,
                    &idl_ends_im
                    );
        }

        keywords_ct++;
    }
    
     // Get the PORES_IM pointer from keyword:
    if (kw.idl_pores_im) {
        // Simple way to delete current memory:
        IDL_StoreScalarZero(kw.idl_pores_im, IDL_TYP_LONG);

        // Create a new array:
        if (!(idl_in_rev->flags & IDL_V_TEMP)) {
            pores_im = (unsigned char*) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,
                    idl_in_rev->value.arr->dim,
                    IDL_ARR_INI_ZERO,
                    &idl_pores_im
                    );
        }

        keywords_ct++;
    }

    // Get the THROATS_IM pointer from keyword:
    if (kw.idl_throats_im) {
        // Simple way to delete current memory:
        IDL_StoreScalarZero(kw.idl_throats_im, IDL_TYP_LONG);

        // Create a new array:
        if (!(idl_in_rev->flags & IDL_V_TEMP)) {
            throats_im = (unsigned char*) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,
                    idl_in_rev->value.arr->dim,
                    IDL_ARR_INI_ZERO,
                    &idl_throats_im
                    );
        }

        keywords_ct++;
    }



    if (idl_in_rev->value.arr->n_dim == 3) {
        // Extract input in C format
        if (idl_in_rev->type == IDL_TYP_BYTE)
            in_rev = (unsigned char *) idl_in_rev->value.arr->data;
        else
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");

        // Extract input in C format
        if (idl_skl_rev->type == IDL_TYP_BYTE)
            skl_rev = (unsigned char *) idl_skl_rev->value.arr->data;
        else
            _p3d_idlPrintNamedError("Input argument SKELETON must be of type BYTE.");


        // Call Pore3D:
        err_code = p3dSkeletonAnalysis_2(
                in_rev,
                skl_rev,
                &stats,
                nodes_im,
                pores_im,
                ends_im,
                throats_im,
                idl_in_rev->value.arr->dim[0],
                idl_in_rev->value.arr->dim[1],
                idl_in_rev->value.arr->dim[2],
                merging_factor, 
                tortuosity_depth,
                resolution,
                _p3d_idlPrintInfo
                );

        if (err_code == P3D_MEM_ERROR) {
            // Print error:
            _p3d_idlPrintNamedError("Error on internal code execution.");
        }
    } else {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");
    }

    // ---- Format struct for IDL output:

    Skel_tags_end_dims[0] = 1;
    Skel_tags_end_dims[1] = (IDL_MEMINT) stats.End_Counter;
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_end_dims[i] = 0;

    Skel_tags_end2end_dims[0] = 1;
    Skel_tags_end2end_dims[1] = (IDL_MEMINT) stats.EndToEnd_Counter;
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_end2end_dims[i] = 0;

    Skel_tags_node_dims[0] = 1;
    Skel_tags_node_dims[1] = (IDL_MEMINT) stats.Node_Counter;
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_node_dims[i] = 0;

    Skel_tags_coordinationnumber_dims[0] = 1;
    Skel_tags_coordinationnumber_dims[1] = (IDL_MEMINT) stats.Node_Counter; 
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_coordinationnumber_dims[i] = 0;    
   
    Skel_tags_tortuosity_dims[0] = 1;
    Skel_tags_tortuosity_dims[1] = (IDL_MEMINT) stats.Tort_Counter;     
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_tortuosity_dims[i] = 0;

    Skel_tags_node2end_dims[0] = 1;
    Skel_tags_node2end_dims[1] = (IDL_MEMINT) stats.NodeToEnd_Counter;
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_node2end_dims[i] = 0;

    Skel_tags_node2node_dims[0] = 1;
    Skel_tags_node2node_dims[1] = (IDL_MEMINT) stats.NodeToNode_Counter;
    for (i = 2; i < IDL_MAX_ARRAY_DIM; i++)
        Skel_tags_node2node_dims[i] = 0;


    // Create an anonymous structure:
    s = IDL_MakeStruct(NULL, SkelStats_tags);
    tmp_dims[0] = 1;
    for (i = 1; i < IDL_MAX_ARRAY_DIM; i++)
		tmp_dims[i] = 0;
    s_data = (char *) IDL_MakeTempStruct(s, 1, tmp_dims, &idl_out_struct, 0);
    
    
    offset = IDL_StructTagInfoByName(s, "CONNECTIVITY_DENSITY", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    *(d_tmp_ptr++) = stats.ConnectivityDensity;

    // Copy values:
    offset = IDL_StructTagInfoByName(s, "COORDINATION_NUMBER", IDL_MSG_LONGJMP, NULL);
    i_tmp_ptr = (int *) (s_data + offset);
    for (i = 0; i < stats.Node_Counter; i++)
        *(i_tmp_ptr++) = stats.CoordinationNumber[i];
    
    offset = IDL_StructTagInfoByName(s, "TORTUOSITY_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.Tort_Counter; 
    
  
    offset = IDL_StructTagInfoByName(s, "TORTUOSITY_X", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.Tort_Counter; i++)
        *(d_tmp_ptr++) = stats.Tort_X[i];
    
    offset = IDL_StructTagInfoByName(s, "TORTUOSITY_Y", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.Tort_Counter; i++)
        *(d_tmp_ptr++) = stats.Tort_Y[i];
    
    offset = IDL_StructTagInfoByName(s, "TORTUOSITY_Z", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.Tort_Counter; i++)
        *(d_tmp_ptr++) = stats.Tort_Z[i];
    
    offset = IDL_StructTagInfoByName(s, "ENDPOINTS_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.End_Counter;

    offset = IDL_StructTagInfoByName(s, "ENDPOINTS_WIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.End_Counter; i++)
        *(d_tmp_ptr++) = stats.End_Width[i];


    offset = IDL_StructTagInfoByName(s, "PORES_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.Node_Counter;

    offset = IDL_StructTagInfoByName(s, "PORES_WIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.Node_Counter; i++)
        *(d_tmp_ptr++) = stats.Node_Width[i];

    offset = IDL_StructTagInfoByName(s, "ENDTOEND_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.EndToEnd_Counter;

    offset = IDL_StructTagInfoByName(s, "ENDTOEND_LENGTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.EndToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.EndToEnd_Length[i];

    offset = IDL_StructTagInfoByName(s, "ENDTOEND_MINWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.EndToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.EndToEnd_MinWidth[i];

    offset = IDL_StructTagInfoByName(s, "ENDTOEND_MEANWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.EndToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.EndToEnd_MeanWidth[i];

    offset = IDL_StructTagInfoByName(s, "ENDTOEND_MAXWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.EndToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.EndToEnd_MaxWidth[i];




    offset = IDL_StructTagInfoByName(s, "NODETOEND_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.NodeToEnd_Counter;

    offset = IDL_StructTagInfoByName(s, "NODETOEND_LENGTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToEnd_Length[i];

    offset = IDL_StructTagInfoByName(s, "NODETOEND_MINWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToEnd_MinWidth[i];

    offset = IDL_StructTagInfoByName(s, "NODETOEND_MEANWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToEnd_MeanWidth[i];

    offset = IDL_StructTagInfoByName(s, "NODETOEND_MAXWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToEnd_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToEnd_MaxWidth[i];



    offset = IDL_StructTagInfoByName(s, "NODETONODE_COUNT", IDL_MSG_LONGJMP, NULL);
    us_tmp_ptr = (unsigned short *) (s_data + offset);
    *(us_tmp_ptr++) = stats.NodeToNode_Counter;


    offset = IDL_StructTagInfoByName(s, "NODETONODE_LENGTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToNode_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToNode_Length[i];

    offset = IDL_StructTagInfoByName(s, "NODETONODE_MINWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToNode_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToNode_MinWidth[i];

    offset = IDL_StructTagInfoByName(s, "NODETONODE_MEANWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToNode_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToNode_MeanWidth[i];

    offset = IDL_StructTagInfoByName(s, "NODETONODE_MAXWIDTH", IDL_MSG_LONGJMP, NULL);
    d_tmp_ptr = (double *) (s_data + offset);
    for (i = 0; i < stats.NodeToNode_Counter; i++)
        *(d_tmp_ptr++) = stats.NodeToNode_MaxWidth[i];

    // Return keyword if required:
    if (kw.idl_nodes_im) IDL_VarCopy(idl_nodes_im, kw.idl_nodes_im);
    if (kw.idl_pores_im) IDL_VarCopy(idl_pores_im, kw.idl_pores_im);
    if (kw.idl_ends_im) IDL_VarCopy(idl_ends_im, kw.idl_ends_im);
    if (kw.idl_throats_im) IDL_VarCopy(idl_throats_im, kw.idl_throats_im);


    // Free C memory: data has been copied and will be handled by IDL:
    if (stats.Node_Width != NULL) free(stats.Node_Width);
    if (stats.End_Width != NULL) free(stats.End_Width);

    if (stats.EndToEnd_Length != NULL) free(stats.EndToEnd_Length);
    if (stats.EndToEnd_MinWidth != NULL) free(stats.EndToEnd_MinWidth);
    if (stats.EndToEnd_MeanWidth != NULL) free(stats.EndToEnd_MeanWidth);
    if (stats.EndToEnd_MaxWidth != NULL) free(stats.EndToEnd_MaxWidth);

    if (stats.NodeToEnd_Length != NULL) free(stats.NodeToEnd_Length);
    if (stats.NodeToEnd_MinWidth != NULL) free(stats.NodeToEnd_MinWidth);
    if (stats.NodeToEnd_MeanWidth != NULL) free(stats.NodeToEnd_MeanWidth);
    if (stats.NodeToEnd_MaxWidth != NULL) free(stats.NodeToEnd_MaxWidth);

    if (stats.NodeToNode_Length != NULL) free(stats.NodeToNode_Length);
    if (stats.NodeToNode_MinWidth != NULL) free(stats.NodeToNode_MinWidth);
    if (stats.NodeToNode_MeanWidth != NULL) free(stats.NodeToNode_MeanWidth);
    if (stats.NodeToNode_MaxWidth != NULL) free(stats.NodeToNode_MaxWidth);

    if (stats.CoordinationNumber != NULL) free(stats.CoordinationNumber);
    
    if (stats.Tort_X != NULL) free(stats.Tort_X);
    if (stats.Tort_Y != NULL) free(stats.Tort_Y);
    if (stats.Tort_Z != NULL) free(stats.Tort_Z);

    // Release keyword resources:
    IDL_KW_FREE;

    // Return output in IDL Format:
    return idl_out_struct;
} 

