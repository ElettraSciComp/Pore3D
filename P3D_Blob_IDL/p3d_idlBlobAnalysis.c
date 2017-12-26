// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dBlob.h"

// REV output struct in IDL format:
static IDL_MEMINT CC_tags_dims[IDL_MAX_ARRAY_DIM];

// Connected components statistics output struct in IDL format.
// NOTE: Order is crucial in reference to the C structure and 
// names should be in upper case format and alphabetical:
static IDL_STRUCT_TAG_DEF BlobStats_tags[] = {
	{ "ASPECT_RATIO", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "COUNT", 0, (void *) IDL_TYP_ULONG},
	{ "EQ_SPHERE", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "EXTENT", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "MAX_AXIS", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "MAX_SPHERE", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "MIN_AXIS", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "SPHERICITY", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ "VOLUME", CC_tags_dims, (void *) IDL_TYP_DOUBLE},
	{ 0}
};

IDL_VPTR p3d_idlBlobAnalysis(int argc, IDL_VPTR argv[], char* argk) {

	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
		IDL_VPTR idl_blob_im;
		IDL_VPTR idl_star_im;
		IDL_LONG conn;
		int cn_there;
		IDL_LONG max_rot;
		int mr_there;
		double resolution;
		int rs_there;
		IDL_LONG skip_borders;
	} KW_RESULT;

	// Alphabetical order is crucial:
	static IDL_KW_PAR kw_pars[] = {
		IDL_KW_FAST_SCAN,
		{ "AXIS_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_star_im)},
		{ "BLOB_IM", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(idl_blob_im)},
		{ "CONN", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(cn_there), (char*) IDL_KW_OFFSETOF(conn) },
		{ "MAX_ROT", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(mr_there), (char*) IDL_KW_OFFSETOF(max_rot) },
		{ "SKIPBORDERS", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(skip_borders)},        
		{ "VOXELSIZE", IDL_TYP_DOUBLE, 1, 0, (int*) IDL_KW_OFFSETOF(rs_there), (char*) IDL_KW_OFFSETOF(resolution)},
		{ NULL}
	};

	KW_RESULT kw;

	// IDL array input and IDL struct output:
	IDL_VPTR idl_in_rev, idl_blob_im, idl_star_im, idl_out_struct;
	unsigned char *in_rev;
	unsigned char *blob_im = NULL;
	unsigned char *star_im = NULL;
	//BlobStats *stats = (BlobStats*) malloc(sizeof(BlobStats));
	BlobStats stats;
	int keywords_ct = 0;

	double resolution = 1.0; // default value
	int max_rot = 1024; // default value

	IDL_MEMINT tmp_dims[IDL_MAX_ARRAY_DIM];
	IDL_MEMINT offset;
	double* d_tmp_ptr;
	//unsigned short* us_tmp_ptr;
	unsigned int* ui_tmp_ptr;
	//unsigned char* uc_tmp_ptr;
	char* s_data;

	unsigned int i; 
	int err_code, skip_borders;
	IDL_StructDefPtr s;
	int conn = CONN6;

	// Get input data in IDL format:
	idl_in_rev = argv[0];

	IDL_ENSURE_SIMPLE(idl_in_rev);
	IDL_ENSURE_ARRAY(idl_in_rev); // Get input data in IDL format:		


	// Process keywords:
	IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);

	skip_borders = (kw.skip_borders == 0) ? P3D_FALSE : P3D_TRUE;

	// Get the BLOB_IM pointer from keyword:
	if (kw.idl_blob_im) {
		// Simple way to delete current memory:
		IDL_StoreScalarZero(kw.idl_blob_im, IDL_TYP_LONG);

		// Create a new array:
		if (!(idl_in_rev->flags & IDL_V_TEMP)) {
			blob_im = (unsigned char*) IDL_MakeTempArray(
				IDL_TYP_BYTE,
				idl_in_rev->value.arr->n_dim,
				idl_in_rev->value.arr->dim,
				IDL_ARR_INI_ZERO,
				&idl_blob_im
				);
		}

		keywords_ct++;
	}

	// Get the STAR_IM pointer from keyword:
	if (kw.idl_star_im) {
		// Simple way to delete current memory:
		IDL_StoreScalarZero(kw.idl_star_im, IDL_TYP_LONG);

		// Create a new array:
		if (!(idl_in_rev->flags & IDL_V_TEMP)) {
			star_im = (unsigned char*) IDL_MakeTempArray(
				IDL_TYP_BYTE,
				idl_in_rev->value.arr->n_dim,
				idl_in_rev->value.arr->dim,
				IDL_ARR_INI_ZERO,
				&idl_star_im
				);
		}

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

	// Get the CONN input argument:
	if (kw.cn_there)
	{
		if ( idl_in_rev->value.arr->n_dim == 3 )
		{
			// Check values:
			if ( ( kw.conn != 6 ) && ( kw.conn != 18 ) && ( kw.conn != 26 ) )
				_p3d_idlPrintNamedError("CONN must be a value of the set {6,18,26}.");

			// Get values:
			if ( kw.conn == 6 )
				conn = CONN6;
			else if ( kw.conn == 18 )
				conn = CONN18;
			else if ( kw.conn == 26 )
				conn = CONN26;
		}
		// else: error on input arguments with further handling.

		keywords_ct++;
	}

	// Get the MAX_ROT input argument:
	if (kw.mr_there)
	{

		// Check values:
		if ( kw.max_rot < 0 )                
			_p3d_idlPrintNamedError("MAX_ROT must be greater than 0.");


		max_rot = kw.max_rot;		

		keywords_ct++;
	}



	if (idl_in_rev->value.arr->n_dim == 3) {
		// Extract input in C format
		if (idl_in_rev->type == IDL_TYP_BYTE)
			in_rev = (unsigned char *) idl_in_rev->value.arr->data;
		else
			_p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");

		// Call Pore3D:
		err_code = p3dBlobAnalysis(
			in_rev,
			&stats,
			blob_im,
			star_im,
			(int) idl_in_rev->value.arr->dim[0],
			(int) idl_in_rev->value.arr->dim[1],
			(int) idl_in_rev->value.arr->dim[2],
			resolution,
			conn,
			max_rot,
			skip_borders,
			_p3d_idlPrintInfo
			);

		if (err_code == P3D_ERROR) {
			// Print error:
			_p3d_idlPrintNamedError("Error on internal code execution.");			
		}

	} else {
		_p3d_idlPrintNamedError("Input argument IMAGE must be a 3D matrix.");		
	}


	// ---- Format struct for IDL output:

	// Now that we know dimension we can create output struct:
	CC_tags_dims[0] = 1;
	CC_tags_dims[1] = (IDL_MEMINT) (stats.blobCount);

	s = IDL_MakeStruct(NULL, BlobStats_tags);

	// Use temporary variable: doing so IDL will handle memory:
	tmp_dims[0] = 1;
	s_data = (char *) IDL_MakeTempStruct(s, 1, tmp_dims, &idl_out_struct, 0);



	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "ASPECT_RATIO", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.aspect_ratio[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "SPHERICITY", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.sphericity[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "COUNT", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	ui_tmp_ptr = (unsigned int *) (s_data + offset);
	// Store the double value into field: 
	*(ui_tmp_ptr) = stats.blobCount;

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "EQ_SPHERE", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.eq_sph[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "MAX_SPHERE", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.max_sph[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "MIN_AXIS", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.l_min[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "MAX_AXIS", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.l_max[i];

	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "EXTENT", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.extent[i];


	// Get the field of the structure:
	offset = IDL_StructTagInfoByName(s, "VOLUME", IDL_MSG_LONGJMP, NULL);
	// Get a pointer to that location:
	d_tmp_ptr = (double *) (s_data + offset);
	// Store the double values into array: 
	for (i = 0; i < stats.blobCount; i++)
		*(d_tmp_ptr++) = stats.volume[i];


	// Return keyword if required:
	if (kw.idl_blob_im) IDL_VarCopy(idl_blob_im, kw.idl_blob_im);
	if (kw.idl_star_im) IDL_VarCopy(idl_star_im, kw.idl_star_im);

	// Free C memory: data has been copied and will be handled by IDL:
	if (stats.volume != NULL) free(stats.volume);
	if (stats.aspect_ratio != NULL) free(stats.aspect_ratio);
	if (stats.eq_sph != NULL) free(stats.eq_sph);
	if (stats.max_sph != NULL) free(stats.max_sph);
	if (stats.l_min != NULL) free(stats.l_min);
	if (stats.l_max != NULL) free(stats.l_max);
	if (stats.extent != NULL) free(stats.extent);
	if (stats.sphericity != NULL) free(stats.sphericity);


	// Release keyword resources:
	IDL_KW_FREE;

	// Return output in IDL Format:
	return idl_out_struct;
}

