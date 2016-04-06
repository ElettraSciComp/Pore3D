// From C library:
#include <stdio.h>
#include <stdlib.h>

// From IDL export library:
#include "p3d_idlAdapter.h"


int IDL_Load(void)
{
	// These tables contain information on the functions and procedures
	// that make up the P3D_IDL.DLM. The information contained in these
	// tables must be identical to that contained in P3D_IDL.DLM. 
	// The syntax is
	//     FUNCTION RtnName [MinArgs] [MaxArgs] [Options...] 
	//     PROCEDURE RtnName [MinArgs] [MaxArgs] [Options...] 
	// where
	//     RtnName: The IDL user level name for the routine;
	//     MinArgs: The minimum number of arguments accepted by this 
	//              routine. If not supplied, 0 is assumed. 
	//     MaxArgs: The maximum number of arguments accepted by this 
	//              routine. If not supplied, 0 is assumed and the 
	//              special value IDL_MAXPARAMS can be used.
	//				Remarks: specifing correct values for MinArgs and
	//                       MaxArgs results in a compile-time error
	//                       checking. Any other error checking is 
	//                       performed at run-time.
	//     Options: Use it function/procedure accepts IDL keywords 
	//              (see documentation details).

	//
	// Functions:
	//
	static IDL_SYSFUN_DEF2 function_addr[] = {
		
		//
		// ---- Input / Output: ----
		//
                { (IDL_FUN_RET) p3d_idlBasicAnalysis, "P3DBASICANALYSIS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlAnisotropyAnalysis, "P3DANISOTROPYANALYSIS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlBlobLabeling, "P3DBLOBLABELING", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
		{ (IDL_FUN_RET) p3d_idlBlobAnalysis, "P3DBLOBANALYSIS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlChamferDT, "P3DCHAMFERDT",	1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlGetMaxVolumeBlob, "P3DGETMAXVOLUMEBLOB", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
		{ (IDL_FUN_RET) p3d_idlGetMinVolumeBlob, "P3DGETMINVOLUMEBLOB", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlMinVolumeFilter, "P3DMINVOLUMEFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlMorphometricAnalysis, "P3DMORPHOMETRICANALYSIS", 1, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlREVEstimation, "P3DREVESTIMATION", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlTextureAnalysis, "P3DTEXTUREANALYSIS", 1, 1, 0, 0 }
	};

	//
	// Procedures:
	//
	/*static IDL_SYSFUN_DEF2 procedure_addr[] = {  
		{ (IDL_SYSRTN_GENERIC) p3d_idlWriteRaw, "P3DWRITERAW", 2, 2, 0, 0} ,
	}; */

	
	// Register routines. The routines must be specified exactly the same as in DLM:	
	return IDL_SysRtnAdd(function_addr, TRUE, IDL_CARRAY_ELTS(function_addr));
		//&& IDL_SysRtnAdd(procedure_addr, FALSE, IDL_CARRAY_ELTS(procedure_addr));

}
