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
		{ (IDL_FUN_RET) p3d_idlReadRaw8, "P3DREADRAW8", 2, 2, 0, 0 },
		{ (IDL_FUN_RET) p3d_idlReadRaw16, "P3DREADRAW16", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlGaussianFilter, "P3DGAUSSIANFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlMeanFilter, "P3DMEANFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlMedianFilter, "P3DMEDIANFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlBilateralFilter, "P3DBILATERALFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlAnisotropicDiffusionFilter, "P3DANISOTROPICDIFFUSIONFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlBoinHaibelRingRemover, "P3DBOINHAIBELRINGREMOVER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },                
                { (IDL_FUN_RET) p3d_idlSijbersPostnovRingRemover, "P3DSIJBERSPOSTNOVRINGREMOVER", 1, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                
                { (IDL_FUN_RET) p3d_idlClearBorderFilter, "P3DCLEARBORDERFILTER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlGetRegionByCoords, "P3DGETREGIONBYCOORDS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                
      		{ (IDL_FUN_RET) p3d_idlCreateBinaryCircle, "P3DCREATEBINARYCIRCLE", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlCreateBinaryCylinder, "P3DCREATEBINARYCYLINDER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                { (IDL_FUN_RET) p3d_idlCreateBinarySphere, "P3DCREATEBINARYSPHERE", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
                
                { (IDL_FUN_RET) p3d_idlFrom16To8, "P3DFROM16TO8", 2, 2, 0, 0 },
                
                { (IDL_FUN_RET) p3d_idlAutoThresholding, "P3DAUTOTHRESHOLDING", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
	};

	//
	// Procedures:
	//
	static IDL_SYSFUN_DEF2 procedure_addr[] = {  		
               
                { (IDL_SYSRTN_GENERIC) p3d_idlWriteRaw, "P3DWRITERAW", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0} ,
	}; 

	
	// Register routines. The routines must be specified exactly the same as in DLM:	
	return IDL_SysRtnAdd(function_addr, TRUE, IDL_CARRAY_ELTS(function_addr))
		&& IDL_SysRtnAdd(procedure_addr, FALSE, IDL_CARRAY_ELTS(procedure_addr));

}
