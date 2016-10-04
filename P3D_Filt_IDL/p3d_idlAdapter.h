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

// From IDL export library:
#include <idl_export.h>

IDL_VPTR p3d_idlReadRaw8(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlReadRaw16(int, IDL_VPTR*, char* );
void p3d_idlWriteRaw(int, IDL_VPTR*, char* );


IDL_VPTR p3d_idlGaussianFilter(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlMeanFilter(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlMedianFilter(int, IDL_VPTR*, char* );

IDL_VPTR p3d_idlAnisotropicDiffusionFilter(int, IDL_VPTR*, char*);
IDL_VPTR p3d_idlBilateralFilter(int, IDL_VPTR*, char* );

IDL_VPTR p3d_idlBoinHaibelRingRemover(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlSijbersPostnovRingRemover(int, IDL_VPTR*, char* );

IDL_VPTR p3d_idlClearBorderFilter(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlCreateBinaryCircle(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlCreateBinaryCylinder(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlCreateBinarySphere(int, IDL_VPTR*, char* );
IDL_VPTR p3d_idlGetRegionByCoords(int, IDL_VPTR*, char* );

IDL_VPTR p3d_idlAutoThresholding(int, IDL_VPTR*, char* );

IDL_VPTR p3d_idlFrom16To8(int, IDL_VPTR*, char* );