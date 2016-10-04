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

#include "_p3d_idlCommon.h"

//
// Internals for error handling and verbose mode
//
int _p3d_idlPrintInfo ( const char* msg, ... )
{
	// Print a message to IDL console without interrupting code execution
	// (printing is inhibited if IDL environment variable !QUIET is set):
	va_list fmtargs;
	char buffer[1024];

	va_start ( fmtargs, msg );
	vsnprintf ( buffer, sizeof(buffer) - 1, msg, fmtargs );
	va_end ( fmtargs );

	IDL_Message (IDL_M_GENERIC, IDL_MSG_INFO, buffer ); 

	return 1;
}

int _p3d_idlPrintError ( const char* msg, ... )
{
	va_list fmtargs;
	char buffer[1024];

	va_start ( fmtargs, msg );
	vsnprintf ( buffer, sizeof(buffer) - 1, msg, fmtargs );
	va_end ( fmtargs );

	// Print a message to IDL console interrupting code execution:
	IDL_Message (IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer ); 

	return 1;
}

int _p3d_idlPrintNamedError ( const char* msg, ... )
{
	va_list fmtargs;
	char buffer[1024];

	va_start ( fmtargs, msg );
	vsnprintf ( buffer, sizeof(buffer) - 1, msg, fmtargs );
	va_end ( fmtargs );

	// Print a message to IDL console interrupting code execution:
	IDL_Message (IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, buffer ); 

	return 1;
}

double getIdlDoubleValue(IDL_VPTR id)
{

  double *idptr;
  IDL_MEMINT n;  
  IDL_VPTR idCvt;
  double val; 

  
  /* Extract data.  This way we don't demand scalar input, just simple 
     input.  We take the first input value */

  if (id->type != IDL_TYP_DOUBLE) {
    idCvt = IDL_BasicTypeConversion(1, &id, IDL_TYP_DOUBLE);
    IDL_VarGetData(idCvt, &n, (char **) &idptr, IDL_TRUE);
    val = idptr[0];
    IDL_Deltmp(idCvt);
  } else {
    IDL_VarGetData(id, &n, (char **) &idptr, IDL_TRUE);
    val = idptr[0];    
  }

  return val;  
}
