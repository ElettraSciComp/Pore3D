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
