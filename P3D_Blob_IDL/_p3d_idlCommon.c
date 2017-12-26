// From C library:
#include <stdio.h>
#include <stdlib.h>

// From IDL export library:
#include <idl_export.h>

#include "_p3d_idlCommon.h"

//
// Internals for error handling and verbose mode
//
int _p3d_idlPrintInfo ( const char* msg, ... )
{
	// Print a message to IDL console without interrupting code execution
	// (printing is inhibited if IDL environment variable !QUIET is set):
	va_list fmtargs;
	char buffer[4096];

	va_start ( fmtargs, msg );
	vsnprintf ( buffer, sizeof(buffer) - 1, msg, fmtargs );
	va_end ( fmtargs );

	IDL_Message (IDL_M_GENERIC, IDL_MSG_INFO, buffer ); 

	return 1;
}

int _p3d_idlPrintError ( const char* msg, ... )
{
	va_list fmtargs;
	char buffer[4096];

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
	char buffer[4096];

	va_start ( fmtargs, msg );
	vsnprintf ( buffer, sizeof(buffer) - 1, msg, fmtargs );
	va_end ( fmtargs );

	// Print a message to IDL console interrupting code execution:
	IDL_Message (IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, buffer ); 

	return 1;
}