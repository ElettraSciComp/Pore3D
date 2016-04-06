// From C library:
#include <stdio.h>
#include <stdlib.h>

// Locals:
#include "_p3d_idlCommon.h"
#include "p3d_idlAdapter.h"

#include "p3dFilt.h"

IDL_VPTR p3d_idlGetRegionByCoords(int argc, IDL_VPTR argv[], char* argk)
{ 
    typedef struct {  
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_LONG conn;
        int cn_there; 
        IDL_LONG cen_data[3];  
        int cen_there;  
        IDL_MEMINT cen_n;  
    } KW_RESULT;  
 
    static IDL_KW_ARR_DESC_R center = { (char*) IDL_KW_OFFSETOF(cen_data), 2, 3, (IDL_LONG*) IDL_KW_OFFSETOF(cen_n) };
 
  
    // Alphabetical order is crucial:
    static IDL_KW_PAR kw_pars[] = {  
        IDL_KW_FAST_SCAN,  
        { "CONN", IDL_TYP_LONG, 1, 0, (int*) IDL_KW_OFFSETOF(cn_there), (char*) IDL_KW_OFFSETOF(conn) }, 
        { "COORDS", IDL_TYP_LONG, 1, IDL_KW_ARRAY, (int*) IDL_KW_OFFSETOF(cen_there), (char*) (&center) },
        { NULL }  
     };  
 
    KW_RESULT kw;
     
    IDL_VPTR idl_out_rev, idl_in_rev;
    unsigned char *in_rev8, *out_rev8;       
    int keywords_ct = 0;
 
    //int conn2D = CONN4;
    int conn3D = CONN6;
 
    unsigned int centerX, centerY, centerZ;
         
    int err_code;
 
     
    // Process keywords:
    IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw); 
 
     
    // Get input data in IDL format:
    idl_in_rev = argv[0]; 
 
    IDL_ENSURE_SIMPLE(idl_in_rev);  
    IDL_ENSURE_ARRAY(idl_in_rev);  
 
 
    // Get the CONN input argument:
    if (kw.cn_there)
    {
        /*if ( idl_in_rev->value.arr->n_dim == 2 )
        {
            // Check values:
            if ( ( kw.conn != 4 ) && ( kw.conn != 8) )
                _p3d_idlPrintNamedError("CONN must be a value of the set {4,8}.");
             
            // Get values:
            if ( kw.conn == 4 )
                conn2D = CONN4;
            else if ( kw.conn == 8 )
                conn2D = CONN8;
        }
        else **/if ( idl_in_rev->value.arr->n_dim == 3 )
        {
            // Check values:
            if ( ( kw.conn != 6 ) && ( kw.conn != 18 ) && ( kw.conn != 26 ) )
                _p3d_idlPrintNamedError("CONN must be a value of the set {6,18,26}.");
             
            // Get values:
            if ( kw.conn == 6 )
                conn3D = CONN6;
            else if ( kw.conn == 18 )
                conn3D = CONN18;
            else if ( kw.conn == 26 )
                conn3D = CONN26;
        }
        // else: error on input arguments with further handling.
 
        keywords_ct++;
    }
 
    // Get COORDS input argument:
    if (kw.cen_there)
    {
        // Check values:
        /*if (kw.cen_n == 2)
        {
            if ( ( kw.cen_data[0] < 0 ) || ( kw.cen_data[0] > idl_in_rev->value.arr->dim[0] ) )
                _p3d_idlPrintNamedError("X value of input argument COORDS must be within IMAGE dimensions.");
 
            if ( ( kw.cen_data[1] < 0 ) || ( kw.cen_data[1] > idl_in_rev->value.arr->dim[1] ) )
                _p3d_idlPrintNamedError("Y value of input argument COORDS must be within IMAGE dimensions.");
 
            // Get values:
            centerX = (unsigned int) kw.cen_data[0];
            centerY = (unsigned int) kw.cen_data[1];
        }
        else */if (kw.cen_n == 3)
        {
            if ( ( kw.cen_data[0] < 0 ) || ( kw.cen_data[0] > idl_in_rev->value.arr->dim[0] ) )
                _p3d_idlPrintNamedError("X value of input argument COORDS must be within IMAGE dimensions.");
 
            if ( ( kw.cen_data[1] < 0 ) || ( kw.cen_data[1] > idl_in_rev->value.arr->dim[1] ) )
                _p3d_idlPrintNamedError("Y value of input argument COORDS must be within IMAGE dimensions.");
 
            if ( ( kw.cen_data[2] < 0 ) || ( kw.cen_data[2] > idl_in_rev->value.arr->dim[2] ) )
                _p3d_idlPrintNamedError("Z value of input argument COORDS must be within IMAGE dimensions.");
 
            // Get values:
            centerX = (unsigned int) kw.cen_data[0];
            centerY = (unsigned int) kw.cen_data[1];
            centerZ = (unsigned int) kw.cen_data[2];
        }
        else
        {
            _p3d_idlPrintNamedError("Input argument COORDS must contain two [ X, Y ] elements.");
        }
 
        keywords_ct++;
    }
    else
    {
        /*if ( idl_in_rev->value.arr->n_dim == 2 )
        {
            // Set default values for centerX and centerY:
            centerX = (unsigned int) idl_in_rev->value.arr->dim[0] / 2;
            centerY = (unsigned int) idl_in_rev->value.arr->dim[1] / 2;
        }
        else*/ if ( idl_in_rev->value.arr->n_dim == 3 )
        {
            // Set default values for centerX, centerY and centerZ:
            centerX = (unsigned int) idl_in_rev->value.arr->dim[0] / 2;
            centerY = (unsigned int) idl_in_rev->value.arr->dim[1] / 2;
            centerZ = (unsigned int) idl_in_rev->value.arr->dim[2] / 2;
        }
        // else: error with further handling.
    }
 
 
    // Call Pore3D depending on input arguments:
    /*if ( idl_in_rev->value.arr->n_dim == 2 )
    {     
        // Extract first input (volume to filter) in C format:
        if (idl_in_rev->type == IDL_TYP_BYTE)             
        {
            in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;  
 
            // Allocate memory for output:
            if (!(idl_in_rev->flags & IDL_V_TEMP))  
                out_rev8 = (unsigned char *) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,  
                    idl_in_rev->value.arr->dim,  
                    IDL_ARR_INI_NOP, 
                    &idl_out_rev
                    );  
 
 
            // Call Pore3D:
            err_code = p3dGetRegionByCoords2D ( 
                in_rev8,
                out_rev8,
                (unsigned int) idl_in_rev->value.arr->dim[0],
                (unsigned int) idl_in_rev->value.arr->dim[1],
                centerX,
                centerY,
                conn2D,
                _p3d_idlPrintInfo
            );
 
            // On exception print error:
            if (err_code == P3D_MEM_ERROR)
                _p3d_idlPrintNamedError("Error on code execution.");  
        }         
        else
        {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");
        }         
    }
    else */if ( idl_in_rev->value.arr->n_dim == 3 )
    {
        // Extract first input (volume to filter) in C format:
        if (idl_in_rev->type == IDL_TYP_BYTE)             
        {
            in_rev8 = (unsigned char *) idl_in_rev->value.arr->data;  
 
            // Allocate memory for output:
            if (!(idl_in_rev->flags & IDL_V_TEMP))  
                out_rev8 = (unsigned char *) IDL_MakeTempArray(
                    IDL_TYP_BYTE,
                    idl_in_rev->value.arr->n_dim,  
                    idl_in_rev->value.arr->dim,  
                    IDL_ARR_INI_NOP, 
                    &idl_out_rev
                    ); 
 
             
            // Call Pore3D:
            err_code = p3dGetRegionByCoords3D ( 
                in_rev8,
                out_rev8,
                (unsigned int) idl_in_rev->value.arr->dim[0],
                (unsigned int) idl_in_rev->value.arr->dim[1],
                (unsigned int) idl_in_rev->value.arr->dim[2],
                centerX,
                centerY,
                centerZ,
                conn3D,
                _p3d_idlPrintInfo
            );     
 
            // On exception print error:
            if (err_code == P3D_MEM_ERROR)
                _p3d_idlPrintNamedError("Error on code execution."); 
        }         
        else
        {
            _p3d_idlPrintNamedError("Input argument IMAGE must be of type BYTE.");
        }
    }
    else
    {
        _p3d_idlPrintNamedError("Input argument IMAGE must be a 2D or 3D matrix.");
    }
 
 
    // Free resources:
    IDL_KW_FREE;
 
    // Return output in IDL Format
    return(idl_out_rev);  
}  