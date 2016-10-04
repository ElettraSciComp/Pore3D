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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits.h>

#include "p3dFilt.h"
#include "p3dTime.h"

#include "Common/p3dCoordsQueue.h"

// Temporary label:
#define _GETREGIONBYCOORDS_TEMP_LABEL 2
#define _GETREGIONBYCOORDS_DEF_LABEL  1
 
void _p3dGetRegionByCoords3D_putCoordsInQueue(
        unsigned char* out_rev,
        int dimx,
        int dimy,
        int dimz,
        int i,
        int j,
        int k,
        int ct,
        int win_size,
        unsigned short m,
        coords_queue_t* queue // FIFO data structure for coords storage
        ) {
    coords_t tmp_coords; // Temporary coords


    // If we're on an object voxel:
    if (out_rev[ I(i, j, k, dimx, dimy) ] == OBJECT) {
        // Check if this is the last step:
        if (ct == (win_size - 2)) {
            // Border voxels are set to m:
            out_rev[ I(i, j, k, dimx, dimy) ] = m;

            // A border voxel is pushed into queue:
            tmp_coords.x = i;
            tmp_coords.y = j;
            tmp_coords.z = k;
            coords_queue_push(queue, tmp_coords);
        } else {
            // Mark voxel for further scan:
            out_rev[ I(i, j, k, dimx, dimy) ] = _GETREGIONBYCOORDS_TEMP_LABEL;
        }
    }
}

void _p3dGetRegionByCoords3D_conn_fun6(
        unsigned char* out_rev,
        int dimx,
        int dimy,
        int dimz,
        int i,
        int j,
        int k,
        int ct,
        int win_size,
        unsigned short m,
        coords_queue_t* queue // FIFO data structure for coords storage
        ) {
    int a, b, c;

    c = k - 1;
    b = j;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j - 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j + 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

}

void _p3dGetRegionByCoords3D_conn_fun18(
        unsigned char* out_rev,
        int dimx,
        int dimy,
        int dimz,
        int i,
        int j,
        int k,
        int ct,
        int win_size,
        unsigned short m,
        coords_queue_t* queue // FIFO data structure for coords storage
        ) {
    int a, b, c;

    // Perform 6-connectivity:
    _p3dGetRegionByCoords3D_conn_fun6(out_rev, dimx, dimy, dimz, i, j, k, ct, win_size, m, queue);

    // Do other 12 tests:

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j + 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j + 1;
    a = i;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);


    // On k plane:
    c = k;
    b = j - 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j - 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j + 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k;
    b = j + 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);


}

void _p3dGetRegionByCoords3D_conn_fun26(
        unsigned char* out_rev,
        int dimx,
        int dimy,
        int dimz,
        int i,
        int j,
        int k,
        int ct,
        int win_size,
        unsigned short m,
        coords_queue_t* queue // FIFO data structure for coords storage
        ) {
    int a, b, c;

    // Perform 18-connectivity:
    _p3dGetRegionByCoords3D_conn_fun18(out_rev, dimx, dimy, dimz, i, j, k, ct, win_size, m, queue);

    // Do other 8 tests:

    // On k-1 plane:
    c = k - 1;
    b = j - 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j - 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j + 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k - 1;
    b = j + 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);


    // On k+1 plane:
    c = k + 1;
    b = j - 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j - 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j + 1;
    a = i - 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

    c = k + 1;
    b = j + 1;
    a = i + 1;
    _p3dGetRegionByCoords3D_putCoordsInQueue(out_rev, dimx, dimy, dimz, a, b, c, ct, win_size, m, queue);

}
 

int p3dGetRegionByCoords3D (   
    unsigned char* in_rev, 
    unsigned char* out_rev, 
    const int dimx,
    const int dimy, 
    const int dimz, 
    const int x,
    const int y,
    const int z,
    int conn,
    int (*wr_log)(const char*, ...)
    )
{     
     // Size of local window (3 in current implementation).
    // Modify this value with 5,7,... if memory problem
    // occurs. Uncomment also some code (see further).
    int win_size;
         
    coords_queue_t queue;        // FIFO data structure for coords storage
    coords_t coords;            // Current coords
 
    int ct, offset;       // Counter for local subvolume size
    int i,j,k;
    int a,b,c;
 
    // Dimensions for padded/cropped volumes:
    int a_dimx, a_dimy, a_dimz;
    int a_rad;
 
    // Padded and cropped temporary input and output:
    unsigned char* tmp_in_rev;
    unsigned char* tmp_out_rev;
 
    // Pointer to a function for kind of connectivity:
    void (*conn_fun) ( 
        unsigned char*, 
        int,  
        int,  
        int,  
        int, 
        int, 
        int, 
        int,   
        int,   
        unsigned short, 
        coords_queue_t*            // FIFO data structure for coords storage
        );
 
    // The win_size parameters controls the tradeoff between computational time
    // and memory occupation. The fastest execution requires win_size = 3.
    win_size = 3;
 
 
    // Set the correct type of connectivity based on input parameter:
    if ( conn == CONN6 )
        conn_fun = _p3dGetRegionByCoords3D_conn_fun6;
    else if ( conn == CONN18 )
        conn_fun = _p3dGetRegionByCoords3D_conn_fun18;
    else // default:
        conn_fun = _p3dGetRegionByCoords3D_conn_fun26;
 
    /*char auth_code;

    //
    // Authenticate:
    //
    auth_code = authenticate("p3dGetRegionByCoords3D");
    if (auth_code == '0') goto AUTH_ERROR;*/
 
      // Start tracking computational time:
    if (wr_log != NULL) {
        p3dResetStartTime();
        wr_log("Pore3D - Extracting blob including voxel [%d,%d,%d]...",x,y,z);
        if (conn == CONN6)
            wr_log("\t6-connectivity used. ");
        else if (conn == CONN18)
            wr_log("\t18-connectivity used. ");
        else // default:
            wr_log("\t26-connectivity used. ");
    }
    
 
    // Apply algorithm:   
 
    // Create temporary input replicate-padded:
    a_rad = 1;
 
    // Compute dimensions of padded REV:
    a_dimx = dimx + a_rad*2;
    a_dimy = dimy + a_rad*2;
    a_dimz = dimz + a_rad*2;
 
    // Initialize input:
    tmp_in_rev = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) );
    if ( tmp_in_rev == NULL ) goto MEM_ERROR;
 
    p3dZeroPadding3D_8 ( in_rev, tmp_in_rev, dimx, dimy, dimz, a_rad, NULL, NULL );
 
    // Initialize output volume with a copy of input volume:
    tmp_out_rev = (unsigned char*) malloc( a_dimx*a_dimy*a_dimz*sizeof(unsigned char) );
    if ( tmp_out_rev == NULL ) goto MEM_ERROR;
 
    memcpy(tmp_out_rev,tmp_in_rev,a_dimx*a_dimy*a_dimz*sizeof(unsigned char));
 
 
    // Initialize variables:
    coords_queue_init(&queue);
 
 
    // Push the user specified object voxel into queue:
    coords.z = z; 
    coords.y = y;
    coords.x = x;
 
    coords_queue_push(&queue, coords);
 
    // While the queue is not empty:
    while(coords_queue_isempty(queue) == P3D_FALSE)
    {
        // Pop the first element of the queue:
        coords = coords_queue_pop(&queue);
         
        // Mark the extracted element:
        tmp_out_rev[ I(coords.x,coords.y,coords.z,a_dimx,a_dimy) ] = _GETREGIONBYCOORDS_TEMP_LABEL;
 
        // Perform sub-volume iterative scanning:
        ct = 1;
        while (ct < win_size)
        {
            offset = ct / 2; // integer division
 
            for (k = (coords.z - offset); k <= (coords.z + offset); k++)
                for (j = (coords.y - offset); j <= (coords.y + offset); j++)
                    for (i = (coords.x - offset); i <= (coords.x + offset); i++)
                    {             
                         
                        if (tmp_out_rev[ I(i,j,k,a_dimx,a_dimy) ] == _GETREGIONBYCOORDS_TEMP_LABEL)
                        {
                            // Perform 6-, 18-, or 26-connectivity:
                            conn_fun ( tmp_out_rev, a_dimx, a_dimy, a_dimx, i,
                                j, k, ct, win_size, BACKGROUND, &queue );
                             
                            // A voxels previously marked is set to background:
                            tmp_out_rev[ I(i,j,k,a_dimx,a_dimy) ] = _GETREGIONBYCOORDS_DEF_LABEL;     
                        }
                    }
            // Increase local sub-volume window:
            ct = ct + 2;
        }
    }     
 
 
    // Crop output:
    #pragma omp parallel for private(a,b)
    for( c = 0; c < dimz; c++ )
        for( b = 0; b < dimy; b++ )
            for( a = 0; a < dimx; a++ )            
            {  
                // Assign OBJECT value to output:
                out_rev[ I(a,b,c,dimx,dimy) ] = ( tmp_out_rev[ I( a + a_rad, b + a_rad, c + a_rad, a_dimx, a_dimy ) ] 
                    == _GETREGIONBYCOORDS_DEF_LABEL ) ? (unsigned char) (OBJECT) : (unsigned char) (BACKGROUND);                 
            }
 
     // Print out the elapsed time:
    if (wr_log != NULL) {
        wr_log("Pore3D - Region extracted successfully in %dm%0.3fs.", p3dGetElapsedTime_min(), p3dGetElapsedTime_sec());
    }
 
    // Free memory:
    if ( tmp_in_rev != NULL ) free(tmp_in_rev);
    if ( tmp_out_rev != NULL ) free(tmp_out_rev);
 
    // Return OK:
    return P3D_SUCCESS;
 
 
MEM_ERROR:

    // Log a ERROR message:
    if (wr_log != NULL) {
        wr_log("Pore3D - Not enough (contiguous) memory. Program will exit.");
    }



    // Free memory if previous malloc were successfully:
    if (tmp_in_rev != NULL) free(tmp_in_rev);
    if (tmp_out_rev != NULL) free(tmp_in_rev);

    // Return error code and exit:
    return P3D_MEM_ERROR;

/*AUTH_ERROR:

    if (wr_log != NULL) {
        wr_log("Pore3D - Authentication error: %s. Program will exit.", auth_code);
    }

    return P3D_AUTH_ERROR;*/
} 







