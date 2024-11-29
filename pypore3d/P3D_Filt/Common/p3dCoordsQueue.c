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

#include "../p3dFilt.h"	
#include "p3dCoordsQueue.h"

void coords_queue_init(coords_queue_t *queue) {
    queue->head = NULL;
    queue->tail = NULL;
}

void coords_queue_push(coords_queue_t *queue, coords_t elem) {

    coords_queue_elem_t *new_q;

    // Alloc memory for the new item:
    new_q = (coords_queue_elem_t *) malloc(sizeof (coords_queue_elem_t));

    // Push item into queue:
    new_q->item = elem;
    new_q->next = NULL;

    // Handle first element pushed:
    if (queue->tail != NULL) {
        queue->tail->next = new_q;
        queue->tail = queue->tail->next;
    } else {
        queue->tail = new_q;
        queue->head = queue->tail;
    }
}

coords_t coords_queue_pop(coords_queue_t *queue) {
    coords_t elem;
    coords_queue_elem_t *temp;

    // Pop item from queue:
    elem = queue->head->item;

    // Free memory:
    temp = queue->head;
    queue->head = queue->head->next;
    free(temp);

    // Handle last element popped:
    if (queue->head == NULL) {
        queue->tail = NULL;
    }

    // Return the item previously popped:
    return elem;
}

int coords_queue_isempty(coords_queue_t queue) {
    return ( ((queue.head == NULL) && (queue.tail == NULL)) ? P3D_TRUE : P3D_FALSE);
}