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


/********************************************************************************
 * File:			p3dCoords_list.h
 *
 * Description:		Header file for the implementation of the dynamic structure
 *					with LIFO (Last-In-First-Out) policy (list or stack)
 *					containing CoordsT elements (see p3dCoordsT.h)
 *					Duplicate elements allowed.
 *
 * Interfaces:		coords_list_init
 *					coords_list_push
 *					coords_list_pop
 *					coords_list_isempty
 *					coords_list_toarray
 *
 * Author:			FB
 *
 * Last Modified:
 *
 * Copyright (C) 2009 Sincrotrone Trieste S.C.p.A. - All rights reserved.
 ********************************************************************************/
#include "p3dCoordsT.h"

/*
	Type definitions:
*/ 

#ifndef COORDS_LIST_T_DEFINED
	#define COORDS_LIST_T_DEFINED

	struct coords_lelem_t {
		coords_t elem;
		struct coords_lelem_t	*next;
	};

	typedef struct coords_lelem_t coords_list_elem_t;
		
	typedef coords_list_elem_t* coords_list_t;

#endif

/*
	Interfaces:
*/

/********************************************************************************
 * Function:		coords_list_init
 * 
 * Description:		Initialize the specified parameter to NULL.
 *
 * Input(s):		coords_list_t*		- The list to initialize
 *	
 * Output:			No return type
 ********************************************************************************/
void coords_list_init ( coords_list_t* );

/********************************************************************************
 * Function:		coords_list_push
 * 
 * Description:		Push the specified element into the list.
 *
 * Input(s):		coords_list_t*		- The list to extend
 *					coords_t			- The element to push into the list
 *	
 * Output:			P3D_SUCCESS if element successfully pushed into the list
 *					P3D_MEM_ERROR if there is not enough memory for the additional 
 *                  element
 ********************************************************************************/
int coords_list_push ( coords_list_t*, coords_t );

/********************************************************************************
 * Function:		coords_list_pop
 * 
 * Description:		Pop an element from the list according to the LIFO policy.
 *
 * Input(s):		coords_list_t*		- The list to reduce
 *	
 * Output:			coords_t			- The popped element
 ********************************************************************************/
coords_t coords_list_pop ( coords_list_t* );

/********************************************************************************
 * Function:		coords_list_isempty
 * 
 * Description:		Return P3D_TRUE if the specified list is empty
 *
 * Input(s):		coords_list_t*		- The list to query
 *	
 * Output:			P3D_TRUE if the list is empty
 *					P3D_FALSE if there is at least one element into the list
 ********************************************************************************/
int coords_list_isempty ( coords_list_t );

/********************************************************************************
 * Function:		coords_list_toarray
 * 
 * Description:		Convert the dynamic structure to a static array. The length
 *					of the array should be known a-priori and specified in input.
 *					List is deleted after this operation. Run-time error occurs if 
 *					the caller specify a number of elements greater than the real
 *					value. Memory leak occurs if the caller specify a number of 
 *					elements lower than the real value.
 *
 * Input(s):		coords_list_t*		- The list to convert
 *					int					- The number of list elements 
 *	
 * Output:			The pointer to the array or NULL if there is not enough 
 *					(contiguous) memory for the array (list is deleted in any 
 *					case).
 ********************************************************************************/
coords_t* coords_list_toarray ( coords_list_t*, int );


