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

#include <stdio.h>

#include "p3dTime.h"

int    p3dTime_min = 0;
double p3dTime_sec = 0.0;

#ifndef _WINDOWS
        #include <sys/time.h>

	static struct timeval p3dTime_tv;
	static unsigned long p3dTime_startTime, p3dTime_crtTime;
#else
        #include <windows.h>
	#include <time.h>

	DWORD	p3dTime_startTime, p3dTime_crtTime;	
#endif


void p3dResetStartTime() 
{
	#ifdef _WINDOWS
		p3dTime_startTime = GetTickCount();
		p3dTime_crtTime = GetTickCount();
	#else
		gettimeofday(&p3dTime_tv, NULL);
		p3dTime_startTime = (p3dTime_tv.tv_sec * 1000) + (p3dTime_tv.tv_usec / 1000);
	#endif	

	p3dTime_crtTime = p3dTime_startTime;
}

double p3dGetElapsedTime () 
{    
	#ifdef _WINDOWS
		p3dTime_crtTime = GetTickCount();
	#else
		gettimeofday(&p3dTime_tv, NULL);
		p3dTime_crtTime = (p3dTime_tv.tv_sec * 1000) + (p3dTime_tv.tv_usec / 1000);
	#endif
	
	// Return time in seconds:
	return ((p3dTime_crtTime - p3dTime_startTime)/1000.0);	             
}

int p3dGetElapsedTime_min () 
{    
	#ifdef _WINDOWS
		p3dTime_crtTime = GetTickCount();
	#else
		gettimeofday(&p3dTime_tv, NULL);
		p3dTime_crtTime = (p3dTime_tv.tv_sec * 1000) + (p3dTime_tv.tv_usec / 1000);
	#endif
	
	// Return time in seconds:
	p3dTime_sec = ((p3dTime_crtTime - p3dTime_startTime)/1000.0);	
        p3dTime_min = (int) (p3dTime_sec / 60.0);
        p3dTime_sec = p3dTime_sec - p3dTime_min*60.0;
        
        return p3dTime_min;        
}

double p3dGetElapsedTime_sec () 
{           
        return p3dTime_sec;        
}
