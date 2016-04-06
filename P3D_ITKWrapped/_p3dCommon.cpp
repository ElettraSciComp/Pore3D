#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "_p3dCommon.h"


std::string _p3dTimeToString ( double sec )
{
        int    min;   
	std::ostringstream oss;

        // Return time in seconds:	
        min = (int) (sec / 60);
        sec = sec - min*60;                
        
	if (!(oss << min << "m" << std::fixed << std::setprecision(3) << sec << "s"))
		return "";   
	return oss.str();
 } 

std::string _p3dIntToString ( int value )
{
	std::ostringstream oss;

	if (!(oss << value))
		return "";   
	return oss.str();
} 

std::string _p3dDoubleToString ( double value )
{
	std::ostringstream oss;

	if (!(oss << value))
		return "";   
	return oss.str();
 }