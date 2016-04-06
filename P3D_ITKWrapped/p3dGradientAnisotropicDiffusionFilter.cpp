#include <limits.h>

#include <itkImage.h>
#include <itkImportImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

#include <itkTimeProbe.h>

#include "_p3dCommon.h"
#include "p3dITKWrapped.h"


int p3dGradientAnisotropicDiffusionFilter2D_8(   
	unsigned char* in_im, 
	unsigned char* out_im, 
	const int dimx,
	const int dimy, 
	const int iterations,
	const double conductance,
	int (*wr_log)(const char*, ...)
	)
{	

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 2;

	typedef unsigned char PixelType;
	typedef float        RealPixelType;  
	
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::Image< RealPixelType, Dimension > RealImageType;

	typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;

	ImportFilterType::Pointer importFilter = ImportFilterType::New();      

	// This filter requires the user to specify the size of the image to be
	// produced as output. The image size should exactly match the number of 
	// pixels available in the locally allocated buffer. 
	ImportFilterType::SizeType  size;

	size[0]  = dimx;  // size along X
	size[1]  = dimy;  // size along Y

	ImportFilterType::IndexType start;
	start.Fill( 0 );

	ImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );

	importFilter->SetRegion( region );

	double origin[ Dimension ];
	origin[0] = 0.0;    // X coordinate 
	origin[1] = 0.0;    // Y coordinate

	importFilter->SetOrigin( origin );

	double spacing[ Dimension ];
	spacing[0] = 1.0;    // along X direction 
	spacing[1] = 1.0;    // along Y direction

	importFilter->SetSpacing( spacing );

	const unsigned long numberOfPixels =  size[0] * size[1];
	
	PixelType * pixelInData = static_cast< PixelType * > ( in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// This filter requires casting to a real type_
	typedef itk::CastImageFilter< ImageType, RealImageType> CastToRealFilterType;
	CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();

	// Set parameters for the specified filter:
	typedef itk::GradientAnisotropicDiffusionImageFilter< RealImageType, RealImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();

	// The input parameter is the input of the "filter" that performs casting:
	toReal->SetInput( importFilter->GetOutput() );

	// The input of the filter is the output of the "filter" that performs casting:
	filter->SetInput( toReal->GetOutput() );

	// Setting filter parameters:
	filter->SetNumberOfIterations( iterations );
	filter->SetTimeStep( 0.125 ); // Minimum stable step
	filter->SetConductanceParameter( conductance );

	// This filter requires also rescaling:
	typedef itk::RescaleIntensityImageFilter< 
		RealImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	// 8-bit range:
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( UCHAR_MAX );

	// The input of the rescaler is the output of the filter:
	rescaler->SetInput( filter->GetOutput() );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return rescaled filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > ( out_im );
	rescaler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	rescaler->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying gradient anisotropic diffusion filter 2D...";		
                        wr_log ( str.c_str() );
                        str  = "\tK: ";
                        str += _p3dDoubleToString ( conductance );	
                        str += ".";		
        		wr_log ( str.c_str() );
                        str  = "\tIterations: ";
                        str += _p3dIntToString ( iterations );	
                        str += ".";  
                        wr_log ( str.c_str() );
		}

		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		rescaler->Update();

		// Stop timer:
		timer.Stop();
	}
	catch( itk::ExceptionObject & exp ) 
	{
		// Print exception:			
		if (wr_log != NULL)
		{
			str  = "Pore3D - Error: ";
			str += exp.GetDescription();
			str += " .";

			wr_log ( str.c_str() );
		}
		
		
		// Return ERROR:
		return P3D_ERROR;
	}

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Gradient anisotropic diffusion filter 2D applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );			
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}


int p3dGradientAnisotropicDiffusionFilter2D_16(   
	unsigned short* in_im, 
	unsigned short* out_im, 
	const int dimx,
	const int dimy, 
	const int iterations,
	const double conductance,
	int (*wr_log)(const char*, ...)
	)
{	

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 2;

	typedef unsigned short PixelType;
	typedef float        RealPixelType;  
	
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::Image< RealPixelType, Dimension > RealImageType;

	typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;

	ImportFilterType::Pointer importFilter = ImportFilterType::New();      

	// This filter requires the user to specify the size of the image to be
	// produced as output. The image size should exactly match the number of 
	// pixels available in the locally allocated buffer. 
	ImportFilterType::SizeType  size;

	size[0]  = dimx;  // size along X
	size[1]  = dimy;  // size along Y

	ImportFilterType::IndexType start;
	start.Fill( 0 );

	ImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );

	importFilter->SetRegion( region );

	double origin[ Dimension ];
	origin[0] = 0.0;    // X coordinate 
	origin[1] = 0.0;    // Y coordinate

	importFilter->SetOrigin( origin );

	double spacing[ Dimension ];
	spacing[0] = 1.0;    // along X direction 
	spacing[1] = 1.0;    // along Y direction

	importFilter->SetSpacing( spacing );

	const unsigned long numberOfPixels =  size[0] * size[1];
	
	PixelType * pixelInData = static_cast< PixelType * > ( in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// This filter requires casting to a real type_
	typedef itk::CastImageFilter< ImageType, RealImageType> CastToRealFilterType;
	CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();

	// Set parameters for the specified filter:
	typedef itk::GradientAnisotropicDiffusionImageFilter< RealImageType, RealImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();

	// The input parameter is the input of the "filter" that performs casting:
	toReal->SetInput( importFilter->GetOutput() );

	// The input of the filter is the output of the "filter" that performs casting:
	filter->SetInput( toReal->GetOutput() );

	// Setting filter parameters:
	filter->SetNumberOfIterations( iterations );
	filter->SetTimeStep( 0.125 ); // Minimum stable step
	filter->SetConductanceParameter( conductance );

	// This filter requires also rescaling:
	typedef itk::RescaleIntensityImageFilter< 
		RealImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	// 8-bit range:
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( USHRT_MAX );

	// The input of the rescaler is the output of the filter:
	rescaler->SetInput( filter->GetOutput() );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return rescaled filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > ( out_im );
	rescaler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	rescaler->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying gradient anisotropic diffusion filter...";
                        wr_log ( str.c_str() );
                        str  = "\tK: ";
                        str += _p3dDoubleToString ( conductance );	
                        str += ".";		
        		wr_log ( str.c_str() );
                        str  = "\tIterations: ";
                        str += _p3dIntToString ( iterations );	
                        str += ".";  
                        wr_log ( str.c_str() );			
		}

		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		rescaler->Update();

		// Stop timer:
		timer.Stop();
	}
	catch( itk::ExceptionObject & exp ) 
	{
		// Print exception:			
		if (wr_log != NULL)
		{
			str  = "Pore3D - Error: ";
			str += exp.GetDescription();
			str += " .";

			wr_log ( str.c_str() );
		}
		
		
		// Return ERROR:
		return P3D_ERROR;
	}

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Gradient anisotropic diffusion filter applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );			
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dGradientAnisotropicDiffusionFilter3D_8(   
	unsigned char* in_rev, 
	unsigned char* out_rev, 
	const int dimx,
	const int dimy, 
	const int dimz,
	const int numberOfIterations,
	const double conductance,
	int (*wr_log)(const char*, ...)
	)
{	

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 3;

	typedef unsigned char PixelType;
	typedef float        RealPixelType;  
	
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::Image< RealPixelType, Dimension > RealImageType;

	typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;

	ImportFilterType::Pointer importFilter = ImportFilterType::New();      

	// This filter requires the user to specify the size of the image to be
	// produced as output. The image size should exactly match the number of 
	// pixels available in the locally allocated buffer. 
	ImportFilterType::SizeType  size;

	size[0]  = dimx;  // size along X
	size[1]  = dimy;  // size along Y
	size[2]  = dimz;  // size along Z

	ImportFilterType::IndexType start;
	start.Fill( 0 );

	ImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );

	importFilter->SetRegion( region );

	double origin[ Dimension ];
	origin[0] = 0.0;    // X coordinate 
	origin[1] = 0.0;    // Y coordinate
	origin[2] = 0.0;    // Z coordinate

	importFilter->SetOrigin( origin );

	double spacing[ Dimension ];
	spacing[0] = 1.0;    // along X direction 
	spacing[1] = 1.0;    // along Y direction
	spacing[2] = 1.0;    // along Z direction

	importFilter->SetSpacing( spacing );

	const unsigned long numberOfPixels =  size[0] * size[1] * size[2];
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// This filter requires casting to a real type_
	typedef itk::CastImageFilter< ImageType, RealImageType> CastToRealFilterType;
	CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();

	// Set parameters for the specified filter:
	typedef itk::GradientAnisotropicDiffusionImageFilter< RealImageType, RealImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();

	// The input parameter is the input of the "filter" that performs casting:
	toReal->SetInput( importFilter->GetOutput() );

	// The input of the filter is the output of the "filter" that performs casting:
	filter->SetInput( toReal->GetOutput() );

	// Setting filter parameters:
	filter->SetNumberOfIterations( numberOfIterations );
	filter->SetTimeStep( 0.0625 ); // Minimum stable step
	filter->SetConductanceParameter( conductance );

	// This filter requires also rescaling:
	typedef itk::RescaleIntensityImageFilter< 
		RealImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	// 8-bit range:
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( UCHAR_MAX );

	// The input of the rescaler is the output of the filter:
	rescaler->SetInput( filter->GetOutput() );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return rescaled filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > (out_rev );
	rescaler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	rescaler->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying gradient anisotropic diffusion filter...";		
			wr_log ( str.c_str() );
                        str  = "\tK: ";
                        str += _p3dDoubleToString ( conductance );	
                        str += ".";		
        		wr_log ( str.c_str() );
                        str  = "\tIterations: ";
                        str += _p3dIntToString ( numberOfIterations );	
                        str += ".";  
                        wr_log ( str.c_str() );
		}
                

		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		rescaler->Update();

		// Stop timer:
		timer.Stop();
	}
	catch( itk::ExceptionObject & exp ) 
	{
		// Print exception:			
		if (wr_log != NULL)
		{
			str  = "Pore3D - Error: ";
			str += exp.GetDescription();
			str += " .";

			wr_log ( str.c_str() );
		}
		
		
		// Return ERROR:
		return P3D_ERROR;
	}

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Gradient anisotropic diffusion filter applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );			
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dGradientAnisotropicDiffusionFilter3D_16(   
	unsigned short* in_rev, 
	unsigned short* out_rev, 
	const int dimx,
	const int dimy, 
	const int dimz,
	const int numberOfIterations,
	const double conductance,
	int (*wr_log)(const char*, ...)
	)
{	

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 3;

	typedef unsigned short PixelType;
	typedef float        RealPixelType;  
	
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::Image< RealPixelType, Dimension > RealImageType;

	typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;

	ImportFilterType::Pointer importFilter = ImportFilterType::New();      

	// This filter requires the user to specify the size of the image to be
	// produced as output. The image size should exactly match the number of 
	// pixels available in the locally allocated buffer. 
	ImportFilterType::SizeType  size;

	size[0]  = dimx;  // size along X
	size[1]  = dimy;  // size along Y
	size[2]  = dimz;  // size along Z

	ImportFilterType::IndexType start;
	start.Fill( 0 );

	ImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );

	importFilter->SetRegion( region );

	double origin[ Dimension ];
	origin[0] = 0.0;    // X coordinate 
	origin[1] = 0.0;    // Y coordinate
	origin[2] = 0.0;    // Z coordinate

	importFilter->SetOrigin( origin );

	double spacing[ Dimension ];
	spacing[0] = 1.0;    // along X direction 
	spacing[1] = 1.0;    // along Y direction
	spacing[2] = 1.0;    // along Z direction

	importFilter->SetSpacing( spacing );

	const unsigned long numberOfPixels =  size[0] * size[1] * size[2];
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// This filter requires casting to a real type_
	typedef itk::CastImageFilter< ImageType, RealImageType> CastToRealFilterType;
	CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();

	// Set parameters for the specified filter:
	typedef itk::GradientAnisotropicDiffusionImageFilter< RealImageType, RealImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();

	// The input parameter is the input of the "filter" that performs casting:
	toReal->SetInput( importFilter->GetOutput() );

	// The input of the filter is the output of the "filter" that performs casting:
	filter->SetInput( toReal->GetOutput() );

	// Setting filter parameters:
	filter->SetNumberOfIterations( numberOfIterations );
	filter->SetTimeStep( 0.0625 ); // Minimum stable step
	filter->SetConductanceParameter( conductance );

	// This filter requires also rescaling:
	typedef itk::RescaleIntensityImageFilter< 
		RealImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	// 16-bit range:
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( USHRT_MAX );

	// The input of the rescaler is the output of the filter:
	rescaler->SetInput( filter->GetOutput() );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return rescaled filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > (out_rev );
	rescaler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	rescaler->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying gradient anisotropic diffusion filter...";		
			wr_log ( str.c_str() );
                        str  = "\tK: ";
                        str += _p3dDoubleToString ( conductance );	
                        str += ".";		
        		wr_log ( str.c_str() );
                        str  = "\tIterations: ";
                        str += _p3dIntToString ( numberOfIterations );	
                        str += ".";  
                        wr_log ( str.c_str() );
		}

		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		rescaler->Update();

		// Stop timer:
		timer.Stop();
	}
	catch( itk::ExceptionObject & exp ) 
	{
		// Print exception:			
		if (wr_log != NULL)
		{
			str  = "Pore3D - Error: ";
			str += exp.GetDescription();
			str += " .";

			wr_log ( str.c_str() );
		}
		
		
		// Return ERROR:
		return P3D_ERROR;
	}

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Gradient anisotropic diffusion filter applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );			
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}
