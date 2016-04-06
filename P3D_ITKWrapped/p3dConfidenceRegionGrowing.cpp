#include <limits.h>

#include <itkImage.h>
#include <itkImportImageFilter.h>

#include <itkConfidenceConnectedImageFilter.h>

#include <itkTimeProbe.h>

#include "_p3dCommon.h"
#include "p3dITKWrapped.h"


int p3dConfidenceConnectedRegionGrowing2D_8 (   
	unsigned char* in_im, 
	unsigned char* out_im, 
	const unsigned int dimx,
	const unsigned int dimy, 
	const double multiplier,
	const unsigned int winsize,
	const unsigned int iterations,
	const unsigned int seed_x,
	const unsigned int seed_y,
	int (*wr_log)( const char*, ... )
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
	typedef itk::Image< PixelType, Dimension > ImageType;

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

	// Set parameters for the specified filter:
	typedef itk::ConfidenceConnectedImageFilter<ImageType, ImageType>   FilterType;	
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	filter->SetMultiplier( multiplier );
	filter->SetNumberOfIterations( iterations );

	ImageType::IndexType  index;
  
	index[0] = seed_x;
	index[1] = seed_y;

	filter->SetSeed( index );

	filter->SetReplaceValue( UCHAR_MAX );
	// Set the radius of the neighborhood over which the statistics 
	// are evaluated:
	filter->SetInitialNeighborhoodRadius( winsize / 2 );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	PixelType * pixelOutData = static_cast< PixelType * > ( out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying confidence connected segmentation 2D...";		
			wr_log ( str.c_str() );
                        
                        str  = "\t Seed point: [";
			str += _p3dIntToString ( (int) seed_x );		
			str += ",";
                        str += _p3dIntToString ( (int) seed_y );		
                        str += "].";
                        str  = "\t Width of seed neighborhood: ";
			str += _p3dIntToString ( (int) winsize );		
			str += ".";
                        str  = "\t K: ";
			str += _p3dDoubleToString ( (int) multiplier );		
			str += ".";                        
                        str  = "\t Iterations: ";
			str += _p3dIntToString ( (int) iterations );		
			str += ".";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		filter->Update();

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
		str  = "Pore3D - Confidence connected segmentation 2D applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}



int p3dConfidenceConnectedRegionGrowing2D_16 (   
	unsigned short* in_im, 
	unsigned char* out_im, 
	const unsigned int dimx,
	const unsigned int dimy, 
	const double multiplier,
	const unsigned int winsize,
	const unsigned int iterations,
	const unsigned int seed_x,
	const unsigned int seed_y,
	int (*wr_log)( const char*, ... )
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
	typedef itk::Image< PixelType, Dimension > ImageType;

	typedef unsigned char OutPixelType;
	typedef itk::Image< OutPixelType, Dimension > OutImageType;


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

	// Set parameters for the specified filter:
	typedef itk::ConfidenceConnectedImageFilter<ImageType, OutImageType>   FilterType;	
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	filter->SetMultiplier( multiplier );
	filter->SetNumberOfIterations( iterations );

	ImageType::IndexType  index;
  
	index[0] = seed_x;
	index[1] = seed_y;

	filter->SetSeed( index );

	filter->SetReplaceValue( UCHAR_MAX );
	filter->SetInitialNeighborhoodRadius( winsize / 2 );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	OutPixelType * pixelOutData = static_cast< OutPixelType * > ( out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying confidence connected segmentation 2D...";		
			wr_log ( str.c_str() );
                        
                        str  = "\t Seed point: [";
			str += _p3dIntToString ( (int) seed_x );		
			str += ",";
                        str += _p3dIntToString ( (int) seed_y );		
                        str += "].";
                        str  = "\t Width of seed neighborhood: ";
			str += _p3dIntToString ( (int) winsize );		
			str += ".";
                        str  = "\t K: ";
			str += _p3dDoubleToString ( (int) multiplier );		
			str += ".";                        
                        str  = "\t Iterations: ";
			str += _p3dIntToString ( (int) iterations );		
			str += ".";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		filter->Update();

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
		str  = "Pore3D - Confidence connected segmentation 2D applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dConfidenceConnectedRegionGrowing3D_8 (   
	unsigned char* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy,
	const unsigned int dimz,
	const double multiplier,
	const unsigned int winsize,
	const unsigned int iterations,
	const unsigned int seed_x,
	const unsigned int seed_y,
	const unsigned int seed_z,
	int (*wr_log)( const char*, ... )
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
	typedef itk::Image< PixelType, Dimension > ImageType;

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
	
	PixelType * pixelInData = static_cast< PixelType * > ( in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// Set parameters for the specified filter:
	typedef itk::ConfidenceConnectedImageFilter<ImageType, ImageType>   FilterType;	
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	filter->SetMultiplier( multiplier );
	filter->SetNumberOfIterations( iterations );

	ImageType::IndexType  index;
  
	index[0] = seed_x;
	index[1] = seed_y;
	index[2] = seed_z;

	filter->SetSeed( index );

	filter->SetReplaceValue( UCHAR_MAX );
	filter->SetInitialNeighborhoodRadius( winsize / 2 );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	PixelType * pixelOutData = static_cast< PixelType * > ( out_rev );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying confidence interval region growing...";		
			wr_log ( str.c_str() );
                        
                        str  = "\t Seed point: [";
			str += _p3dIntToString ( (int) seed_x );		
			str += ",";
                        str += _p3dIntToString ( (int) seed_y );		
                        str += ",";
                        str += _p3dIntToString ( (int) seed_z );		
                        str += "].";
                        str  = "\t Width of seed neighborhood: ";
			str += _p3dIntToString ( (int) winsize );		
			str += ".";
                        str  = "\t K: ";
			str += _p3dDoubleToString ( (int) multiplier );		
			str += ".";                        
                        str  = "\t Iterations: ";
			str += _p3dIntToString ( (int) iterations );		
			str += ".";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		filter->Update();

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
		str  = "Pore3D - Confidence interval region growing applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}

int p3dConfidenceConnectedRegionGrowing3D_16 (   
	unsigned short* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy,
	const unsigned int dimz,
	const double multiplier,
	const unsigned int winsize,
	const unsigned int iterations,
	const unsigned int seed_x,
	const unsigned int seed_y,
	const unsigned int seed_z,
	int (*wr_log)( const char*, ... )
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
	typedef itk::Image< PixelType, Dimension > ImageType;
	
	typedef unsigned char OutPixelType;
	typedef itk::Image< OutPixelType, Dimension > OutImageType;

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
	
	PixelType * pixelInData = static_cast< PixelType * > ( in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */

	// Set parameters for the specified filter:
	typedef itk::ConfidenceConnectedImageFilter<ImageType, OutImageType>   FilterType;	
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	filter->SetMultiplier( multiplier );
	filter->SetNumberOfIterations( iterations );

	ImageType::IndexType  index;
  
	index[0] = seed_x;
	index[1] = seed_y;
	index[2] = seed_z;

	filter->SetSeed( index );

	filter->SetReplaceValue( UCHAR_MAX );
	filter->SetInitialNeighborhoodRadius( winsize / 2 );

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	OutPixelType * pixelOutData = static_cast< OutPixelType * > ( out_rev );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying confidence interval region growing...";		
			wr_log ( str.c_str() );
                        
                        str  = "\t Seed point: [";
			str += _p3dIntToString ( (int) seed_x );		
			str += ",";
                        str += _p3dIntToString ( (int) seed_y );		
                        str += ",";
                        str += _p3dIntToString ( (int) seed_z );		
                        str += "].";
                        str  = "\t Width of seed neighborhood: ";
			str += _p3dIntToString ( (int) winsize );		
			str += ".";
                        str  = "\t K: ";
			str += _p3dDoubleToString ( (int) multiplier );		
			str += ".";                        
                        str  = "\t Iterations: ";
			str += _p3dIntToString ( (int) iterations );		
			str += ".";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		// Apply filter:
		filter->Update();

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
		str  = "Pore3D - Confidence interval region growing applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}







