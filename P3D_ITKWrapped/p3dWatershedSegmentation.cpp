#include <itkImage.h>
#include <itkImportImageFilter.h>

#include <itkMorphologicalWatershedImageFilter.h>

#include <itkTimeProbe.h>

#include "_p3dCommon.h"
#include "p3dITKWrapped.h"

/**
 *
 * References
 * ----------
 * [1] F. Meyer. "Topographic distance and watershed lines". Signal Processing, 
 * 38:113-125, 1994.
 *
 ****************/


int _p3dMeyerWatershedSegmentation2D_8 (				// breadth-first version of the algorithm
	 unsigned char* in_im,
	 unsigned short* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 2;

	typedef unsigned char InPixelType;
	typedef unsigned short OutPixelType;
	typedef itk::Image< InPixelType, Dimension > InImageType;
	typedef itk::Image< OutPixelType, Dimension > OutImageType;

	typedef itk::ImportImageFilter< InPixelType, Dimension > ImportFilterType;

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
	
	InPixelType * pixelInData = static_cast< InPixelType * > (in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */


	// Set parameters for the specified filter:
	typedef itk::MorphologicalWatershedImageFilter< InImageType, OutImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	if ( conn == CONN4) 
		filter->SetFullyConnected(false); // 4-connectivity => false
	else 
		filter->SetFullyConnected(true); // 8-connectivity => true

	filter->SetMarkWatershedLine(true); // mark watersheds or not


	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	OutPixelType* pixelOutData = static_cast< OutPixelType* > (out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Apply filter:
		filter->Update();
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

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dMeyerWatershedSegmentation2D_16 (				// breadth-first version of the algorithm
	 unsigned short* in_im,
	 unsigned short* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 2;

	typedef unsigned short PixelType;
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */


	// Set parameters for the specified filter:
	typedef itk::MorphologicalWatershedImageFilter< ImageType, ImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	if ( conn == CONN4) 
		filter->SetFullyConnected(false); // 4-connectivity => false
	else 
		filter->SetFullyConnected(true); // 8-connectivity => true

	filter->SetMarkWatershedLine(true); // mark watersheds or not
    

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > (out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Apply filter:
		filter->Update();	
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


	// Return OK:
	return P3D_SUCCESS;
}

int p3dWatershedSegmentation2D_8 (				
	 unsigned char* in_im,
	 unsigned char* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	int ct, err_code;
	unsigned short* tmp_out_im = (unsigned short*) malloc(dimx*dimy*sizeof(unsigned short));

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Start tracking computational time:
	timer.Start();

	// Log a message:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Applying watershed segmentation 2D...";
	
		wr_log ( str.c_str() );
	}

	// Compute watersheds regions:
	err_code = _p3dMeyerWatershedSegmentation2D_8 ( in_im, tmp_out_im, dimx, dimy, conn, wr_log );

	// Extract watersheds:
	if ( err_code == P3D_SUCCESS )
	{
		#pragma omp parallel for private( ct )
		for ( ct = 0; ct < (dimx*dimy); ct++ )
		{
			out_im[ ct ] = (unsigned char) ( (tmp_out_im[ ct ] == 0) ? UCHAR_MAX : 0 );
		}

		// Release resources:
		free(tmp_out_im);
	}
	else
	{
		// Release resources:
		free(tmp_out_im);

		// Return ERROR:
		return P3D_ERROR;
	}

	// Stop timer:
	timer.Stop();

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Watershed segmentation 2D applied successfully in ";
		//str += p3dDoubleToString ( timer.GetMeanTime() );		
		//str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Return error code:
	return err_code;
}

int p3dWatershedSegmentation2D_16 (				
	 unsigned short* in_im,
	 unsigned char* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	int ct, err_code;
	unsigned short* tmp_out_im = (unsigned short*) malloc(dimx*dimy*sizeof(unsigned short));

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Start tracking computational time:
	timer.Start();

	// Log a message:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Applying watershed segmentation 2D...";
	
		wr_log ( str.c_str() );
	}

	// Compute watersheds regions:
	err_code = _p3dMeyerWatershedSegmentation2D_16 ( in_im, tmp_out_im, dimx, dimy, conn, wr_log );

	// Extract watersheds:
	if ( err_code == P3D_SUCCESS )
	{
		#pragma omp parallel for private( ct )
		for ( ct = 0; ct < (dimx*dimy); ct++ )
		{
			out_im[ ct ] = (unsigned char) ( (tmp_out_im[ ct ] == 0) ? UCHAR_MAX : 0 );
		}

		// Release resources:
		free(tmp_out_im);
	}
	else
	{
		// Release resources:
		free(tmp_out_im);

		// Return ERROR:
		return P3D_ERROR;
	}	

	// Stop timer:
	timer.Stop();

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Watershed segmentation 2D applied successfully in ";
		//str += p3dDoubleToString ( timer.GetMeanTime() );		
		//str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Release resources:
	free(tmp_out_im);

	// Return error code:
	return err_code;
}

int _p3dMeyerWatershedSegmentation3D_8 (				// breadth-first version of the algorithm
	 unsigned char* in_im,
	 unsigned short* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const unsigned int dimz,	
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 3;

	typedef unsigned char InPixelType;
	typedef unsigned short OutPixelType;
	typedef itk::Image< InPixelType, Dimension > InImageType;
	typedef itk::Image< OutPixelType, Dimension > OutImageType;

	typedef itk::ImportImageFilter< InPixelType, Dimension > ImportFilterType;

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
	
	InPixelType * pixelInData = static_cast< InPixelType * > (in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */


	// Set parameters for the specified filter:
	typedef itk::MorphologicalWatershedImageFilter< InImageType, OutImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	if ( conn == CONN6) 
		filter->SetFullyConnected(false); // 6-connectivity => false
	else 
		filter->SetFullyConnected(true); // 26-connectivity => true

	filter->SetMarkWatershedLine(true); // mark watersheds or not


	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	OutPixelType* pixelOutData = static_cast< OutPixelType* > (out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Apply filter:
		filter->Update();
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

	// Return OK:
	return P3D_SUCCESS;
}

int _p3dMeyerWatershedSegmentation3D_16 (				// breadth-first version of the algorithm
	 unsigned short* in_im,
	 unsigned short* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const unsigned int dimz,	
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	// Declare a C++ string stream used for log messages:
	std::string str;

	// Next, we select the data type to use to represent the image pixels.  We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int Dimension = 3;

	typedef unsigned short PixelType;
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_im );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
	 */


	// Set parameters for the specified filter:
	typedef itk::MorphologicalWatershedImageFilter< ImageType, ImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	
	filter->SetInput( importFilter->GetOutput() );

	if ( conn == CONN6) 
		filter->SetFullyConnected(false); // 6-connectivity => false
	else 
		filter->SetFullyConnected(true); // 26-connectivity => true

	filter->SetMarkWatershedLine(true); // mark watersheds or not
    

	/*
	 * End of section	 
	 **********************************************************************/


    // Return filtered output:
	PixelType* pixelOutData = static_cast< PixelType* > (out_im );
	filter->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	filter->GetOutput()->Allocate();


	// Execute filter catching exceptions:
	try
	{
		// Apply filter:
		filter->Update();	
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


	// Return OK:
	return P3D_SUCCESS;
}

int p3dWatershedSegmentation3D_8 (				
	 unsigned char* in_im,
	 unsigned char* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const unsigned int dimz,	
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	int ct, err_code;
	unsigned short* tmp_out_im = (unsigned short*) malloc(dimx*dimy*dimz*sizeof(unsigned short));

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Start tracking computational time:
	timer.Start();

	// Log a message:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Applying watershed segmentation 3D...";
	
		wr_log ( str.c_str() );
	}

	// Compute watersheds regions:
	err_code = _p3dMeyerWatershedSegmentation3D_8 ( in_im, tmp_out_im, dimx, dimy, dimz, conn, wr_log );

	// Extract watersheds:
	if ( err_code == P3D_SUCCESS )
	{
		#pragma omp parallel for private( ct )
		for ( ct = 0; ct < (dimx*dimy*dimz); ct++ )
		{
			out_im[ ct ] = (unsigned char) ( (tmp_out_im[ ct ] == 0) ? UCHAR_MAX : 0 );
		}

		// Release resources:
		free(tmp_out_im);
	}
	else
	{
		// Release resources:
		free(tmp_out_im);

		// Return ERROR:
		return P3D_ERROR;
	}

	// Stop timer:
	timer.Stop();

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Watershed segmentation 3D applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Return error code:
	return err_code;
}

int p3dWatershedSegmentation3D_16 (				
	 unsigned short* in_im,
	 unsigned char* out_im,
	 const unsigned int dimx,
	 const unsigned int dimy, 
	 const unsigned int dimz,	
	 const int conn,
	 int (*wr_log)(const char*, ...)
	 )
{
	int ct, err_code;
	unsigned short* tmp_out_im = (unsigned short*) malloc(dimx*dimy*dimz*sizeof(unsigned short));

	// Declare a timer:
	itk::TimeProbe timer;

	// Declare a C++ string stream used for log messages:
	std::string str;

	// Start tracking computational time:
	timer.Start();

	// Log a message:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Applying watershed segmentation 3D...";
	
		wr_log ( str.c_str() );
	}

	// Compute watersheds regions:
	err_code = _p3dMeyerWatershedSegmentation3D_16 ( in_im, tmp_out_im, dimx, dimy, dimz, conn, wr_log );

	// Extract watersheds:
	if ( err_code == P3D_SUCCESS )
	{
		#pragma omp parallel for private( ct )
		for ( ct = 0; ct < (dimx*dimy*dimz); ct++ )
		{
			out_im[ ct ] = (unsigned char) ( (tmp_out_im[ ct ] == 0) ? UCHAR_MAX : 0 );
		}

		// Release resources:
		free(tmp_out_im);
	}
	else
	{
		// Release resources:
		free(tmp_out_im);

		// Return ERROR:
		return P3D_ERROR;
	}	

	// Stop timer:
	timer.Stop();

	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Watershed segmentation 3D applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += ".";
		
		wr_log ( str.c_str() );
	}

	// Release resources:
	free(tmp_out_im);

	// Return error code:
	return err_code;
}