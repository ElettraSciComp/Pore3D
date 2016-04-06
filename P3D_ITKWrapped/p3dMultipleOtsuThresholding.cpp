#include <limits.h>

#include <itkImage.h>
#include <itkImportImageFilter.h>

#include <itkThresholdLabelerImageFilter.h> 
#include <itkScalarImageToHistogramGenerator.h>
#include <itkOtsuMultipleThresholdsCalculator.h>

#include <itkNumericTraits.h>

#include <itkTimeProbe.h>

#include "_p3dCommon.h"
#include "p3dITKWrapped.h"

// The number of thresholds is the number of regions - 1

int p3dMultipleOtsuThresholding2D_8 (   
	unsigned char* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy, 
	const unsigned int regions,
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
     */

	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > 
		ScalarImageToHistogramGeneratorType;

	typedef ScalarImageToHistogramGeneratorType::HistogramType    HistogramType;

	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType >   CalculatorType;

	typedef itk::ThresholdLabelerImageFilter< ImageType, ImageType > ThresholdLabelerType;

	
	
	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = 
		ScalarImageToHistogramGeneratorType::New();

	CalculatorType::Pointer calculator = CalculatorType::New();

	ThresholdLabelerType::Pointer labeler = ThresholdLabelerType::New();
	
	
	// Parameters settings:
	scalarImageToHistogramGenerator->SetNumberOfBins( (UCHAR_MAX + 1) / 2 );
	calculator->SetNumberOfThresholds( regions - 1 );	

	// Setting pipeline:	
	scalarImageToHistogramGenerator->SetInput( importFilter->GetOutput() );
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
  
	
	// Return filtered output:
	PixelType * pixelOutData = static_cast< PixelType * > (out_rev );
	labeler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	labeler->GetOutput()->Allocate();

	/*
	 *
	 *****************************************************************/

	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying multiple Otsu's thresholding...";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		//Invoke pipeline manually in order:
	    importFilter->Update();
  
		scalarImageToHistogramGenerator->Compute();
  
		calculator->Update();

		labeler->SetInput( importFilter->GetOutput() );
	    labeler->SetRealThresholds( calculator->GetOutput()) ;
    
		labeler->Update();

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

	  
	// Print out the computed set of thresholds:

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	//const CalculatorType::OutputType &thresholdVector = labeler->GetRealThresholds(); 
	CalculatorType::OutputType::const_iterator itNum;

	if (wr_log != NULL)
	{
		for(itNum = thresholdVector.begin(); itNum < thresholdVector.end(); itNum++) 
		{			
			str  = "\tPore3D - Threshold  ";
			str += _p3dIntToString( (unsigned int) (itNum - thresholdVector.begin()) ) ;
			str += " proposed by Otsu's method is: ";
			str += _p3dIntToString ( static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>(*itNum) );
			str += " .";	
		
			wr_log ( str.c_str() );
		}  
    }


	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Multiple Otsu's thresholding applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}









int p3dMultipleOtsuThresholding2D_16 (   
	unsigned short* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy, 	
	const unsigned int numberOfRegions,
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
     */

	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > 
		ScalarImageToHistogramGeneratorType;

	typedef ScalarImageToHistogramGeneratorType::HistogramType    HistogramType;

	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType >   CalculatorType;

	typedef itk::ThresholdLabelerImageFilter< ImageType, OutImageType > ThresholdLabelerType;

	
	
	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = 
		ScalarImageToHistogramGeneratorType::New();

	CalculatorType::Pointer calculator = CalculatorType::New();

	ThresholdLabelerType::Pointer labeler = ThresholdLabelerType::New();
	
	
	// Parameters settings:
	scalarImageToHistogramGenerator->SetNumberOfBins( (UCHAR_MAX + 1) / 2 );
	calculator->SetNumberOfThresholds( numberOfRegions - 1);	

	// Setting pipeline:	
	scalarImageToHistogramGenerator->SetInput( importFilter->GetOutput() );
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
  
	
	// Return filtered output:
	OutPixelType * pixelOutData = static_cast< OutPixelType * > ( out_rev );
	labeler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	labeler->GetOutput()->Allocate();

	/*
	 *
	 *****************************************************************/

	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying multiple Otsu's thresholding...";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		//Invoke pipeline manually in order:
	    importFilter->Update();
  
		scalarImageToHistogramGenerator->Compute();
  
		calculator->Update();

		labeler->SetInput( importFilter->GetOutput() );
	    labeler->SetRealThresholds( calculator->GetOutput()) ;
    
		labeler->Update();

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

	  
	// Print out the computed set of thresholds:

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	//const CalculatorType::OutputType &thresholdVector = labeler->GetRealThresholds(); 
	CalculatorType::OutputType::const_iterator itNum;

	if (wr_log != NULL)
	{
		for(itNum = thresholdVector.begin(); itNum < thresholdVector.end(); itNum++) 
		{			
			str  = "\tPore3D - Threshold  ";
			str += _p3dIntToString( (unsigned int) (itNum - thresholdVector.begin()) ) ;
			str += " proposed by Otsu's method is: ";
			str += _p3dIntToString ( static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>(*itNum) );
			str += " .";	
		
			wr_log ( str.c_str() );
		}
    }


	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Multiple Otsu's thresholding applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}


int p3dMultipleOtsuThresholding3D_8 (   
	unsigned char* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy, 
	const unsigned int dimz,
	const unsigned int numberOfRegions,
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
     */

	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > 
		ScalarImageToHistogramGeneratorType;

	typedef ScalarImageToHistogramGeneratorType::HistogramType    HistogramType;

	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType >   CalculatorType;

	typedef itk::ThresholdLabelerImageFilter< ImageType, ImageType > ThresholdLabelerType;

	
	
	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = 
		ScalarImageToHistogramGeneratorType::New();

	CalculatorType::Pointer calculator = CalculatorType::New();

	ThresholdLabelerType::Pointer labeler = ThresholdLabelerType::New();
	
	
	// Parameters settings:
	scalarImageToHistogramGenerator->SetNumberOfBins( (UCHAR_MAX + 1) / 2 );
	calculator->SetNumberOfThresholds( numberOfRegions - 1);
        //calculator->

	// Setting pipeline:	
	scalarImageToHistogramGenerator->SetInput( importFilter->GetOutput() );
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
  
	
	// Return filtered output:
	PixelType * pixelOutData = static_cast< PixelType * > (out_rev );
	labeler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	labeler->GetOutput()->Allocate();

	/*
	 *
	 *****************************************************************/

	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying multiple Otsu's thresholding...";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		//Invoke pipeline manually in order:
	    importFilter->Update();
  
		scalarImageToHistogramGenerator->Compute();
  
		calculator->Update();

		labeler->SetInput( importFilter->GetOutput() );
	    labeler->SetRealThresholds( calculator->GetOutput()) ;
    
		labeler->Update();

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

	  
	// Print out the computed set of thresholds:

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	//const CalculatorType::OutputType &thresholdVector = labeler->GetRealThresholds(); 
	CalculatorType::OutputType::const_iterator itNum;

	if (wr_log != NULL)
	{
		for(itNum = thresholdVector.begin(); itNum < thresholdVector.end(); itNum++) 
		{			
			str  = "\tPore3D - Threshold  ";
			str += _p3dIntToString( (unsigned int) (itNum - thresholdVector.begin()) ) ;
			str += " proposed by Otsu's method is: ";
			str += _p3dIntToString ( static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>(*itNum) );
			str += " .";	
		
			wr_log ( str.c_str() );
		} 
    }


	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Multiple Otsu's thresholding applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}









int p3dMultipleOtsuThresholding3D_16 (   
	unsigned short* in_rev, 
	unsigned char* out_rev, 
	const unsigned int dimx,
	const unsigned int dimy, 
	const unsigned int dimz,
	const unsigned int numberOfRegions,
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
	
	PixelType * pixelInData = static_cast< PixelType * > (in_rev );
	importFilter->SetImportPointer( 
		pixelInData, 
		numberOfPixels, 
		false // filter won't delete buffer after use
	);
	

	/*******************************************************************
	 * Modify following section for code extensions:	 
     */

	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > 
		ScalarImageToHistogramGeneratorType;

	typedef ScalarImageToHistogramGeneratorType::HistogramType    HistogramType;

	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType >   CalculatorType;

	typedef itk::ThresholdLabelerImageFilter< ImageType, OutImageType > ThresholdLabelerType;

	
	
	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = 
		ScalarImageToHistogramGeneratorType::New();

	CalculatorType::Pointer calculator = CalculatorType::New();

	ThresholdLabelerType::Pointer labeler = ThresholdLabelerType::New();
	
	
	// Parameters settings:
	scalarImageToHistogramGenerator->SetNumberOfBins( (UCHAR_MAX + 1) / 2 );
	calculator->SetNumberOfThresholds( numberOfRegions - 1 );	

	// Setting pipeline:	
	scalarImageToHistogramGenerator->SetInput( importFilter->GetOutput() );
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
  
	
	// Return filtered output:
	OutPixelType * pixelOutData = static_cast< OutPixelType * > ( out_rev );
	labeler->GetOutput()->GetPixelContainer()->SetImportPointer(
		pixelOutData,
		numberOfPixels,
		false // filter won't delete buffer after use
		);
	labeler->GetOutput()->Allocate();

	/*
	 *
	 *****************************************************************/

	// Execute filter catching exceptions:
	try
	{	
		// Log a message:
		if (wr_log != NULL)
		{
			str  = "Pore3D - Applying multiple Otsu's thresholding...";
		
			wr_log ( str.c_str() );
		}
		
		// Start tracking computational time:
		timer.Start();

		//Invoke pipeline manually in order:
	    importFilter->Update();
  
		scalarImageToHistogramGenerator->Compute();
  
		calculator->Update();

		labeler->SetInput( importFilter->GetOutput() );
	    labeler->SetRealThresholds( calculator->GetOutput()) ;
    
		labeler->Update();

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

	  
	// Print out the computed set of thresholds:

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	//const CalculatorType::OutputType &thresholdVector = labeler->GetRealThresholds(); 
	CalculatorType::OutputType::const_iterator itNum;

	if (wr_log != NULL)
	{
		for(itNum = thresholdVector.begin(); itNum < thresholdVector.end(); itNum++) 
		{			
			str  = "\tPore3D - Threshold  ";
			str += _p3dIntToString( (unsigned int) (itNum - thresholdVector.begin()) ) ;
			str += " proposed by Otsu's method is: ";
			str += _p3dIntToString ( static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>(*itNum) );
			str += " .";	
		
			wr_log ( str.c_str() );
		}   
    }


	// Print out the elapsed time:
	if (wr_log != NULL)
	{
		str  = "Pore3D - Multiple Otsu's thresholding applied successfully in ";
		str += _p3dTimeToString ( timer.GetMean() );		
		str += " seconds.";
		
		wr_log ( str.c_str() );
	}

	// Return OK:
	return P3D_SUCCESS;
}







