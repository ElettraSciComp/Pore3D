#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "p3dBlob.h"
#include "p3dFilt.h"
#include "p3dSkel.h"

int customPrint_nolf(const char *msg, ...) {
	va_list fmtargs;
	char buffer[1024];
	
	
	//FILE *foutVF;

	va_start(fmtargs, msg);
	vsnprintf(buffer, sizeof (buffer) - 1, msg, fmtargs);
	va_end(fmtargs);

	return printf("%s", buffer);

	/*
	//printf( "%s\n", buffer );
	// open the file
	if ((foutVF = fopen("C:\\loglog.txt","a")) == NULL)
	{
	printf("Cannot open output file loglog.txt for appending");
	return -1;
	}

	fprintf(foutVF, "%s\n", buffer );

	// close the file
	fclose(foutVF);
	return 0;*/
}

int customPrint(const char *msg, ...) {
	va_list fmtargs;
	char buffer[1024];
	//FILE *foutVF;

	va_start(fmtargs, msg);
	vsnprintf(buffer, sizeof (buffer) - 1, msg, fmtargs);
	va_end(fmtargs);

	return printf("%s\n", buffer);

	/*
	//printf( "%s\n", buffer );
	// open the file
	if ((foutVF = fopen("C:\\loglog.txt","a")) == NULL)
	{
	printf("Cannot open output file loglog.txt for appending");
	return -1;
	}

	fprintf(foutVF, "%s\n", buffer );

	// close the file
	fclose(foutVF);
	return 0;*/
}

int customProgress(const int msg, ...) {
	return customPrint_nolf("\tProgress: %d.\r", msg);

}

int main(int argc, char* argv[]) {
	const int dimx = atoi(argv[3]);
	const int dimy = atoi(argv[4]);
	const int dimz = atoi(argv[5]);	

	/*unsigned char* nodes_im = NULL;
	unsigned char* ends_im = NULL;
	unsigned char* throats_im = NULL;*/

	//struct SkeletonStats skl_stats;
	//BlobStats* blob_stats;

	// Allocate memory for input and output images:
	unsigned char* in_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char));
	unsigned char* out_im = (unsigned char*) malloc(dimx * dimy * dimz * sizeof (unsigned char));

	//skl_stats = (SkeletonStats*) malloc(sizeof(SkeletonStats));
	//blob_stats = (BlobStats*) malloc(sizeof(BlobStats));

	//omp_set_num_threads(8);

	// Read the slice:
	p3dReadRaw8 ( argv[1], in_im, dimx, dimy, dimz, customPrint, NULL);

	//p3dReadRaw8(argv[1], in_im, dimx, dimy, dimz, customPrint, NULL);


	//p3dSingleRegionGrowing3D_8(in_im, skl_im, dimx, dimy, dimz, 75, 75, 10, 0.2, CONN6, customPrint, NULL);


	//p3dMunchEtAlRingRemover2D_8(in_im, out_im, dimx, dimy, 991.5, 991.5, 8, 2.5, 1, 1.5, customPrint , NULL);
	//p3dAnisotropicDiffusionFilter3D_8(in_im, out_im, dimx, dimy, dimz, 1, 0.01, 0.01, 50, customPrint, customProgress);
	p3dBilateralFilter3D_8(in_im, out_im, dimx, dimy, dimz, 3, 1.0, 10.0, 50, customPrint, customProgress );
	//p3dGaussianFilter3D_8(in_im, out_im, dimx, dimy, dimz, 3, 30.0, customPrint, customProgress );

	//p3dBlobLabeling_ushort(in_im, out_im, dimx, dimy, dimz, CONN6, P3D_TRUE, P3D_FALSE, customPrint);
	//p3dBlobAnalysis(in_im, blob_stats, out_im, out_im, dimx, dimy, dimz, 1.0, CONN26, 1024, P3D_TRUE, customPrint);


	//p3dSkeletonAnalysis2(in_im, skl_im, &skl_stats, nodes_im, ends_im, throats_im, dimx, dimy, dimz, 0.00298, customPrint);
	//p3dSkeletonAnalysis2(in_im, skl_im, &skl_stats, NULL, NULL, NULL, dimx, dimy, dimz, 0.00298, customPrint);


	//p3dWriteRaw16(out_im, argv[2], dimx, dimy, dimz, P3D_TRUE, P3D_FALSE, customPrint, NULL);
	p3dWriteRaw8(out_im, argv[2], dimx, dimy, dimz);


	// Release resources:	
	if (in_im != NULL) free(in_im);
	if (out_im != NULL) free(out_im);
	/*if (nodes_im != NULL) free(nodes_im);
	if (ends_im != NULL) free(ends_im);
	if (throats_im != NULL) free(throats_im);*/

	// Free C memory:
	/* if (skl_stats.Node_Width != NULL) free ( skl_stats.Node_Width );
	if (skl_stats.End_Width != NULL) free ( skl_stats.End_Width );

	if (skl_stats.EndToEnd_Length != NULL) free ( skl_stats.EndToEnd_Length );
	if (skl_stats.EndToEnd_MinWidth != NULL) free ( skl_stats.EndToEnd_MinWidth );
	if (skl_stats.EndToEnd_MeanWidth != NULL) free ( skl_stats.EndToEnd_MeanWidth );
	if (skl_stats.EndToEnd_MaxWidth != NULL) free ( skl_stats.EndToEnd_MaxWidth );

	if (skl_stats.NodeToEnd_Length != NULL) free ( skl_stats.NodeToEnd_Length );
	if (skl_stats.NodeToEnd_MinWidth != NULL) free ( skl_stats.NodeToEnd_MinWidth );
	if (skl_stats.NodeToEnd_MeanWidth != NULL) free ( skl_stats.NodeToEnd_MeanWidth );
	if (skl_stats.NodeToEnd_MaxWidth != NULL) free ( skl_stats.NodeToEnd_MaxWidth );

	if (skl_stats.NodeToNode_Length != NULL) free ( skl_stats.NodeToNode_Length );
	if (skl_stats.NodeToNode_MinWidth != NULL) free ( skl_stats.NodeToNode_MinWidth );
	if (skl_stats.NodeToNode_MeanWidth != NULL) free ( skl_stats.NodeToNode_MeanWidth );
	if (skl_stats.NodeToNode_MaxWidth != NULL) free ( skl_stats.NodeToNode_MaxWidth );

	if (skl_stats.CoordinationNumber != NULL ) free ( skl_stats.CoordinationNumber ); */

	//if (skl_stats != NULL ) free ( skl_stats );


	//return 0;
}