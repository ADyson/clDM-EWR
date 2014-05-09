// MyDlg2.cpp : implementation file
//
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include "stdafx.h"
#include "MyDlg2.h"
#include "NumEdit.h"
#include "fftw3.h"
#include "boost/lexical_cast.hpp"
#include <omp.h>

#define PLATFORM 0

/*
#define STRINGIFY(src) #src

inline const char* Kernels() {
  static const char* kernels =
    #include "pcpcf.cl"
    ;
  return kernels;
}
*/

const char* pcpcfSource = 
"__kernel void clPCPCF(__global const double2* fft1, __global const double2* fft2, __global double2* fftresult, __global const double* CLxFrequencies, __global const double* CLyFrequencies, int sizeX, int sizeY, float focalstep, float wavelength) \n"
"{		\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]); \n"
"		double compensation = cos(3.1415926534 * focalstep * wavelength * frequency * frequency); \n"
"		double c1r = fft1[Index].x; \n"
"		double c1i = fft1[Index].y; \n"
"		double c2r = fft2[Index].x; \n"
"		double c2i = fft2[Index].y; \n"
"		double denom2 = sqrt( compensation*compensation   * (c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i) \n" 
"		+ 0.000000001 * 0.000000001 + 2 * 0.000000001 * compensation * (c1r*c2r + c1i*c2i)); \n"
"		fftresult[Index].x = (frequency<=4)*compensation*(c1r*c2r + c1i*c2i)/denom2; \n"
"		fftresult[Index].y = (frequency<=4)*compensation*(c2i*c1r - c1i*c2r)/denom2; \n"
"	}	\n"
"}		\n"
;

void PCPCF(float* seriesdata,int im1,int im2 ,fftw_plan fftplan,fftw_plan fftplanInv,fftw_complex* complexImageOne,fftw_complex* complexImageTwo,int &xShiftAmount,int &yShiftAmount, int sizeX, int sizeY, float* xFrequencies, float* yFrequencies, float wavelength, float fStep, int iLeft,int iTop, int xDim, int yDim)
{
	fftw_complex* dataOne = new fftw_complex[sizeX*sizeY];
	fftw_complex* dataTwo = new fftw_complex[sizeX*sizeY];

	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			complexImageOne[i+j*sizeX][0] = seriesdata[im1*xDim*yDim + i + iLeft + (j+iTop)*xDim];
			complexImageOne[i+j*sizeX][1] = 0;
		}

	fftw_execute(fftplan);

	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			dataOne[i+j*sizeX][0] = complexImageTwo[i+j*sizeX][0]; 
			dataOne[i+j*sizeX][1] = complexImageTwo[i+j*sizeX][1];
		}

	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			complexImageOne[i+j*sizeX][0] = seriesdata[im2*xDim*yDim + i + iLeft + (j+iTop)*xDim];
			complexImageOne[i+j*sizeX][1] = 0;
		}

	fftw_execute(fftplan);

	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			dataTwo[i+j*sizeX][0] = complexImageTwo[i+j*sizeX][0]; 
			dataTwo[i+j*sizeX][1] = complexImageTwo[i+j*sizeX][1];
		}


		// PCPCF time...

	for(uint i = 0; i < sizeX;i++)
		for(uint j = 0 ; j < sizeY; j++)
		{
			double frequency = sqrt(xFrequencies[i]*xFrequencies[i] + yFrequencies[j]*yFrequencies[j]);
			double compensation = cos(3.1415926534 * fStep * wavelength * frequency * frequency);
			
			double c1r = dataOne[i + j*sizeX][0];
			double c1i = dataOne[i + j*sizeX][1];
			double c2r = dataTwo[i + j*sizeX][0];
			double c2i = dataTwo[i + j*sizeX][1];
			
			double denom2 = sqrt( compensation*compensation   * (c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i) 
			+ 0.000000001 * 0.000000001 + 2 * 0.000000001 * compensation * (c1r*c2r + c1i*c2i));
			
			complexImageTwo[i + j*sizeX][0] = (frequency<=4)*compensation*(c1r*c2r + c1i*c2i)/denom2;
			complexImageTwo[i + j*sizeX][1] = (frequency<=4)*compensation*(c2i*c1r - c1i*c2r)/denom2;	
		}

	fftw_execute(fftplanInv);

	// Display PCPCF
	//DM::Image pcpcf = DM::RealImage("PCPCF",4,sizeX,sizeY);
	//Gatan::PlugIn::ImageDataLocker pcpcfLocker(pcpcf);
	//float* pcpcfdata = (float*) pcpcfLocker.get();

	float maxHeight1 = 0.0f;
	int maxPosition1 = 0;
	for(uint j = 0; j< sizeX*sizeY;j++)
	{	
		float val =  sqrt(complexImageOne[j][0]*complexImageOne[j][0] + complexImageOne[j][1]*complexImageOne[j][1]);
		//pcpcfdata[j] = val;

		if(val > maxHeight1)
		{
			maxHeight1 = val;
			maxPosition1 = j;
		}
	}

	// Find shift from max height position
	int xShift;
	int yShift;

	// Translate linear array index into row and column.
	int maxindexr = floor(float((maxPosition1)/(sizeX)))+1;
	int maxindexc = maxPosition1 + 1 - sizeX*floor(float((maxPosition1)/(sizeX)));


	// Shift is positive or negative depending on image quadrant it appears in.
	if(maxindexr > floor(sizeY/2 +.5))
		yShift = maxindexr - (sizeY) -1;
	else
		yShift = maxindexr - 1;

	if(maxindexc > floor(sizeX/2 +.5))
		xShift = maxindexc - (sizeX) -1;
	else
		xShift = maxindexc - 1;

	//pcpcfLocker.~ImageDataLocker();
	//DM::ImageDocument  pcpcfDoc = pcpcf.GetOrCreateImageDocument();
	//pcpcfDoc.Show();
	//pcpcf.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcpcf.GetImageDisplay(0).SetSurveyTechnique(4);

	xShiftAmount = xShift;
	yShiftAmount = yShift;

	delete[] dataOne;
	delete[] dataTwo;

}

void OpenCLPCPCF(float* seriesdata,int im1,int im2 ,cl_command_queue cmdQueue, clAmdFftPlanHandle fftplan,cl_mem fftInputBuffer,cl_mem fftOutputBuffer,cl_mem dataBuffer1,cl_mem dataBuffer2,cl_mem clMedBuffer,cl_kernel clPCPCF, int &xShiftAmount,int &yShiftAmount, int sizeX, int sizeY, cl_mem CLxFrequencies, cl_mem CLyFrequencies, float wavelength, float fStep, int iLeft,int iTop, int xDim, int yDim)
{
	std::vector< std::complex< double > > dataOne( sizeX*sizeY );
	std::vector< std::complex< double > > dataTwo( sizeX*sizeY );

	cl_uint status;
	clAmdFftStatus fftStatus;
	
	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			dataOne[i+j*sizeX] = seriesdata[im1*xDim*yDim + i + iLeft + (j+iTop)*xDim];
		}

	

	clEnqueueWriteBuffer( cmdQueue, fftInputBuffer, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &dataOne[ 0 ], 
				0, NULL, NULL );

	fftStatus = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&fftInputBuffer, &dataBuffer1, clMedBuffer );


	

	for(uint j = 0; j< sizeY;j++)
		for(uint i = 0; i< sizeX;i++)
		{
			dataTwo[i+j*sizeX] = seriesdata[im2*xDim*yDim + i + iLeft + (j+iTop)*xDim];
		}

	clEnqueueWriteBuffer( cmdQueue, fftInputBuffer, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &dataTwo[ 0 ], 
				0, NULL, NULL );

	fftStatus = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&fftInputBuffer, &dataBuffer2, clMedBuffer );

	
	status = clSetKernelArg(clPCPCF,0,sizeof(cl_mem),&dataBuffer1);
	status |= clSetKernelArg(clPCPCF,1,sizeof(cl_mem),&dataBuffer2);
	status |= clSetKernelArg(clPCPCF,2,sizeof(cl_mem),&fftInputBuffer);
	status |= clSetKernelArg(clPCPCF,3,sizeof(cl_mem),&CLxFrequencies);
	status |= clSetKernelArg(clPCPCF,4,sizeof(cl_mem),&CLyFrequencies);
	status |= clSetKernelArg(clPCPCF,5,sizeof(int),&sizeX);
	status |= clSetKernelArg(clPCPCF,6,sizeof(int),&sizeY);
	status |= clSetKernelArg(clPCPCF,7,sizeof(float),&fStep);
	status |= clSetKernelArg(clPCPCF,8,sizeof(float),&wavelength);

	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[2];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;


	// Messed up when I tried this
	size_t* localWorkSize = new size_t[2];
	localWorkSize[0] = 32;
	localWorkSize[1] = 16;

	status = clEnqueueNDRangeKernel(cmdQueue,clPCPCF,2,NULL,globalWorkSize,NULL,0,NULL,NULL);

	// Now Inverse FFT
	fftStatus = clAmdFftEnqueueTransform( fftplan, CLFFT_BACKWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&fftInputBuffer, &fftOutputBuffer, clMedBuffer );

	//Now read back
	clFinish(cmdQueue);

	if(fftStatus != CL_SUCCESS)
	{
		DM::Result("FFT Failed");
	}

	if(status != CL_SUCCESS)
	{
		DM::Result("OpenCL Error");
	}

	clEnqueueReadBuffer( cmdQueue, fftOutputBuffer, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &dataOne[ 0 ], 
				0, NULL, NULL );



	// Display PCPCF
	//DM::Image pcpcf = DM::RealImage("PCPCF",8,sizeX,sizeY);
	//Gatan::PlugIn::ImageDataLocker pcpcfLocker(pcpcf);
	//double* pcpcfdata = (double*) pcpcfLocker.get();

	float maxHeight1 = 0.0f;
	int maxPosition1 = 0;
	for(uint j = 0; j< sizeX*sizeY;j++)
	{	
		double val =  sqrt(dataOne[j].real()*dataOne[j].real() + dataOne[j].imag()*dataOne[j].imag());
	//	pcpcfdata[j] = val;

		if(val > maxHeight1)
		{
			maxHeight1 = val;
			maxPosition1 = j;
		}
	}

	// Find shift from max height position
	int xShift;
	int yShift;

	// Translate linear array index into row and column.
	int maxindexr = floor(float((maxPosition1)/(sizeX)))+1;
	int maxindexc = maxPosition1 + 1 - sizeX*floor(float((maxPosition1)/(sizeX)));


	// Shift is positive or negative depending on image quadrant it appears in.
	if(maxindexr > floor(sizeY/2 +.5))
		yShift = maxindexr - (sizeY) -1;
	else
		yShift = maxindexr - 1;

	if(maxindexc > floor(sizeX/2 +.5))
		xShift = maxindexc - (sizeX) -1;
	else
		xShift = maxindexc - 1;

	//pcpcfLocker.~ImageDataLocker();
	//DM::ImageDocument  pcpcfDoc = pcpcf.GetOrCreateImageDocument();
	//pcpcfDoc.Show();
	//pcpcf.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcpcf.GetImageDisplay(0).SetSurveyTechnique(4);

	xShiftAmount = xShift;
	yShiftAmount = yShift;

	dataOne.clear();
	dataTwo.clear();
}

float GetSum(float* data, int length, int start)
{
	float sum=0;

	for(int i = start; i < start + length; i++)
		sum+=data[i];

	return sum;
}


void wavetransferfunction (fftw_complex* input, float * k0x, float * k0y, int nkx, int nky, float wavel, float beta,
									  float delta, float A1r, float A1i, float A2, float B2, float C1, float C3,float objap)
{


	for(int xIndex = 0; xIndex < nkx; xIndex++)
		for(int yIndex = 0; yIndex < nky; yIndex++)
		{
			int Index = (yIndex * nkx + xIndex);
			
			float k = sqrt( k0x[xIndex] * k0x[xIndex] + k0y[yIndex] * k0y[yIndex] );
			float kx = k0x[xIndex];
			float ky = k0y[yIndex];
			
			float SpatCohPart = fabs( wavel * 2 * 3.141592654 * C1 * k + 2 * 3.141592654 * wavel * wavel * wavel * C3 * k * k * k );
			
			float CohEnv = exp( -(( beta * beta )/( 4 * wavel * wavel )) * SpatCohPart * SpatCohPart ) *
						   exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavel * k * k ) * fabs( 3.141592654 * wavel * k * k )));
						
			
			float gamma = ( .166666666 * 3.1415927 * wavel ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky 
				  + 6 * C1 * ky * ky + 3 * C3 * wavel * wavel * kx * kx * kx * kx + 6 * C3 * wavel * wavel * kx * kx * ky * ky
				  + 3 * C3 * wavel * wavel * ky * ky * ky * ky + 12 * A1i * kx * ky) ;
						
			input[Index][0] = (k<=objap) * CohEnv * cosf(gamma); 
			input[Index][1] = (k<=objap) * CohEnv * -1 * sinf(gamma);

			if(Index==0)
			{
				input[Index][0] = 0.5; 
				input[Index][1] = 0;
			}
		}
}

const char* wavetransferfunctionsource = 
"__kernel void clWaveTransferFunction(__global double2* clw, __global const double* CLxFrequencies, __global const double* CLyFrequencies, int sizeX, int sizeY, float wavelength, float beta, float delta, float A1r, float A1i, float C1, float C3, float objap)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]);	\n"
"		double kx = CLxFrequencies[xid];	\n"
"		double ky = CLyFrequencies[yid];	\n"
"		double SpatCohPart = fabs( wavelength * 2 * 3.141592654 * C1 * frequency + 2 * 3.141592654 * wavelength * wavelength * wavelength * C3 * frequency * frequency * frequency );	\n"
"		double CohEnv = exp( -(( beta * beta )/( 4 * wavelength * wavelength )) * SpatCohPart * SpatCohPart ) *	\n"
"						exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavelength * frequency * frequency ) * fabs( 3.141592654 * wavelength * frequency * frequency )));	\n"
"		double gamma = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky	\n"
"				  + 6 * C1 * ky * ky + 3 * C3 * wavelength * wavelength * kx * kx * kx * kx + 6 * C3 * wavelength * wavelength * kx * kx * ky * ky	\n"
"				  + 3 * C3 * wavelength * wavelength * ky * ky * ky * ky + 12 * A1i * kx * ky) ;	\n"
"		clw[Index].x = (frequency<=objap) * CohEnv * cos(gamma);	\n"
"		clw[Index].y = (frequency<=objap) * CohEnv * -1 * sin(gamma);	\n"
"		if(Index==0)	\n"
"		{\n"
"			clw[Index].x = 0.5; \n" 
"			clw[Index].y = 0; \n"
"		}\n"
"	}	\n"
"}	\n"
;

const char* wavetransferfunctionminussource = 
"__kernel void clWaveTransferFunctionMinus(__global double2* clw, __global const double* CLxFrequencies, __global const double* CLyFrequencies, int sizeX, int sizeY, float wavelength, float beta, float delta, float A1r, float A1i, float C1, float C3, float objap)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]);	\n"
"		double kx = -CLxFrequencies[xid];	\n"
"		double ky = -CLyFrequencies[yid];	\n"
"		double SpatCohPart = fabs( wavelength * 2 * 3.141592654 * C1 * frequency + 2 * 3.141592654 * wavelength * wavelength * wavelength * C3 * frequency * frequency * frequency );	\n"
"		double CohEnv = exp( -(( beta * beta )/( 4 * wavelength * wavelength )) * SpatCohPart * SpatCohPart ) *	\n"
"						exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavelength * frequency * frequency ) * fabs( 3.141592654 * wavelength * frequency * frequency )));	\n"
"		float gamma = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky	\n"
"				  + 6 * C1 * ky * ky + 3 * C3 * wavelength * wavelength * kx * kx * kx * kx + 6 * C3 * wavelength * wavelength * kx * kx * ky * ky	\n"
"				  + 3 * C3 * wavelength * wavelength * ky * ky * ky * ky + 12 * A1i * kx * ky) ;	\n"
"		clw[Index].x = (frequency<=objap) * CohEnv * cos(gamma);	\n"
"		clw[Index].y = (frequency<=objap) * CohEnv * -1 * sin(gamma);	\n"
"		if(Index==0)	\n"
"		{\n"
"			clw[Index].x = 0.5; \n" 
"			clw[Index].y = 0; \n"
"		}\n"
"	}	\n"
"}	\n"
;


void wavetransferfunctionminus (fftw_complex* input, float * k0x, float * k0y, int nkx, int nky, float wavel, float beta,
									  float delta, float A1r, float A1i, float A2, float B2, float C1, float C3,float objap)
{


	for(int xIndex = 0; xIndex < nkx; xIndex++)
		for(int yIndex = 0; yIndex < nky; yIndex++)
		{
			int Index = (yIndex * nkx + xIndex);
			
			float k = sqrt( k0x[xIndex] * k0x[xIndex] + k0y[yIndex] * k0y[yIndex] );
			float kx = -k0x[xIndex];
			float ky = -k0y[yIndex];
			
			float SpatCohPart = fabs( wavel * 2 * 3.141592654 * C1 * k + 2 * 3.141592654 * wavel * wavel * wavel * C3 * k * k * k );
			
			float CohEnv = exp( -(( beta * beta )/( 4 * wavel * wavel )) * SpatCohPart * SpatCohPart ) *
						   exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavel * k * k ) * fabs( 3.141592654 * wavel * k * k )));
						
			float gamma = ( .166666666 * 3.141592654 * wavel ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky 
				  + 6 * C1 * ky * ky + 3 * C3 * wavel * wavel * kx * kx * kx * kx + 6 * C3 * wavel * wavel * kx * kx * ky * ky
				  + 3 * C3 * wavel * wavel * ky * ky * ky * ky + 12 * A1i * kx * ky) ;
						
			input[Index][0] = (k<=objap) * CohEnv * cosf(gamma); 
			input[Index][1] = (k<=objap) * CohEnv * -1 * sinf(gamma);

			if(Index==0)
			{
				input[Index][0] = 0.5; 
				input[Index][1] = 0;
			}
		}
}

void wavetransferfunctionmtf (fftw_complex* input, float * k0x, float * k0y, int nkx, int nky, float wavel, float beta,
									  float delta, float A1r, float A1i, float A2, float B2, float C1, float C3,float objap, float* mtfdata, int mtflength, float scaleMTFx, float scaleMTFy)
{


	for(int xIndex = 0; xIndex < nkx; xIndex++)
		for(int yIndex = 0; yIndex < nky; yIndex++)
		{
			int Index = (yIndex * nkx + xIndex);
			
			float k = sqrt( k0x[xIndex] * k0x[xIndex] + k0y[yIndex] * k0y[yIndex] );
			float kx = k0x[xIndex];
			float ky = k0y[yIndex];

			int xposmtf;
			int yposmtf;
			
			if(xIndex < nkx/2)
				xposmtf = xIndex;
			else
				xposmtf = (nkx-xIndex);
				
			if(yIndex < nky/2)
				yposmtf = yIndex;
			else
				yposmtf = nky-yIndex;

			float mtfposition = sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);
				
			int mtfpos = floor(sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));
			
			float interp = mtfposition - float(mtfpos);
			
			if(mtfpos>mtflength-2)
				mtfpos=mtflength-2;
			
			float mtf = (1-interp)*mtfdata[mtfpos] + interp*mtfdata[mtfpos+1];
			
			float SpatCohPart = fabs( wavel * 2 * 3.141592654 * C1 * k + 2 * 3.141592654 * wavel * wavel * wavel * C3 * k * k * k );
			
			float CohEnv = exp( -(( beta * beta )/( 4 * wavel * wavel )) * SpatCohPart * SpatCohPart ) *
						   exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavel * k * k ) * fabs( 3.141592654 * wavel * k * k )));
						
			float gamma = ( .166666666 * 3.141592654 * wavel ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky 
				  + 6 * C1 * ky * ky + 3 * C3 * wavel * wavel * kx * kx * kx * kx + 6 * C3 * wavel * wavel * kx * kx * ky * ky
				  + 3 * C3 * wavel * wavel * ky * ky * ky * ky + 12 * A1i * kx * ky) ;
						
			input[Index][0] = (k<=objap) * mtf * CohEnv * cosf(gamma); 
			input[Index][1] = (k<=objap) * mtf * CohEnv * -1 * sinf(gamma);

			if(Index==0)
			{
				input[Index][0] = 0.5; 
				input[Index][1] = 0;
			}
		}
}

const char* wavetransferfunctionmtfsource = 
"__kernel void clWaveTransferFunctionMTF(__global double2* clw, __global const double* CLxFrequencies, __global const double* CLyFrequencies, int sizeX, int sizeY, \n"
" float wavelength, float beta, float delta, float A1r, float A1i, float C1, float C3, float objap, __global const float* clMTF, int mtflength, float scaleMTFx, float scaleMTFy)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; 																			\n"
"		double frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]);	\n"
"		double kx = CLxFrequencies[xid];																		\n"
"		double ky = CLyFrequencies[yid];																		\n"
"																												\n"
"		int xposmtf;																							\n"
"		int yposmtf;																							\n"
"																												\n"
"		if(xid < sizeX/2)																						\n"
"			xposmtf = xid;																						\n"
"		else																									\n"
"			xposmtf = (sizeX-xid);																				\n"
"																												\n"
"		if(yid < sizeY/2)																						\n"
"			yposmtf = yid;																						\n"
"		else																									\n"
"			yposmtf = sizeY-yid;																				\n"
"																												\n"
"		float mtfposition = sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);	\n"
"																												\n"
"		int mtfpos = floor(sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));	\n"
"																												\n"		
"		float interp = mtfposition - floor(mtfposition);														\n"
"																												\n"		
"		if(mtfpos>mtflength-2)																					\n"
"			mtfpos=mtflength-2;																					\n"
"																												\n"
"		float mtf = (1-interp)*clMTF[mtfpos] + interp*clMTF[mtfpos+1];											\n"
"		double SpatCohPart = fabs( wavelength * 2 * 3.141592654 * C1 * frequency + 2 * 3.141592654 * wavelength * wavelength * wavelength * C3 * frequency * frequency * frequency );	\n"
"		double CohEnv = exp( -(( beta * beta )/( 4 * wavelength * wavelength )) * SpatCohPart * SpatCohPart ) *	\n"
"						exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavelength * frequency * frequency ) * fabs( 3.141592654 * wavelength * frequency * frequency )));	\n"
"		float gamma = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky	\n"
"				  + 6 * C1 * ky * ky + 3 * C3 * wavelength * wavelength * kx * kx * kx * kx + 6 * C3 * wavelength * wavelength * kx * kx * ky * ky	\n"
"				  + 3 * C3 * wavelength * wavelength * ky * ky * ky * ky + 12 * A1i * kx * ky) ;	\n"
"		clw[Index].x = (frequency<=objap) * mtf * CohEnv * cos(gamma);	\n"
"		clw[Index].y = (frequency<=objap) * mtf * CohEnv * -1 * sin(gamma);	\n"
"		if(Index==0)	\n"
"		{\n"
"			clw[Index].x = 0.5; \n" 
"			clw[Index].y = 0;	\n"
"		}\n"
"	}	\n"
"}	\n"
;

void wavetransferfunctionminusmtf (fftw_complex* input, float * k0x, float * k0y, int nkx, int nky, float wavel, float beta,
									  float delta, float A1r, float A1i, float A2, float B2, float C1, float C3,float objap, float* mtfdata, int mtflength, float scaleMTFx, float scaleMTFy)
{


	for(int xIndex = 0; xIndex < nkx; xIndex++)
		for(int yIndex = 0; yIndex < nky; yIndex++)
		{
			int Index = (yIndex * nkx + xIndex);
			
			float k = sqrt( k0x[xIndex] * k0x[xIndex] + k0y[yIndex] * k0y[yIndex] );
			float kx = -k0x[xIndex];
			float ky = -k0y[yIndex];

			int xposmtf;
			int yposmtf;
			
			if(xIndex < nkx/2)
				xposmtf = xIndex;
			else
				xposmtf = (nkx-xIndex);
				
			if(yIndex < nky/2)
				yposmtf = yIndex;
			else
				yposmtf = nky-yIndex;
				
			float mtfposition = sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);
				
			int mtfpos = floor(sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));
			
			float interp = mtfposition - float(mtfpos);
			
			if(mtfpos>mtflength-2)
				mtfpos=mtflength-2;
			
			float mtf = (1-interp)*mtfdata[mtfpos] + interp*mtfdata[mtfpos+1];
			
			float SpatCohPart = fabs( wavel * 2 * 3.141592654 * C1 * k + 2 * 3.141592654 * wavel * wavel * wavel * C3 * k * k * k );
			
			float CohEnv = exp( -(( beta * beta )/( 4 * wavel * wavel )) * SpatCohPart * SpatCohPart ) *
						   exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavel * k * k ) * fabs( 3.141592654 * wavel * k * k )));
						
			float gamma = ( .166666666 * 3.1415927 * wavel ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky 
				  + 6 * C1 * ky * ky + 3 * C3 * wavel * wavel * kx * kx * kx * kx + 6 * C3 * wavel * wavel * kx * kx * ky * ky
				  + 3 * C3 * wavel * wavel * ky * ky * ky * ky + 12 * A1i * kx * ky) ;
						
			input[Index][0] = (k<=objap) * mtf * CohEnv * cosf(gamma); 
			input[Index][1] = (k<=objap) * mtf * CohEnv * -1 * sinf(gamma);

			if(Index==0)
			{
				input[Index][0] = 0.5; 
				input[Index][1] = 0;
			}
		}
}

const char* wavetransferfunctionminusmtfsource = 
"__kernel void clWaveTransferFunctionMinusMTF(__global double2* clw, __global double* CLxFrequencies, __global double* CLyFrequencies, int sizeX, int sizeY, float wavelength, float beta, float delta, float A1r, float A1i, float C1, float C3, float objap, __global float* clMTF, int mtflength, float scaleMTFx, float scaleMTFy)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]);	\n"
"		double kx = -CLxFrequencies[xid];	\n"
"		double ky = -CLyFrequencies[yid];	\n"
"											\n"
"		int xposmtf;						\n"
"		int yposmtf;						\n"
"											\n"
"		if(xid < sizeX/2)					\n"
"			xposmtf = xid;					\n"
"		else								\n"
"			xposmtf = (sizeX-xid);			\n"
"											\n"
"		if(yid < sizeY/2)					\n"
"			yposmtf = yid;					\n"
"		else								\n"
"			yposmtf = sizeY-yid;			\n"
"											\n"
"		float mtfposition = sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);	\n"
"																												\n"
"		int mtfpos = floor(sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));	\n"
"																												\n"		
"		float interp = mtfposition - floor(mtfposition);														\n"
"																												\n"		
"		if(mtfpos>mtflength-2)																					\n"
"			mtfpos=mtflength-2;																					\n"
"																												\n"
"		float mtf = (1-interp)*clMTF[mtfpos] + interp*clMTF[mtfpos+1];											\n"
"		double SpatCohPart = fabs( wavelength * 2 * 3.141592654 * C1 * frequency + 2 * 3.141592654 * wavelength * wavelength * wavelength * C3 * frequency * frequency * frequency );	\n"
"		double CohEnv = exp( -(( beta * beta )/( 4 * wavelength * wavelength )) * SpatCohPart * SpatCohPart ) *	\n"
"						exp( -(0.5 * delta * delta)*(fabs( 3.141592654 * wavelength * frequency * frequency ) * fabs( 3.141592654 * wavelength * frequency * frequency )));	\n"
"		float gamma = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky	\n"
"				  + 6 * C1 * ky * ky + 3 * C3 * wavelength * wavelength * kx * kx * kx * kx + 6 * C3 * wavelength * wavelength * kx * kx * ky * ky	\n"
"				  + 3 * C3 * wavelength * wavelength * ky * ky * ky * ky + 12 * A1i * kx * ky) ;	\n"
"		clw[Index].x = (frequency<=objap) * mtf * CohEnv * cos(gamma);	\n"
"		clw[Index].y = (frequency<=objap) * mtf * CohEnv * -1 * sin(gamma);	\n"
"		if(Index==0)	\n"
"		{\n"
"			clw[Index].x = 0.5; \n" 
"			clw[Index].y = 0; \n"
"		}\n"
"	}	\n"
"}\n"
;


const char* wienerwsource = 
"__kernel void clWienerW(__global double* clW, __global double2* clw, int sizeX, int sizeY, int init)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double absw = sqrt(clw[Index].x*clw[Index].x + clw[Index].y*clw[Index].y);	\n"
"		if(init==1)	\n"
"		{\n"
"			clW[Index] = 0; \n" 
"		}\n"
"		clW[Index] += absw*absw;	\n"
"	}	\n"
"}	\n"
;


const char* wienervsource = 
"__kernel void clWienerV(__global double2* clV, __global double2* clw, __global double2* clwminus, int sizeX, int sizeY, int init)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double real = clw[Index].x * clwminus[Index].x - clw[Index].y * clwminus[Index].y;	\n"
"		double imag = clw[Index].x * clwminus[Index].y + clw[Index].y * clwminus[Index].x;	\n"
"		if(init==1)	\n"
"		{\n"
"			clV[Index].x = 0; \n" 
"			clV[Index].y = 0; \n" 
"		}\n"
"		clV[Index].x += real;	\n"
"		clV[Index].y += imag;	\n"
"	}	\n"
"}	\n"
;


const char* wienertsource = 
"__kernel void clWienerT(__global double2* clT, __global double2* clw, __global double2* fft, int sizeX, int sizeY, int init)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double real = clw[Index].x*fft[Index].x + clw[Index].y*fft[Index].y;	\n"
"		double imag = clw[Index].x*fft[Index].y - clw[Index].y*fft[Index].x;	\n"
"		if(init==1)	\n"
"		{\n"
"			clT[Index].x = 0; \n" 
"			clT[Index].y = 0; \n" 
"		}\n"
"		clT[Index].x += real;	\n"
"		clT[Index].y += imag;	\n"
"	}	\n"
"}	\n"
;



const char* wienerusource = 
"__kernel void clWienerU(__global double2* clU, __global double2* clwminus, __global double2* fft, int sizeX, int sizeY, int init)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double real = clwminus[Index].x*fft[Index].x - clwminus[Index].y*fft[Index].y;	\n"
"		double imag = clwminus[Index].x*fft[Index].y + clwminus[Index].y*fft[Index].x;	\n"
"		if(init==1)	\n"
"		{\n"
"			clU[Index].x = 0; \n" 
"			clU[Index].y = 0; \n" 
"		}\n"
"		clU[Index].x += real;	\n"
"		clU[Index].y += imag;	\n"
"	}	\n"
"}	\n"
;


const char* makerestoredsource = 
"__kernel void clMakeRestored(__global double* clW,__global double* clWminus, __global double2* clV,__global double2* clT, __global double2* clU, __global double2* clRestored, int sizeX, int sizeY)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double real = (clWminus[Index] + 0.2)*clT[Index].x - clV[Index].x * clU[Index].x - clV[Index].y * clU[Index].y;	\n"
"		double imag = (clWminus[Index] + 0.2)*clT[Index].y - clV[Index].x * clU[Index].y + clV[Index].y * clU[Index].x;	\n"
"		double denom = ((clWminus[Index] + 0.2) * (clW[Index] + 0.2)) - ( clV[Index].x * clV[Index].x + clV[Index].y * clV[Index].y ); \n"
"		clRestored[Index].x = real/denom;	\n"
"		clRestored[Index].y = imag/denom;	\n"
"	}	\n"
"}	\n"
;

void MakeRestoredWavefunctionmtfnps (double * W, double * Wminus, fftw_complex * V, fftw_complex * T, fftw_complex * U, fftw_complex * Restored, int nkx, int nky, float* mtfdata, float* npsdata, int mtflength, int npslength, float scaleMTFx, float scaleMTFy)
{
	for(int xIndex = 0; xIndex < nkx; xIndex++)
		for(int yIndex = 0; yIndex < nky; yIndex++)
			{
				int Index = (yIndex * nkx + xIndex);
				
				int xposmtf;
				int yposmtf;
				
				if(xIndex < nkx/2)
					xposmtf = xIndex;
				else
					xposmtf = nkx-xIndex;
					
				if(yIndex < nky/2)
					yposmtf = yIndex;
				else
					yposmtf = nky-yIndex;
					
				float mtfposition = sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);
				
				int mtfpos = floor(sqrtf(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));
			
				float interp = mtfposition - float(mtfpos);
			
				if(mtfpos>mtflength-2)
					mtfpos=mtflength-2;
			
				float mtf = (1-interp)*mtfdata[mtfpos] + interp*mtfdata[mtfpos+1];
				float nps = (1-interp)*npsdata[mtfpos] + interp*npsdata[mtfpos+1];
					
				if(mtf<=.001)
					mtf = 0.001;
		
				if(nps<=.001)
					nps=.001;
			
				double snr = nps/(100*mtf);
					
				// Until i can fix this bit...
				//snr = 0.2;

				//if(Index==0)
				//{
				//	DM::Result(t_to_string<float>(mtf));
				//	DM::Result(t_to_string<float>(nps));
				//	DM::Result(t_to_string<float>(snr));
				//}

				//if(snr > 10 ||snr < 0.001)
				//{
				//	DM::Result("Ooops");
				//}

				
				double realnumerator = (Wminus[Index] + snr) * T[Index][0] - V[Index][0]*U[Index][0] - V[Index][1]*U[Index][1];
				double imaginarynumerator = (Wminus[Index] + snr) * T[Index][1] - V[Index][0]*U[Index][1] + V[Index][1]*U[Index][0];
				double denominator = ((Wminus[Index] + snr) * (W[Index] + snr)) - ( V[Index][0] * V[Index][0] + V[Index][1] * V[Index][1] );

				if(denominator==0)
				{
					DM::Result("Denominator fail");
				}
				
				// Factor of 2 :)
				Restored[Index][0] = realnumerator/denominator;
				Restored[Index][1] = imaginarynumerator/denominator;
			}
}

const char* makerestoredMTFNPSsource = 
"__kernel void clMakeRestoredMTFNPS(__global double* clW,__global double* clWminus, __global double2* clV,__global double2* clT, __global double2* clU, __global double2* clRestored, int sizeX, int sizeY, __global float* clMTF, __global float* clNPS, int mtflength, float scaleMTFx, float scaleMTFy)	\n"
"{											\n"
"	//Get the work items ID					\n"
"	int xid = get_global_id(0);				\n"
"	int yid = get_global_id(1);				\n"
"											\n"
"	if(xid<sizeX&&yid<sizeY)				\n"
"	{										\n"
"		int Index = xid + yid*sizeX;		\n"
"											\n"
"		int xposmtf;						\n"
"		int yposmtf;						\n"
"											\n"
"		if(xid < sizeX/2)					\n"
"			xposmtf = xid;					\n"
"		else								\n"
"			xposmtf = (sizeX-xid);			\n"
"											\n"
"		if(yid < sizeY/2)					\n"
"			yposmtf = yid;					\n"
"		else								\n"
"			yposmtf = sizeY-yid;			\n"
"											\n"
"		float mtfposition = sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf);			\n"
"																														\n"
"		int mtfpos = floor(sqrt(scaleMTFy*yposmtf*scaleMTFy*yposmtf + scaleMTFx*xposmtf*scaleMTFx*xposmtf));			\n"
"																														\n"		
"		float interp = mtfposition - floor(mtfposition);																\n"
"																														\n"		
"		if(mtfpos>mtflength-2)																							\n"
"			mtfpos=mtflength-2;																							\n"
"																														\n"
"		float mtf = (1-interp)*clMTF[mtfpos] + interp*clMTF[mtfpos+1];													\n"
"		if(mtf < 0.001)																									\n"
"			mtf = 0.001;																								\n"
"		float nps = (1-interp)*clNPS[mtfpos] + interp*clNPS[mtfpos+1];													\n"
"		if(nps < 0.001)																									\n"
"			nps = 0.001;																								\n"
"		float snr = nps/(100*mtf);																						\n"
"		double real = (clWminus[Index] + snr) * clT[Index].x - clV[Index].x*clU[Index].x - clV[Index].y*clU[Index].y;	\n"
"		double imag = (clWminus[Index] + snr) * clT[Index].y - clV[Index].x*clU[Index].y + clV[Index].y*clU[Index].x;	\n"
"		double denom = ((clWminus[Index] + snr) * (clW[Index] + snr)) - ( clV[Index].x * clV[Index].x + clV[Index].y * clV[Index].y ); \n"
"		clRestored[Index].x = real/denom;	\n"
"		clRestored[Index].y = imag/denom;	\n"
"	}	\n"
"}	\n"
;

const char* getQsource = 
"__kernel void clCalculateQ(__global double* clW, __global double* clWminus, __global double2* clV, __global double2* clQ, int sizeX, int sizeY)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double wtf  = (clWminus[Index]+0.2)*clW[Index] -clV[Index].x*clV[Index].x - clV[Index].y*clV[Index].y;\n"
"		double wtfCr = 0.2*clV[Index].x;\n"
"		double wtfCi = 0.2*clV[Index].y;\n"
"		double denom = (clWminus[Index]+0.2)*(clW[Index]+0.2)-clV[Index].x*clV[Index].x - clV[Index].y*clV[Index].y; \n"
"		clQ[Index].x = wtf*wtf/(denom*denom) - wtfCr*wtfCr/(denom*denom) - wtfCi*wtfCi/(denom*denom);\n"
"		clQ[Index].y = 0;\n"
"	}\n"
"}	\n"
;

const char* minuswavefunctionsource = 
"__kernel void clMinusWavefunction( __global double2* clFFT, __global double2* clFFTminus, int sizeX, int sizeY)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		int Index2 = sizeX-xid + (sizeY-yid)*sizeX; \n"
"		if(xid==0) \n"
"			Index2-=sizeX; \n"
"		if(yid==0) \n"
"			Index2-=sizeX*sizeY; \n"
"		clFFTminus[Index] = clFFT[Index2]; \n"
"	}\n"
"}	\n"
;

const char* getPCIsource = 
"__kernel void clCalculatePCI(__global double2* clQ, __global double2* clFFT, __global double2* clFFTminus, \n"
"__global double* CLxFrequencies, __global double* CLyFrequencies, __global double2* clPCI, __global double2* clPCIC, \n"
" __global double2* clPCIM, __global double2* clPCIS, int sizeX, int sizeY, float C1, float C3, float A1r, float A1i, float wavelength, float objap)	\n"
"{	\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	\n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		double freq = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]);	\n"
"		double kx = -CLxFrequencies[xid];	\n"
"		double ky = -CLyFrequencies[yid];	\n"
"		double gamma = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx * kx + 6 * C1 * kx * kx - 6 * A1r * ky * ky	\n"
"				  + 6 * C1 * ky * ky + 3 * C3 * wavelength * wavelength * kx * kx * kx * kx + 6 * C3 * wavelength * wavelength * kx * kx * ky * ky	\n"
"				  + 3 * C3 * wavelength * wavelength * ky * ky * ky * ky + 12 * A1i * kx * ky) ;	\n"
"		double arg1 = atan2(clFFT[Index].y,clFFT[Index].x);					\n"
"		double arg2 = atan2(clFFTminus[Index].y,clFFTminus[Index].x);		\n"
"		double angle = atan2(ky,kx);										\n"
"		clPCI[Index].x = (freq<=objap)*clQ[Index].x*-1 *cos(arg1+arg2+2*gamma); \n"
"		clPCI[Index].y = 0; \n"
"		clPCIC[Index].x = (freq<=objap)*clQ[Index].x*-1 *cos(arg1+arg2+2*gamma)*cos(2*angle); \n"
"		clPCIC[Index].y = 0; \n"
"		clPCIS[Index].x = (freq<=objap)*clQ[Index].x*-1 *cos(arg1+arg2+2*gamma)*sin(2*angle); \n"
"		clPCIS[Index].y = 0; \n"
"		clPCIM[Index].x = (freq<=objap)*(clPCI[Index].x + sqrt(clPCIC[Index].x * clPCIC[Index].x + clPCIS[Index].x * clPCIS[Index].x)); \n"
"		clPCIM[Index].y = 0; \n"
"	}\n"
"}	\n"
;

// MyDlg2 dialog

IMPLEMENT_DYNAMIC(MyDlg2, CDialog)

MyDlg2::MyDlg2(CWnd* pParent /*=NULL*/)
	: CDialog(MyDlg2::IDD, pParent)
	, m_mtfpath(_T(""))
	, m_npspath(_T(""))
	, m_infLimitVal(_T("8"))
	, m_fstartVal(_T("-100"))
	, m_fendVal(_T("100"))
	, m_fstepsVal(_T("50"))
	, m_SphericalVal(_T("0"))
{
	gotmtf=false;
	gotnps=false;
	OpenCLAvailable = false;
	context = NULL;
	numDevices = 0;
	devices = NULL;


	// Maybe Can Do OpenCL setup and device registering here - Print to Ouput with device data?

	

	// Discover and initialize available platforms
	cl_uint numPlatforms = 0;
	cl_platform_id * platforms = NULL;

	// Use clGetPlatformIds() to retrieve the number of platforms
	status = clGetPlatformIDs(0,NULL,&numPlatforms);

	// Allocate enough space for each platform
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));

	// Fill in platforms with clGetPlatformIDs()
	status = clGetPlatformIDs(numPlatforms,platforms,NULL);

	// Discover and initialize available devices

	
	// use clGetDeviceIDs() to retrieve number of devices present
	status = clGetDeviceIDs(platforms[PLATFORM],CL_DEVICE_TYPE_ALL,0,NULL,&numDevices);

	// Allocate enough space for each device
	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));

	// Fill in devices with clGetDeviceIDs()
	status = clGetDeviceIDs(platforms[PLATFORM],CL_DEVICE_TYPE_ALL,numDevices,devices,NULL);

	

	// Create a context using clCreateContext() and associate with the devices
	context = clCreateContext(NULL,numDevices,devices,NULL,NULL,&status);

	

	// Create a command queue using clCreateCommandQueue(),
	// and associate it with the device you want to run on
	cmdQueue = clCreateCommandQueue(context,devices[0],0,&status);

	// Most of initialisation is done, would be nice to print device information...
	//Getting the device name

	size_t deviceNameLength = 4096;
	size_t actualSize;
	char* tempDeviceName = (char*)malloc(4096);
	char* deviceName;

	status |= clGetDeviceInfo(devices[0], CL_DEVICE_NAME, deviceNameLength, tempDeviceName, &actualSize);

	if(status == CL_SUCCESS)
	{
		deviceName = (char*)malloc(actualSize);
		memcpy(deviceName, tempDeviceName, actualSize);
		free(tempDeviceName);
	}


	if(status!=CL_SUCCESS)
	{
		DM::Result("Could not setup OpenCL on this computer, run clinfo to check availability");
	}
	else
	{
		std::string devName(deviceName);
		DM::Result("Using OpenCL on device "+devName);

		OpenCLAvailable = true;
	}


}

MyDlg2::~MyDlg2()
{
}

void MyDlg2::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT2, m_mtfpath);
	DDX_Text(pDX, IDC_EDIT1, m_npspath);
	DDX_Control(pDX, IDC_EDIT2, m_mtflocation);
	DDX_Control(pDX, IDC_EDIT1, m_npslocation);
	DDX_Control(pDX, IDC_EDIT6, m_infLimit);
	DDX_Control(pDX, IDC_EDIT5, m_fstart);
	DDX_Control(pDX, IDC_EDIT4, m_fend);
	DDX_Control(pDX, IDC_EDIT11, m_Spherical);
	DDX_Text(pDX, IDC_EDIT6, m_infLimitVal);
	DDX_Text(pDX, IDC_EDIT5, m_fstartVal);
	DDX_Text(pDX, IDC_EDIT4, m_fendVal);
	DDX_Text(pDX, IDC_EDIT3, m_fstepsVal);
	DDX_Text(pDX, IDC_EDIT11, m_SphericalVal);
}


BEGIN_MESSAGE_MAP(MyDlg2, CDialog)
	ON_BN_CLICKED(IDC_BUTTON1, &MyDlg2::OnBnClickedButton1Alt)
	ON_BN_CLICKED(IDC_BUTTON2, &MyDlg2::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &MyDlg2::OnBnClickedButton3)
END_MESSAGE_MAP()


// MyDlg2 message handlers

// MTF Button
void MyDlg2::OnBnClickedButton2()
{
	// TODO: Add your control notification handler code here
	DM::String mtflocation;
	DM::OpenDialog(mtflocation);

	//Check file exists
	if(!DM::DoesFileExist(mtflocation))
	{
		DM::OpenAndSetProgressWindow("error","loading","file");
		return;
	}

	mtffile = DM::NewImageFromFile(mtflocation);
	gotmtf = true;


	// These don't seem to be working to display the text...
	std::string mtfpath = mtflocation;
	//m_mtfpath = mtfpath.c_str();

	m_mtflocation.SetWindowText(mtfpath.c_str());
	//UpdateData(TRUE);

	// Check it is at least loading the file
	mtffile.GetOrCreateImageDocument().Show();

}

// NPS Button
void MyDlg2::OnBnClickedButton3()
{
	// TODO: Add your control notification handler code here
	DM::String npslocation;
	DM::OpenDialog(npslocation);

	//Check file exists
	if(!DM::DoesFileExist(npslocation))
	{
		DM::OpenAndSetProgressWindow("error","loading","file");
		return;
	}

	npsfile = DM::NewImageFromFile(npslocation);
	gotnps = true;

	// These don't seem to be working to display the text...
	std::string npspath = npslocation;
	//m_npspath = npspath.c_str();

	m_npslocation.SetWindowText(npspath.c_str());
	//UpdateData(TRUE);

	// Check it is at least loading the file
	npsfile.GetOrCreateImageDocument().Show();
}





// Alternate Version for any number of images 
void MyDlg2::OnBnClickedButton1Alt()
{
	// At this point we want the front image to be a stack of any size;
	DM::Image frontImage = DM::GetFrontImage();

	// Get dimensions so we know how big to make the slices;
	unsigned long xDim;
	unsigned long yDim;
	unsigned long numberOfImages;

	frontImage.GetDimensionSizes(xDim,yDim,numberOfImages);

	// Check for ROI and if not create one.

	int numROI;
	DM::ImageDisplay imDisp;
	imDisp = frontImage.GetImageDisplay(0);
	
	numROI = imDisp.CountROIs();

	DM::ROI regionOfInterest;

	if(numROI==0)
	{
		// Create image display
		regionOfInterest = DM::NewROI();
		regionOfInterest.SetRectangle(0,0,yDim,xDim);
		imDisp.AddROI(regionOfInterest);
	}

	else
	{
		regionOfInterest = imDisp.GetROI(0);
	}

	// Get Region Positions
	float top;
	float left;
	float bottom;
	float right;

	// Check its a rectangle ROI
	if(!regionOfInterest.IsRectangle())
	{
		DM::OpenAndSetProgressWindow("ROI is not a rectangle","","");
		return;
	}
	

	regionOfInterest.GetRectangle(&top,&left,&bottom,&right);

	// Do Some Checks
	int iTop = round(top);
	int iLeft = round(left);
	int iBottom = round(bottom);
	int iRight = round(right);

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	// Check for power of 2
	//TODO: Add padding option to do other sizes
	if(!(sizeX==64||sizeX==128||sizeX==256||sizeX==512||sizeX==1024||sizeX==2048||sizeX==4096))
	{
		DM::Result("ROI must be power of 2 size");
		return;
	}

	if(!(sizeY==64||sizeY==128||sizeY==256||sizeY==512||sizeY==1024||sizeY==2048||sizeY==4096))
	{
		DM::Result("ROI must be power of 2 size");
		return;
	}


	// Copy slice data into seperate images
	Gatan::PlugIn::ImageDataLocker seriesLocker(frontImage);

	float* seriesdata = (float*) seriesLocker.get();

	// Need pixelscale, and voltage, and focalstep to get right parameters.
	DM::TagGroup imagetags = frontImage.GetTagGroup();

	float fStep;
	imagetags.GetTagAsFloat("Focal Series:Adjusted focalstep",&fStep);

	float voltage;
	imagetags.GetTagAsFloat("Microscope Info:Voltage",&voltage);

	float pixelscale = frontImage.GetDimensionScale(0);
	
	// Setup two frequency arrays for this size.

	double* xFrequencies = new double[sizeX];
	double* yFrequencies = new double[sizeY];

	int midX = ceil(float(sizeX)/2);
	int midY = ceil(float(sizeY)/2);
	
	for(int i = 1 ; i <= sizeX ; i++)
	{
		if(i <= midX)
			xFrequencies[i-1] = (i - 1)/(pixelscale * sizeX);
		else xFrequencies[i-1] = (i - 1 - sizeX)/(pixelscale * sizeX);
	}

	for(int i = 1 ; i <= sizeY ; i++)
	{
		if(i <= midY)
			yFrequencies[i-1] = (i - 1)/(pixelscale * sizeY);
		else yFrequencies[i-1] = (i - 1 - sizeY)/(pixelscale * sizeY);
	}

	// Calculate Electron Wavelength in nm

	float e = 1.6e-19;
	float wavelength = 6.63e-34*3e+8*1e+9 / sqrt((e*voltage*(2*9.11e-31*9e+16 + e*voltage)));

	
	// Use a seperate function to take seriesdata and 2 image positions and return an x and y shift value

	int* xShiftVals = new int[numberOfImages];
	int* yShiftVals = new int[numberOfImages];

	xShiftVals[0] = 0;
	yShiftVals[0] = 0;

	clAmdFftStatus fftStatus;
	clAmdFftSetupData fftSetupData;

	clAmdFftInitSetupData(&fftSetupData);
	
	fftStatus = clAmdFftSetup(&fftSetupData);

	cl_event outEvent = NULL;

	//	Local Data
	size_t buffSizeBytesIn = 0;
	size_t buffSizeBytesOut = 0;
	size_t fftVectorSize= 0;
	size_t fftVectorSizePadded = 0;
	size_t fftBatchSize = 0;
	cl_uint nBuffersOut = 0;
	cl_uint profileCount = 0;


	clAmdFftPlanHandle fftplan;
	clAmdFftDim fftdim = CLFFT_2D;
	clAmdFftResultLocation	place = CLFFT_OUTOFPLACE;
	clAmdFftLayout inLayout  = CLFFT_COMPLEX_INTERLEAVED;
	clAmdFftLayout outLayout = CLFFT_COMPLEX_INTERLEAVED;

	size_t clLengths[ 3 ];
	size_t clPadding[ 3 ] = {0, 0, 0 };  // *** TODO
	size_t clStrides[ 4 ];
	size_t batchSize = 1;


	clLengths[0]=sizeX;
	clLengths[1]=sizeY;
	clLengths[2]=1;

	clStrides[ 0 ] = 1;
	clStrides[ 1 ] = clStrides[ 0 ] * (clLengths[ 0 ] + clPadding[ 0 ]);
	clStrides[ 2 ] = clStrides[ 1 ] * (clLengths[ 1 ] + clPadding[ 1 ]);
	clStrides[ 3 ] = clStrides[ 2 ] * (clLengths[ 2 ] + clPadding[ 2 ]);

	fftVectorSize	= clLengths[ 0 ] * clLengths[ 1 ] * clLengths[ 2 ];
	fftVectorSizePadded = clStrides[ 3 ];
	fftBatchSize	= fftVectorSizePadded * batchSize;

	// Probably always just use complex interleaved
	nBuffersOut      = 1;
	buffSizeBytesOut = fftBatchSize * sizeof( std::complex< double > );


	// Setup cl_mem device buffers
	cl_mem fftInputBuffer = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem fftOutputBuffer = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem dataBuffer1 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem dataBuffer2 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);


	fftStatus = clAmdFftCreateDefaultPlan( &fftplan, context, fftdim, clLengths );


	//	Default plan creates a plan that expects an inPlace transform with interleaved complex numbers
	fftStatus = clAmdFftSetResultLocation( fftplan, place );
	fftStatus = clAmdFftSetPlanPrecision(fftplan,CLFFT_DOUBLE);
	fftStatus = clAmdFftSetLayout( fftplan, inLayout, outLayout );
	fftStatus = clAmdFftSetPlanBatchSize( fftplan, batchSize );
	fftStatus = clAmdFftSetPlanScale (fftplan, CLFFT_FORWARD, 1/sqrtf(sizeX*sizeY));
	fftStatus = clAmdFftSetPlanScale (fftplan, CLFFT_BACKWARD, 1/sqrtf(sizeX*sizeY));


	// Not using padding here yet
	if ((clPadding[ 0 ] | clPadding[ 1 ] | clPadding[ 2 ]) != 0) {
		clAmdFftSetPlanInStride  ( fftplan, fftdim, clStrides );
		clAmdFftSetPlanOutStride ( fftplan, fftdim, clStrides );
		clAmdFftSetPlanDistance  ( fftplan, clStrides[ fftdim ], clStrides[ fftdim ]);
	}

	fftStatus =clAmdFftBakePlan( fftplan, 1, &cmdQueue, NULL, NULL );
	
	//get the buffersize
	size_t buffersize=0;
	fftStatus = clAmdFftGetTmpBufSize(fftplan, &buffersize );
		
	//allocate the intermediate buffer	
	cl_mem clMedBuffer=NULL;
		
	if (buffersize)
	{
		cl_int medstatus;
		clMedBuffer = clCreateBuffer ( context, CL_MEM_READ_WRITE, buffersize, 0, &medstatus);
	}
	
	// Upload frequencies to GPU
	cl_mem CLxFrequencies = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeof(double), 0, &status);
	cl_mem CLyFrequencies = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeY * sizeof(double), 0, &status);

	clEnqueueWriteBuffer( cmdQueue, CLxFrequencies, CL_TRUE, 0, sizeX*sizeof(double), &xFrequencies[ 0 ], 
				0, NULL, NULL );
	clEnqueueWriteBuffer( cmdQueue, CLyFrequencies, CL_TRUE, 0, sizeY*sizeof(double), &yFrequencies[ 0 ], 
				0, NULL, NULL );


	// Now make the kernel for PCPCF

	// Create a program using clCreateProgramWithSource()
	const char* pcpcfkernel = pcpcfSource;

	cl_program program = clCreateProgramWithSource(context,1,&pcpcfkernel,NULL,&status);

	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting kernel code");
	}

	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(program,numDevices,devices,NULL,NULL,NULL);

	cl_kernel clPCPCF;

	// Use clCreateKernel() to create a kernel from the function declared earlier
	clPCPCF = clCreateKernel(program,"clPCPCF",&status);

	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling PCPCF kernel");
	}


	// Copy data to relevant buffers and calculate pcpcf
	for(int i = 1 ; i < numberOfImages ; i++)
	{
		int xShiftAmount;
		int yShiftAmount;

		OpenCLPCPCF(seriesdata,i-1,i,cmdQueue,fftplan,fftInputBuffer,fftOutputBuffer,dataBuffer1,dataBuffer2,clMedBuffer,clPCPCF,xShiftAmount,yShiftAmount,sizeX,sizeY,CLxFrequencies,CLyFrequencies,wavelength,fStep, iLeft,iTop,xDim,yDim);

		xShiftVals[i] = xShiftAmount;
		yShiftVals[i] = yShiftAmount;
	}


	// Crop and calculate reconstruction

	// Now we have shifts can determine if ROI region is available in all images, if not modify ROI.

	int maxShiftx = 0;
	int maxShifty = 0;
	int minShiftx = 0;
	int minShifty = 0;

	int* xShifts = new int[numberOfImages];
	int* yShifts = new int[numberOfImages];

	int sumXShift = 0;
	int sumYShift = 0;

	for(int i = 0; i <numberOfImages ; i++)
	{
		sumXShift+= xShiftVals[i];
		sumYShift +=yShiftVals[i];
		xShifts[i] = sumXShift;
		yShifts[i] = sumYShift;
	}

	for(int i = 0; i <numberOfImages ; i++)
	{
		if(xShifts[i] > maxShiftx)
			maxShiftx = xShifts[i];

		if(xShifts[i] < minShiftx)
			minShiftx = xShifts[i];

		if(yShifts[i] > maxShifty)
			maxShifty = yShifts[i];

		if(yShifts[i] < minShifty)
			minShifty = yShifts[i];

	}

	// Test if ROI is in ok place. (Shifts could be other way round)

	if(iLeft < -minShiftx)
		iLeft = -minShiftx;

	if(xDim-iRight < maxShiftx)
		iRight = xDim - maxShiftx;

	if(iTop < -minShifty)
		iTop = -minShifty;

	if(yDim-iBottom < maxShifty)
		iBottom = yDim - maxShifty;

	sizeX = iRight - iLeft;
	sizeY = iBottom - iTop;

	// Make sure can still do power of two size
	// Check for power of 2
	//TODO: Add padding option to do other sizes
	if(!(sizeX==64||sizeX==128||sizeX==256||sizeX==512||sizeX==1024||sizeX==2048||sizeX==4096))
	{
		DM::Result("ROI area not available in all images after drift correction");
		return;
	}

	if(!(sizeY==64||sizeY==128||sizeY==256||sizeY==512||sizeY==1024||sizeY==2048||sizeY==4096))
	{
		DM::Result("ROI area not available in all images after drift correction");
		return;
	}

	

	// Cleanup First wave of Stuff

	delete[] xFrequencies;
	delete[] yFrequencies;

	// OpenCL cleanup

	clReleaseKernel(clPCPCF);
	clReleaseProgram(program);
	
	// Dont clean frequencies or buffers yes, usefull later.

	// Now cut shift corrected ROI into a stack image ready for restoration
	
	DM::Image driftCorrected = DM::RealImage("Drift Corrected Stack",4,sizeX,sizeY,numberOfImages);
	Gatan::PlugIn::ImageDataLocker driftCorrectedLocker(driftCorrected);
	float* driftCorrecteddata = (float*) driftCorrectedLocker.get();

	for(int image = 0; image < numberOfImages; image++ )
		for(int j = 0; j < sizeY; j++)
			for(int i = 0; i < sizeX; i++)
			{
				int offsetX = xShifts[image];
				int offsetY = yShifts[image];
				driftCorrecteddata[image*sizeX*sizeY + j*sizeX + i] = seriesdata[image*xDim*yDim + (j+offsetY+iTop)*xDim + i + iLeft + offsetX];
			}

	driftCorrectedLocker.~ImageDataLocker();
	driftCorrected.SetDimensionCalibration(0,0,pixelscale,"nm",0);
	driftCorrected.SetDimensionCalibration(1,0,pixelscale,"nm",0);
	driftCorrected.GetOrCreateImageDocument().Show();

	// Now we can do a restoration from this driftcorrected series
	// Need to remake frequencies first incase sizes have changed.

	/* Sizes cant change in OpenCL version

	float* xFrequencies2 = new float[sizeX];
	float* yFrequencies2 = new float[sizeY];

	int midX2 = ceil(float(sizeX)/2);
	int midY2 = ceil(float(sizeY)/2);
	
	for(int i = 1 ; i <= sizeX ; i++)
	{
		if(i <= midX2)
			xFrequencies2[i-1] = (i - 1)/(pixelscale * sizeX);
		else xFrequencies2[i-1] = (i - 1 - sizeX)/(pixelscale * sizeX);
	}

	for(int i = 1 ; i <= sizeY ; i++)
	{
		if(i <= midY2)
			yFrequencies2[i-1] = (i - 1)/(pixelscale * sizeY);
		else yFrequencies2[i-1] = (i - 1 - sizeY)/(pixelscale * sizeY);
	}

	*/



	// To do a restoration need to setup the W,W-,T,U,V arrays before hand, then can add each image seperately.

	cl_mem clW = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( double ), 0, &status);
	cl_mem clWminus = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( double ), 0, &status);
	cl_mem clw = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clwminus = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clT = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clU = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clV = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clRestored = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);
	cl_mem clRestoredMinus = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< double > ), 0, &status);

	// Get data from GUI (refresh it first).
	UpdateData(TRUE);

	float fstartrange = boost::lexical_cast<float>(m_fstartVal);
	float fendrange = boost::lexical_cast<float>(m_fendVal);
	int numsteps = boost::lexical_cast<int>(m_fstepsVal);
	float inflimit = boost::lexical_cast<float>(m_infLimitVal);
	float C3 = boost::lexical_cast<float>(m_SphericalVal)*1000;

	float searchstep = (fendrange-fstartrange)/(numsteps-1);

	// NOTE - these shouldnt be fixed, should be options if doing recon PCF.
	float beta = 0.0002; //rad
	float delta = 3; //nm
	float normalise = sqrtf(sizeX*sizeY);
	float objap = inflimit;

	// Used in MTF/NPS

	// Want to work out required adjustment in scale for mtf in x,y directions based on image size.
	// NOTE: This only works for mtf of size 725 from 1024^2 images...
	// Could get 1024 by mtflength*2/sqrt(2)
	float scaleMTFx = 1024.0f/float(sizeX);
	float scaleMTFy = 1024.0f/float(sizeY);

	//DM::Result(t_to_string(scaleMTFx));
	
	float* mtfdata;
	float* npsdata;
	int mtflength;
	int npslength;

	Gatan::PlugIn::ImageDataLocker mtfLocker;
	Gatan::PlugIn::ImageDataLocker npsLocker;

	cl_mem clMTF;
	cl_mem clNPS;

	if(gotmtf)
	{
			mtfLocker = Gatan::PlugIn::ImageDataLocker(mtffile);
			mtfdata = (float*) mtfLocker.get();
			mtflength = mtffile.GetDimensionSize(0);
			clMTF = clCreateBuffer ( context, CL_MEM_READ_WRITE, mtflength * sizeof( float ), 0, &status);
			clEnqueueWriteBuffer( cmdQueue, clMTF, CL_TRUE, 0, mtflength*sizeof(float), &mtfdata[ 0 ], 
				0, NULL, NULL );
	}

	if(gotnps)
	{
			npsLocker = Gatan::PlugIn::ImageDataLocker(npsfile);
			npsdata = (float*) npsLocker.get();
			npslength = npsfile.GetDimensionSize(0);
			clNPS = clCreateBuffer ( context, CL_MEM_READ_WRITE, mtflength * sizeof( float ), 0, &status);
			clEnqueueWriteBuffer( cmdQueue, clNPS, CL_TRUE, 0, mtflength*sizeof(float), &npsdata[ 0 ], 
				0, NULL, NULL );
	}


	mtfLocker.~ImageDataLocker();
	npsLocker.~ImageDataLocker();

	std::vector< std::complex< double > > dataOne( sizeX*sizeY );
	std::vector< std::complex< double > > dataTwo( sizeX*sizeY );

	// Setup all kernels for Restoration

	const char* wavetransferfunctionkernel = wavetransferfunctionsource;
	cl_program wavetransferfunctionprogram = clCreateProgramWithSource(context,1,&wavetransferfunctionkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wtf kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wavetransferfunctionprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwtf;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwtf = clCreateKernel(wavetransferfunctionprogram,"clWaveTransferFunction",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wtf kernel");
	}

	const char* wavetransferfunctionminuskernel = wavetransferfunctionminussource;
	cl_program wavetransferfunctionminusprogram = clCreateProgramWithSource(context,1,&wavetransferfunctionminuskernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wtf minus kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wavetransferfunctionminusprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwtfminus;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwtfminus = clCreateKernel(wavetransferfunctionminusprogram,"clWaveTransferFunctionMinus",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wtf minus kernel");
	}

	const char* wavetransferfunctionmtfkernel = wavetransferfunctionmtfsource;
	cl_program wavetransferfunctionmtfprogram = clCreateProgramWithSource(context,1,&wavetransferfunctionmtfkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wtf mtf kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wavetransferfunctionmtfprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwtfmtf;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwtfmtf = clCreateKernel(wavetransferfunctionmtfprogram,"clWaveTransferFunctionMTF",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wtf mtf kernel");
	}

	const char* wavetransferfunctionminusmtfkernel = wavetransferfunctionminusmtfsource;
	cl_program wavetransferfunctionminusmtfprogram = clCreateProgramWithSource(context,1,&wavetransferfunctionminusmtfkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wtf minus mtf kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wavetransferfunctionminusmtfprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwtfminusmtf;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwtfminusmtf = clCreateKernel(wavetransferfunctionminusmtfprogram,"clWaveTransferFunctionMinusMTF",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wtf minus mtf kernel");
	}

	const char* wienerwkernel = wienerwsource;
	cl_program wienerwprogram = clCreateProgramWithSource(context,1,&wienerwkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wienerW kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wienerwprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwienerw;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwienerw = clCreateKernel(wienerwprogram,"clWienerW",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wienerW kernel");
	}

	const char* wienerwminuskernel = wienerwsource;
	cl_program wienerwminusprogram = clCreateProgramWithSource(context,1,&wienerwminuskernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wienerW minus kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wienerwminusprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwienerwminus;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwienerwminus = clCreateKernel(wienerwminusprogram,"clWienerW",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wienerW minus kernel");
	}

	const char* wienervkernel = wienervsource;
	cl_program wienervprogram = clCreateProgramWithSource(context,1,&wienervkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wienerV kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wienervprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwienerv;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwienerv = clCreateKernel(wienervprogram,"clWienerV",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wienerV kernel");
	}

	const char* wienertkernel = wienertsource;
	cl_program wienertprogram = clCreateProgramWithSource(context,1,&wienertkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wienerT kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wienertprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwienert;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwienert = clCreateKernel(wienertprogram,"clWienerT",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wienerT kernel");
	}

	const char* wienerukernel = wienerusource;
	cl_program wieneruprogram = clCreateProgramWithSource(context,1,&wienerukernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting wienerU kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(wieneruprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clwieneru;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clwieneru = clCreateKernel(wieneruprogram,"clWienerU",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling wienerU kernel");
	}

	const char* makerestoredkernel = makerestoredsource;
	cl_program makerestoredprogram = clCreateProgramWithSource(context,1,&makerestoredkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting makeRestored kernel");
	}

	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(makerestoredprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clmakerestored;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clmakerestored = clCreateKernel(makerestoredprogram,"clMakeRestored",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling makeRestored kernel");
	}

	const char* makerestoredMTFNPSkernel = makerestoredMTFNPSsource;
	cl_program makerestoredMTFNPSprogram = clCreateProgramWithSource(context,1,&makerestoredMTFNPSkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting makeRestored mtfnps kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(makerestoredMTFNPSprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clmakerestoredMTFNPS;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clmakerestoredMTFNPS = clCreateKernel(makerestoredMTFNPSprogram,"clMakeRestoredMTFNPS",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling makeRestored mtfnps kernel");
	}


	if(gotmtf&&gotnps)
	{
		// Set Kernel Arguments
		status = clSetKernelArg(clmakerestoredMTFNPS,0,sizeof(cl_mem),&clW);
		status |= clSetKernelArg(clmakerestoredMTFNPS,1,sizeof(cl_mem),&clWminus);
		status |= clSetKernelArg(clmakerestoredMTFNPS,2,sizeof(cl_mem),&clV);
		status |= clSetKernelArg(clmakerestoredMTFNPS,3,sizeof(cl_mem),&clT);
		status |= clSetKernelArg(clmakerestoredMTFNPS,4,sizeof(cl_mem),&clU);
		status |= clSetKernelArg(clmakerestoredMTFNPS,5,sizeof(cl_mem),&clRestored);
		status |= clSetKernelArg(clmakerestoredMTFNPS,6,sizeof(int),&sizeX);
		status |= clSetKernelArg(clmakerestoredMTFNPS,7,sizeof(int),&sizeY);
		status |= clSetKernelArg(clmakerestoredMTFNPS,8,sizeof(cl_mem),&clMTF);
		status |= clSetKernelArg(clmakerestoredMTFNPS,9,sizeof(cl_mem),&clNPS);
		status |= clSetKernelArg(clmakerestoredMTFNPS,10,sizeof(int),&mtflength);
		status |= clSetKernelArg(clmakerestoredMTFNPS,11,sizeof(float),&scaleMTFx);
		status |= clSetKernelArg(clmakerestoredMTFNPS,12,sizeof(float),&scaleMTFy);
	}
	else
	{
		// Set Kernel Arguments
		status = clSetKernelArg(clmakerestored,0,sizeof(cl_mem),&clW);
		status |= clSetKernelArg(clmakerestored,1,sizeof(cl_mem),&clWminus);
		status |= clSetKernelArg(clmakerestored,2,sizeof(cl_mem),&clV);
		status |= clSetKernelArg(clmakerestored,3,sizeof(cl_mem),&clT);
		status |= clSetKernelArg(clmakerestored,4,sizeof(cl_mem),&clU);
		status |= clSetKernelArg(clmakerestored,5,sizeof(cl_mem),&clRestored);
		status |= clSetKernelArg(clmakerestored,6,sizeof(int),&sizeX);
		status |= clSetKernelArg(clmakerestored,7,sizeof(int),&sizeY);
	}


	Gatan::PlugIn::ImageDataLocker driftCorrectedLocker2(driftCorrected);
	float* driftCorrecteddata2 = (float*) driftCorrectedLocker2.get();

	// Process each Image
	for(int image = 0 ; image < numberOfImages ; image++)
	{
		// Images need to be normalised before FFT for FSR
		float sum = GetSum(driftCorrecteddata2,sizeX*sizeY,image*sizeX*sizeY);
		

		// I'm not sure if I should still be using driftcorrecteddata after destroying image locker?
		for(int j = 0 ; j < sizeY ; j++)
			for(int i = 0 ; i < sizeX ; i++)
			{
				dataOne[i+j*sizeX] = (sizeX*sizeY*driftCorrecteddata2[image*sizeX*sizeY + j*sizeX + i]/(sum)) - 1;		
			}

		clEnqueueWriteBuffer( cmdQueue, fftInputBuffer, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &dataOne[ 0 ], 
				0, NULL, NULL );

		fftStatus = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&fftInputBuffer, &dataBuffer1, clMedBuffer );
		
	
		// Now calculate w,wminus and add to W,W-,T,U,V
		// Because of image order first is most underfocus

		// Changed to try do first restore in plane of last image (usually around focus :))
		float defocus = image*fStep - numberOfImages*fStep;
		float zero = 0;
		int init = 0;

		if(image==0)
		{
			init = 1;

			if(gotmtf)
			{
				// Set Kernel Arguments
				status = clSetKernelArg(clwtfmtf,0,sizeof(cl_mem),&clw);
				status |= clSetKernelArg(clwtfmtf,1,sizeof(cl_mem),&CLxFrequencies);
				status |= clSetKernelArg(clwtfmtf,2,sizeof(cl_mem),&CLyFrequencies);
				status |= clSetKernelArg(clwtfmtf,3,sizeof(int),&sizeX);
				status |= clSetKernelArg(clwtfmtf,4,sizeof(int),&sizeY);
				status |= clSetKernelArg(clwtfmtf,5,sizeof(float),&wavelength);
				status |= clSetKernelArg(clwtfmtf,6,sizeof(float),&beta);
				status |= clSetKernelArg(clwtfmtf,7,sizeof(float),&delta);
				status |= clSetKernelArg(clwtfmtf,8,sizeof(float),&zero);
				status |= clSetKernelArg(clwtfmtf,9,sizeof(float),&zero);

				status |= clSetKernelArg(clwtfmtf,11,sizeof(float),&C3);
				status |= clSetKernelArg(clwtfmtf,12,sizeof(float),&inflimit);
				status |= clSetKernelArg(clwtfmtf,13,sizeof(cl_mem),&clMTF);
				status |= clSetKernelArg(clwtfmtf,14,sizeof(int),&mtflength);
				status |= clSetKernelArg(clwtfmtf,15,sizeof(float),&scaleMTFx);
				status |= clSetKernelArg(clwtfmtf,16,sizeof(float),&scaleMTFy);
			}
			else
			{
				// Set Kernel Arguments
				status = clSetKernelArg(clwtf,0,sizeof(cl_mem),&clw);
				status |= clSetKernelArg(clwtf,1,sizeof(cl_mem),&CLxFrequencies);
				status |= clSetKernelArg(clwtf,2,sizeof(cl_mem),&CLyFrequencies);
				status |= clSetKernelArg(clwtf,3,sizeof(int),&sizeX);
				status |= clSetKernelArg(clwtf,4,sizeof(int),&sizeY);
				status |= clSetKernelArg(clwtf,5,sizeof(float),&wavelength);
				status |= clSetKernelArg(clwtf,6,sizeof(float),&beta);
				status |= clSetKernelArg(clwtf,7,sizeof(float),&delta);
				status |= clSetKernelArg(clwtf,8,sizeof(float),&zero);
				status |= clSetKernelArg(clwtf,9,sizeof(float),&zero);

				status |= clSetKernelArg(clwtf,11,sizeof(float),&C3);
				status |= clSetKernelArg(clwtf,12,sizeof(float),&inflimit);
			}

			if(gotmtf)
			{
				// Set Kernel Arguments
				status = clSetKernelArg(clwtfminusmtf,0,sizeof(cl_mem),&clwminus);
				status |= clSetKernelArg(clwtfminusmtf,1,sizeof(cl_mem),&CLxFrequencies);
				status |= clSetKernelArg(clwtfminusmtf,2,sizeof(cl_mem),&CLyFrequencies);
				status |= clSetKernelArg(clwtfminusmtf,3,sizeof(int),&sizeX);
				status |= clSetKernelArg(clwtfminusmtf,4,sizeof(int),&sizeY);
				status |= clSetKernelArg(clwtfminusmtf,5,sizeof(float),&wavelength);
				status |= clSetKernelArg(clwtfminusmtf,6,sizeof(float),&beta);
				status |= clSetKernelArg(clwtfminusmtf,7,sizeof(float),&delta);
				status |= clSetKernelArg(clwtfminusmtf,8,sizeof(float),&zero);
				status |= clSetKernelArg(clwtfminusmtf,9,sizeof(float),&zero);
				
				status |= clSetKernelArg(clwtfminusmtf,11,sizeof(float),&C3);
				status |= clSetKernelArg(clwtfminusmtf,12,sizeof(float),&inflimit);
				status |= clSetKernelArg(clwtfminusmtf,13,sizeof(cl_mem),&clMTF);
				status |= clSetKernelArg(clwtfminusmtf,14,sizeof(int),&mtflength);
				status |= clSetKernelArg(clwtfminusmtf,15,sizeof(float),&scaleMTFx);
				status |= clSetKernelArg(clwtfminusmtf,16,sizeof(float),&scaleMTFy);
			}
			else
			{
				// Set Kernel Arguments
				status = clSetKernelArg(clwtfminus,0,sizeof(cl_mem),&clwminus);
				status |= clSetKernelArg(clwtfminus,1,sizeof(cl_mem),&CLxFrequencies);
				status |= clSetKernelArg(clwtfminus,2,sizeof(cl_mem),&CLyFrequencies);
				status |= clSetKernelArg(clwtfminus,3,sizeof(int),&sizeX);
				status |= clSetKernelArg(clwtfminus,4,sizeof(int),&sizeY);
				status |= clSetKernelArg(clwtfminus,5,sizeof(float),&wavelength);
				status |= clSetKernelArg(clwtfminus,6,sizeof(float),&beta);
				status |= clSetKernelArg(clwtfminus,7,sizeof(float),&delta);
				status |= clSetKernelArg(clwtfminus,8,sizeof(float),&zero);
				status |= clSetKernelArg(clwtfminus,9,sizeof(float),&zero);
				
				status |= clSetKernelArg(clwtfminus,11,sizeof(float),&C3);
				status |= clSetKernelArg(clwtfminus,12,sizeof(float),&inflimit);
			}

			// Set Kernel Arguments
			status = clSetKernelArg(clwienerw,0,sizeof(cl_mem),&clW);
			status |= clSetKernelArg(clwienerw,1,sizeof(cl_mem),&clw);
			status |= clSetKernelArg(clwienerw,2,sizeof(int),&sizeX);
			status |= clSetKernelArg(clwienerw,3,sizeof(int),&sizeY);
			

			// Set Kernel Arguments
			status = clSetKernelArg(clwienerwminus,0,sizeof(cl_mem),&clWminus);
			status |= clSetKernelArg(clwienerwminus,1,sizeof(cl_mem),&clwminus);
			status |= clSetKernelArg(clwienerwminus,2,sizeof(int),&sizeX);
			status |= clSetKernelArg(clwienerwminus,3,sizeof(int),&sizeY);
			

			// Set Kernel Arguments
			status = clSetKernelArg(clwienerv,0,sizeof(cl_mem),&clV);
			status |= clSetKernelArg(clwienerv,1,sizeof(cl_mem),&clw);
			status |= clSetKernelArg(clwienerv,2,sizeof(cl_mem),&clwminus);
			status |= clSetKernelArg(clwienerv,3,sizeof(int),&sizeX);
			status |= clSetKernelArg(clwienerv,4,sizeof(int),&sizeY);
			

			// Set Kernel Arguments
			status = clSetKernelArg(clwienert,0,sizeof(cl_mem),&clT);
			status |= clSetKernelArg(clwienert,1,sizeof(cl_mem),&clw);
			status |= clSetKernelArg(clwienert,2,sizeof(cl_mem),&dataBuffer1);
			status |= clSetKernelArg(clwienert,3,sizeof(int),&sizeX);
			status |= clSetKernelArg(clwienert,4,sizeof(int),&sizeY);
			

			// Set Kernel Arguments
			status = clSetKernelArg(clwieneru,0,sizeof(cl_mem),&clU);
			status |= clSetKernelArg(clwieneru,1,sizeof(cl_mem),&clwminus);
			status |= clSetKernelArg(clwieneru,2,sizeof(cl_mem),&dataBuffer1);
			status |= clSetKernelArg(clwieneru,3,sizeof(int),&sizeX);
			status |= clSetKernelArg(clwieneru,4,sizeof(int),&sizeY);
			
		}
		
		if(gotmtf)
		{
			// Only these two change
			status |= clSetKernelArg(clwtfmtf,10,sizeof(float),&defocus);
			status |= clSetKernelArg(clwtfminusmtf,10,sizeof(float),&defocus);
		}
		else
		{
			// Only these two change
			status |= clSetKernelArg(clwtf,10,sizeof(float),&defocus);
			status |= clSetKernelArg(clwtfminus,10,sizeof(float),&defocus);
		}

		// Oops, And this one!!
		status |= clSetKernelArg(clwienerw,4,sizeof(int),&init);
		status |= clSetKernelArg(clwienerwminus,4,sizeof(int),&init);
		status |= clSetKernelArg(clwienerv,5,sizeof(int),&init);
		status |= clSetKernelArg(clwienert,5,sizeof(int),&init);
		status |= clSetKernelArg(clwieneru,5,sizeof(int),&init);


		// Define and index space of work items for execution
		// A workgroup size is not required but can be used
		size_t* globalWorkSize = new size_t[2];
		
		// There are 'elements' work items
		globalWorkSize[0] = sizeX;
		globalWorkSize[1] = sizeY;


		// Messed up when I tried this
		size_t* localWorkSize = new size_t[2];
		localWorkSize[0] = 32;
		localWorkSize[1] = 16;

		if(gotmtf)
		{
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfmtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfminusmtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		}
		else
		{
			status = clEnqueueNDRangeKernel(cmdQueue,clwtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfminus,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		}

		status = clEnqueueNDRangeKernel(cmdQueue,clwienerw,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienerwminus,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienerv,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienert,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwieneru,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
	}

	clFinish(cmdQueue);

	// Can Free A Few Objects Now
	clReleaseMemObject(fftInputBuffer);
	clReleaseMemObject(fftOutputBuffer);
	// clReleaseMemObject(dataBuffer1); // need this
	driftCorrectedLocker2.~ImageDataLocker();
	
	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[2];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;


	if(gotmtf&&gotnps)
	{
		status = clEnqueueNDRangeKernel(cmdQueue,clmakerestoredMTFNPS,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	}
	else
	{
		status = clEnqueueNDRangeKernel(cmdQueue,clmakerestored,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	}

	status = clAmdFftEnqueueTransform( fftplan, CLFFT_BACKWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clRestored, &dataBuffer2, clMedBuffer ); 

	clFinish(cmdQueue);

	// Actually want FFT of this :)

	std::vector< std::complex< double > > HostRestored( sizeX*sizeY );
	clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &HostRestored[ 0 ], 
				0, NULL, NULL );


	/*
	if(gotmtf&&gotnps)
	{
		MakeRestoredWavefunctionmtfnps(W,Wminus,V,T,U,fft2,sizeX,sizeY,mtfdata,npsdata,mtflength,npslength,scaleMTFx,scaleMTFy);
	}
	else
	{
		MakeRestoredWavefunction(W,Wminus,V,T,U,fft2,sizeX,sizeY);
	}
	
	fftw_execute(fftplanInv2);
	*/
	
	DM::Image Restored = DM::ComplexImage("RestoredWave",16,sizeX,sizeY);

	Gatan::PlugIn::ImageDataLocker RestoredLocker(Restored);
	complex128* RestoredData = (complex128*) RestoredLocker.get();

	for(int j = 0 ; j < sizeY ; j++)
		for(int i = 0 ; i < sizeX ; i++)
			{
				RestoredData[i+j*sizeX].real() = HostRestored[i+j*sizeX].real()+1;
				RestoredData[i+j*sizeX].imag() = HostRestored[i+j*sizeX].imag();
			}

	
	RestoredLocker.~ImageDataLocker();
	Restored.SetDimensionCalibration(0,0,pixelscale,"nm",0);
	Restored.SetDimensionCalibration(1,0,pixelscale,"nm",0);
	Restored.GetOrCreateImageDocument().Show();
	Restored.GetImageDisplay(0).SetComplexMode(5);

	// More Cleanup
	seriesLocker.~ImageDataLocker();
	HostRestored.clear();
	
	// Now to calculate Q on GPU - might end up doing reductions on CPU til i make/find a kernel for it.
	// Or get sums by taking zero frequency component of FFT :D (Pretty wasteful - temporary measure)

	const char* getQkernel = getQsource;
	cl_program getQprogram = clCreateProgramWithSource(context,1,&getQkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting getQ kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(getQprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clgetQ;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clgetQ = clCreateKernel(getQprogram,"clCalculateQ",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling CalculateQ kernel");
	}

	const char* minuswavefunctionkernel = minuswavefunctionsource;
	cl_program minuswavefunctionprogram = clCreateProgramWithSource(context,1,&minuswavefunctionkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting minuswavefunction kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(minuswavefunctionprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clminuswavefunction;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clminuswavefunction = clCreateKernel(minuswavefunctionprogram,"clMinusWavefunction",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling minuswavefunction kernel");
	}

	const char* calculatePCIkernel = getPCIsource;
	cl_program calculatePCIprogram = clCreateProgramWithSource(context,1,&calculatePCIkernel,NULL,&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after getting calculatePCI kernel");
	}
	// Build the program for the devices with clBuildProgram()
	status = clBuildProgram(calculatePCIprogram,numDevices,devices,NULL,NULL,NULL);
	cl_kernel clcalculatePCI;
	// Use clCreateKernel() to create a kernel from the function declared earlier
	clcalculatePCI = clCreateKernel(calculatePCIprogram,"clCalculatePCI",&status);
	if(status != CL_SUCCESS)
	{
		DM::Result("Problem after compiling calculatePCI kernel");
	}

	cl_mem clQ = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<double>),0,&status);
	cl_mem clPCI = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<double>),0,&status);
	cl_mem clPCIC = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<double>),0,&status);
	cl_mem clPCIM = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<double>),0,&status);
	cl_mem clPCIS = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<double>),0,&status);


	// Set Kernel Arguments
	status  = clSetKernelArg(clgetQ,0,sizeof(cl_mem),&clW);
	status |= clSetKernelArg(clgetQ,1,sizeof(cl_mem),&clWminus);
	status |= clSetKernelArg(clgetQ,2,sizeof(cl_mem),&clV);
	status |= clSetKernelArg(clgetQ,3,sizeof(cl_mem),&clQ);
	status |= clSetKernelArg(clgetQ,4,sizeof(int),&sizeX);
	status |= clSetKernelArg(clgetQ,5,sizeof(int),&sizeY);

	// Set Kernel Arguments
	status = clSetKernelArg(clminuswavefunction,0,sizeof(cl_mem),&clRestored);
	status |= clSetKernelArg(clminuswavefunction,1,sizeof(cl_mem),&clRestoredMinus);
	status |= clSetKernelArg(clminuswavefunction,2,sizeof(int),&sizeX);
	status |= clSetKernelArg(clminuswavefunction,3,sizeof(int),&sizeY);


	status = clEnqueueNDRangeKernel(cmdQueue,clgetQ,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	status = clEnqueueNDRangeKernel(cmdQueue,clminuswavefunction,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clQ, &dataBuffer2, clMedBuffer );


	//Get sum of Q by looking at zero freq component of the FFT - too lazy to write a reduction sum.
	std::complex<double> sumQ = 0;

	clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &sumQ, 
				0, NULL, NULL );


	// Trial lots of C1 Values

	int trials = numsteps;
	
	float* Fpci = new float[trials];
	float* Fpcic = new float[trials];
	float* Fpcis = new float[trials];
	float* Fpcim = new float[trials];
	float* Fpci2 = new float[trials];

	float A1r = 0;
	float A1i = 0;
	float zero2 = 0;
	
	float StartDf = fstartrange;
	float StepDf = searchstep;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	// Set Kernel Arguments
	status  = clSetKernelArg(clcalculatePCI,0,sizeof(cl_mem),&clQ);
	status |= clSetKernelArg(clcalculatePCI,1,sizeof(cl_mem),&clRestored);
	status |= clSetKernelArg(clcalculatePCI,2,sizeof(cl_mem),&clRestoredMinus);
	status |= clSetKernelArg(clcalculatePCI,3,sizeof(cl_mem),&CLxFrequencies);
	status |= clSetKernelArg(clcalculatePCI,4,sizeof(cl_mem),&CLyFrequencies);
	status |= clSetKernelArg(clcalculatePCI,5,sizeof(cl_mem),&clPCI);
	status |= clSetKernelArg(clcalculatePCI,6,sizeof(cl_mem),&clPCIC);
	status |= clSetKernelArg(clcalculatePCI,7,sizeof(cl_mem),&clPCIM);
	status |= clSetKernelArg(clcalculatePCI,8,sizeof(cl_mem),&clPCIS);
	status |= clSetKernelArg(clcalculatePCI,9,sizeof(int),&sizeX);
	status |= clSetKernelArg(clcalculatePCI,10,sizeof(int),&sizeY);
	// Not sure if should be C3 or 0 because its corrected now?
	status |= clSetKernelArg(clcalculatePCI,12,sizeof(float),&zero2);
	status |= clSetKernelArg(clcalculatePCI,13,sizeof(float),&A1r);
	status |= clSetKernelArg(clcalculatePCI,14,sizeof(float),&A1i);
	status |= clSetKernelArg(clcalculatePCI,15,sizeof(float),&wavelength);
	status |= clSetKernelArg(clcalculatePCI,16,sizeof(float),&objap);
	
	


	for(int trial = 0 ; trial < trials ; trial++)
	{
		std::complex<double> fpcitrial = 0;
		std::complex<double> fpcictrial = 0;
		std::complex<double> fpcistrial = 0;
		std::complex<double> fpcimtrial = 0;
 
		float C1 = StartDf + trial*StepDf;

		status |= clSetKernelArg(clcalculatePCI,11,sizeof(float),&C1);

		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.
		status = clEnqueueNDRangeKernel(cmdQueue,clcalculatePCI,2,NULL,globalWorkSize, NULL,0,NULL,NULL);

		// Get sums of PCI(C/S/M) from FFT
		// TODO: proper reduction sum not FFT...
		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCI, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcitrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIC, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcictrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIS, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcistrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIM, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcimtrial, 
				0, NULL, NULL );


		Fpci[trial] = fpcitrial.real()/sumQ.real();
		Fpcic[trial] = fpcictrial.real()/sumQ.real();
		Fpcis[trial] = fpcistrial.real()/sumQ.real();
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcim[trial] = fpcimtrial.real()/sumQ.real();

	}



	/*

	// Trial lots of C1 Values

	int trials = numsteps;
	
	float* Fpci = new float[trials];
	float* Fpcic = new float[trials];
	float* Fpcis = new float[trials];
	float* Fpcim = new float[trials];
	float* Fpci2 = new float[trials];

	float A1r = 0;
	float A1i = 0;
	
	float StartDf = fstartrange;
	float StepDf = searchstep;

	DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	pcidisplay.GetOrCreateImageDocument().Show();
	pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);
	
	


	for(int trial = 0 ; trial < trials ; trial++)
	{
		double fpcitrial = 0;
		double fpcictrial = 0;
		double fpcistrial = 0;
		double fpcimtrial = 0;
 
		float C1 = StartDf + trial*StepDf;

		Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		float* pciData2 = (float*) pciLocker2.get();

		// Declare vars for OpenMP loop...
		int iLoop;
		int indexPos;
		int indexPos2;
		double arg1;
		double arg2;
		float kx1;
		float ky1;
		float freq1;
		float angle1;
		double gamma1;

		#pragma omp parallel for private(iLoop,indexPos,indexPos2,arg1,arg2,kx1,ky1,angle1,gamma1,freq1)
		for(int j = 0 ; j < sizeY ; j++)
			for( iLoop = 0 ; iLoop < sizeX ; iLoop++)
			{
				int indexPos = iLoop+j*sizeX;
				int indexPos2 = sizeX-iLoop + (sizeY-j)*sizeX;

				if(iLoop==0)
					indexPos2-=sizeX;

				if(j==0)
					indexPos2-=sizeX*sizeY;


				arg1 = atan2(fft2[indexPos][1],fft2[indexPos][0]);
				arg2 = atan2(fft2[indexPos2][1],fft2[indexPos2][0]);
				
				kx1 = xFrequencies2[iLoop];
				ky1 = yFrequencies2[j];

				freq1 = sqrt(kx1*kx1 + ky1*ky1);

				angle1 = atan2(ky1,kx1);

				gamma1 = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1r * kx1 * kx1 + 6 * C1 * kx1 * kx1 - 6 * A1r * ky1 * ky1 
				  + 6 * C1 * ky1 * ky1 + 3 * C3 * wavelength * wavelength * kx1 * kx1 * kx1 * kx1 + 6 * C3 * wavelength * wavelength * kx1 * kx1 * ky1 * ky1
				  + 3 * C3 * wavelength * wavelength * ky1 * ky1 * ky1 * ky1 + 12 * A1i * kx1 * ky1);

				#pragma omp atomic
				fpcitrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1);
				#pragma omp atomic
				fpcictrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)*cos(2*angle1);
				#pragma omp atomic
				fpcistrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)*sin(2*angle1);
				#pragma omp atomic
				fpcimtrial += (freq1<=objap)*(fpcitrial + sqrt(fpcictrial*fpcictrial + fpcistrial*fpcistrial));

				pciData2[iLoop+j*sizeX]= (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)/sumQ;
			}
			
			
			pciLocker2.Notify_DataChanged();
			pciLocker2.~ImageDataLocker();

			// Works for a while then stops.
			//DM::UpdateDisplay(pcidisplay,1,sizeY);
			DM::UpdateImage(pcidisplay);
			

			Fpci[trial] = fpcitrial/sumQ;
			Fpcic[trial] = fpcictrial/sumQ;
			Fpcis[trial] = fpcistrial/sumQ;
			// Not sure if this is ok, should do sumQ inside loop.
			Fpcim[trial] = fpcimtrial/sumQ;
	}

	*/
	

	
			

	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < trials ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = StartDf + trial*StepDf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < trials ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial <trials ; trial++)
	{
		if(Fpci2[trial] > maxheight2)
		{
			maxheight2 = Fpci2[trial];
			C1b = StartDf + trial*StepDf;
		}
	}

	float DefocusGuess = (C1a + C1b) /2;
	float A1rGuess = 1.5*cos(2*astigaxis)*(C1a-C1b)/2;
	float A1iGuess = 1.5*sin(2*astigaxis)*(C1a-C1b)/2;

	CString dg = t_to_string(DefocusGuess).c_str();
	CString A1rg = t_to_string(A1rGuess).c_str();
	CString A1ig = t_to_string(A1iGuess).c_str();
	// Print out some results and graphs

	DM::OpenAndSetProgressWindow("defocus - "+dg,"A1(real) - "+A1rg,"A1(imag) - "+A1ig);

	DM::Image pcigraph = DM::RealImage("PCI",4,trials);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial <trials ; trial++)
	{
		pciData[trial]=Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,fstartrange,searchstep,"nm",0);
	DM::UpdateImage(pcigraph);



	// Redo more accurate PCI in correct area 
	// Trial lots of C1 Values

	int trials2 = 50;
	
	float* Fpciref = new float[trials2];
	float* Fpcicref = new float[trials2];
	float* Fpcisref = new float[trials2];
	float* Fpcimref = new float[trials2];
	float* Fpci2ref = new float[trials2];

	A1r = A1rGuess;
	A1i = A1iGuess;
	
	StartDf = DefocusGuess-10;
	StepDf = 21.0f/trials2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	// Set Kernel Arguments - DONE ALREADY NOT CHANGED
	//status  = clSetKernelArg(clcalculatePCI,0,sizeof(cl_mem),&clQ);
	//status |= clSetKernelArg(clcalculatePCI,1,sizeof(cl_mem),&clRestored);
	//status |= clSetKernelArg(clcalculatePCI,2,sizeof(cl_mem),&clRestoredMinus);
	//status |= clSetKernelArg(clcalculatePCI,3,sizeof(cl_mem),&CLxFrequencies);
	//status |= clSetKernelArg(clcalculatePCI,4,sizeof(cl_mem),&CLyFrequencies);
	//status |= clSetKernelArg(clcalculatePCI,5,sizeof(cl_mem),&clPCI);
	//status |= clSetKernelArg(clcalculatePCI,6,sizeof(cl_mem),&clPCIC);
	//status |= clSetKernelArg(clcalculatePCI,7,sizeof(cl_mem),&clPCIM);
	//status |= clSetKernelArg(clcalculatePCI,8,sizeof(cl_mem),&clPCIS);
	//status |= clSetKernelArg(clcalculatePCI,9,sizeof(int),&sizeX);
	//status |= clSetKernelArg(clcalculatePCI,10,sizeof(int),&sizeY);
	// Not sure if should be C3 or 0 because its corrected now?
	//status |= clSetKernelArg(clcalculatePCI,12,sizeof(float),&C3);
	status |= clSetKernelArg(clcalculatePCI,13,sizeof(float),&A1r);
	status |= clSetKernelArg(clcalculatePCI,14,sizeof(float),&A1i);
	//status |= clSetKernelArg(clcalculatePCI,15,sizeof(float),&wavelength);
	//status |= clSetKernelArg(clcalculatePCI,16,sizeof(float),&objap);
	
	


	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		std::complex<double> fpcitrial = 0;
		std::complex<double> fpcictrial = 0;
		std::complex<double> fpcistrial = 0;
		std::complex<double> fpcimtrial = 0;
 
		float C1 = StartDf + trial*StepDf;

		status |= clSetKernelArg(clcalculatePCI,11,sizeof(float),&C1);

		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.
		status = clEnqueueNDRangeKernel(cmdQueue,clcalculatePCI,2,NULL,globalWorkSize, NULL,0,NULL,NULL);

		// Get sums of PCI(C/S/M) from FFT
		// TODO: proper reduction sum not FFT...
		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCI, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcitrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIC, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcictrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIS, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcistrial, 
				0, NULL, NULL );

		status = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clPCIM, &dataBuffer2, clMedBuffer );

		clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeof(std::complex<double>) , &fpcimtrial, 
				0, NULL, NULL );


		Fpciref[trial] = fpcitrial.real()/sumQ.real();
		Fpcicref[trial] = fpcictrial.real()/sumQ.real();
		Fpcisref[trial] = fpcistrial.real()/sumQ.real();
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcimref[trial] = fpcimtrial.real()/sumQ.real();

	}


	/*
	// Now do a more accurate PCI 2 more times to refine astigmatism estimate...

	// Trial lots of C1 Values

	int trials2 = 50;
	
	float* Fpcia = new float[50];
	float* Fpcica = new float[50];
	float* Fpcisa = new float[50];
	float* Fpcima = new float[50];
	float* Fpci2a = new float[50];

	float A1ra = A1rGuess;
	float A1ia = A1iGuess;
	float C3a = 0;
	float StartDfa=DefocusGuess-25;
	float StepDfa = 1;

	DM::Image pcidisplaya = DM::RealImage("PCI Refined",4,sizeX,sizeY);
	pcidisplaya.GetOrCreateImageDocument().Show();
	pcidisplaya.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	pcidisplaya.GetImageDisplay(0).SetSurveyTechnique(1);
	
	


	for(int trial = 0 ; trial < 50 ; trial++)
	{
		double fpcitrial = 0;
		double fpcictrial = 0;
		double fpcistrial = 0;
		double fpcimtrial = 0;
 
		float C1 = StartDfa + trial*StepDfa;

		Gatan::PlugIn::ImageDataLocker pciLocker2a(pcidisplaya);
		float* pciData2 = (float*) pciLocker2a.get();

		// Declare vars for OpenMP loop...
		int iLoop;
		int indexPos;
		int indexPos2;
		double arg1;
		double arg2;
		float kx1;
		float ky1;
		float freq1;
		float angle1;
		double gamma1;

		#pragma omp parallel for private(iLoop,indexPos,indexPos2,arg1,arg2,kx1,ky1,angle1,gamma1,freq1)
		for(int j = 0 ; j < sizeY ; j++)
			for( iLoop = 0 ; iLoop < sizeX ; iLoop++)
			{
				int indexPos = iLoop+j*sizeX;
				int indexPos2 = sizeX-iLoop + (sizeY-j)*sizeX;

				if(iLoop==0)
					indexPos2-=sizeX;

				if(j==0)
					indexPos2-=sizeX*sizeY;


				arg1 = atan2(fft2[indexPos][1],fft2[indexPos][0]);
				arg2 = atan2(fft2[indexPos2][1],fft2[indexPos2][0]);
				
				kx1 = xFrequencies2[iLoop];
				ky1 = yFrequencies2[j];

				freq1 = sqrt(kx1*kx1 + ky1*ky1);

				angle1 = atan2(ky1,kx1);

				gamma1 = ( .166666666 * 3.1415927 * wavelength ) * ( 6 * A1ra * kx1 * kx1 + 6 * C1 * kx1 * kx1 - 6 * A1ra * ky1 * ky1 
				  + 6 * C1 * ky1 * ky1 + 3 * C3 * wavelength * wavelength * kx1 * kx1 * kx1 * kx1 + 6 * C3 * wavelength * wavelength * kx1 * kx1 * ky1 * ky1
				  + 3 * C3 * wavelength * wavelength * ky1 * ky1 * ky1 * ky1 + 12 * A1ia * kx1 * ky1);

				#pragma omp atomic
				fpcitrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1);
				#pragma omp atomic
				fpcictrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)*cos(2*angle1);
				#pragma omp atomic
				fpcistrial += (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)*sin(2*angle1);
				#pragma omp atomic
				fpcimtrial += (freq1<=objap)*(fpcitrial + sqrt(fpcictrial*fpcictrial + fpcistrial*fpcistrial));

				pciData2[iLoop+j*sizeX]= (freq1<=objap)*Q[iLoop+j*sizeX]* -1 *cos(arg1+arg2+2*gamma1)/sumQ;
			}
			
			
			pciLocker2a.Notify_DataChanged();
			pciLocker2a.~ImageDataLocker();

			// Works for a while then stops.
			//DM::UpdateDisplay(pcidisplay,1,sizeY);
			DM::UpdateImage(pcidisplaya);
			

			Fpcia[trial] = fpcitrial/sumQ;
			Fpcica[trial] = fpcictrial/sumQ;
			Fpcisa[trial] = fpcistrial/sumQ;
			// Not sure if this is ok, should do sumQ inside loop.
			Fpcima[trial] = fpcimtrial/sumQ;
	}

	
	*/
	
	// Find max

	float maxheighta = 0;
	float C1aa = 0;
	int besttriala = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpcimref[trial] > maxheighta)
		{
			maxheighta = Fpcimref[trial];
			C1aa = StartDf + trial*StepDf;
			besttriala = trial;
		}
	}

	float astigaxisa = atan2(Fpcisref[besttriala],Fpcicref[besttriala])/2;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		Fpci2ref[trial] = Fpciref[trial] - Fpcicref[trial] * cos(2*astigaxisa) - Fpcisref[trial]*sin(2*astigaxisa);
	}

	float maxheight2a = 0;
	float C1ba = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpci2ref[trial] > maxheight2a)
		{
			maxheight2a = Fpci2ref[trial];
			C1ba = StartDf + trial*StepDf;
		}
	}

	float DefocusGuess2 = (C1aa + C1ba) /2;
	float A1rGuess2 = 1.5*cos(2*astigaxisa)*(C1aa-C1ba)/2 + A1rGuess;
	float A1iGuess2 = 1.5*sin(2*astigaxisa)*(C1aa-C1ba)/2 + A1iGuess;

	CString dga = t_to_string(DefocusGuess2).c_str();
	CString A1rga = t_to_string(A1rGuess2).c_str();
	CString A1iga = t_to_string(A1iGuess2).c_str();
	// Print out some results and graphs

	DM::OpenAndSetProgressWindow("defocus - "+dga,"A1(real) - "+A1rga,"A1(imag) - "+A1iga);

	DM::Image pcigraph2 = DM::RealImage("PCI - Refined",4,trials2);
	Gatan::PlugIn::ImageDataLocker pciLocker2(pcigraph2);
	float* pciData2 = (float*) pciLocker2.get();

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		pciData2[trial]=Fpciref[trial];
	}

	pciLocker2.~ImageDataLocker();
	pcigraph2.GetOrCreateImageDocument().Show();
	pcigraph2.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);


	// Then do restoration again at correct focal value

	Gatan::PlugIn::ImageDataLocker driftCorrectedLocker3(driftCorrected);
	float* driftCorrecteddata3 = (float*) driftCorrectedLocker3.get();

	// Process each Image
	for(int image = 0 ; image < numberOfImages ; image++)
	{
		// Images need to be normalised before FFT for FSR
		float sum = GetSum(driftCorrecteddata3,sizeX*sizeY,image*sizeX*sizeY);
		
		for(int j = 0 ; j < sizeY ; j++)
			for(int i = 0 ; i < sizeX ; i++)
			{
				dataOne[i+j*sizeX] = (sizeX*sizeY*driftCorrecteddata3[image*sizeX*sizeY + j*sizeX + i]/(sum)) - 1;		
			}

		clEnqueueWriteBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &dataOne[ 0 ], 
				0, NULL, NULL );

		fftStatus = clAmdFftEnqueueTransform( fftplan, CLFFT_FORWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&dataBuffer2, &dataBuffer1, clMedBuffer );
		
	
		// Now calculate w,wminus and add to W,W-,T,U,V
		// Because of image order first is most underfocus

		// Changed to try do first restore in plane of last image (usually around focus :))
		float defocus = image*fStep - numberOfImages*fStep + DefocusGuess2;
		float zero = 0;
		int init = 0;

		if(image==0)
		{
			init = 1;

			if(gotmtf)
			{
			status |= clSetKernelArg(clwtfmtf,8,sizeof(float),&A1rGuess2);
			status |= clSetKernelArg(clwtfmtf,9,sizeof(float),&A1iGuess2);
			
			// Set Kernel Arguments
			status |= clSetKernelArg(clwtfminusmtf,8,sizeof(float),&A1rGuess2);
			status |= clSetKernelArg(clwtfminusmtf,9,sizeof(float),&A1iGuess2);
			}
			else
			{
			status |= clSetKernelArg(clwtf,8,sizeof(float),&A1rGuess2);
			status |= clSetKernelArg(clwtf,9,sizeof(float),&A1iGuess2);
			
			// Set Kernel Arguments
			status |= clSetKernelArg(clwtfminus,8,sizeof(float),&A1rGuess2);
			status |= clSetKernelArg(clwtfminus,9,sizeof(float),&A1iGuess2);
			}
			
		}
		
		// Only these two change
		if(gotmtf)
		{
			// Only these two change
			status |= clSetKernelArg(clwtfmtf,10,sizeof(float),&defocus);
			status |= clSetKernelArg(clwtfminusmtf,10,sizeof(float),&defocus);
		}
		else
		{
			// Only these two change
			status |= clSetKernelArg(clwtf,10,sizeof(float),&defocus);
			status |= clSetKernelArg(clwtfminus,10,sizeof(float),&defocus);
		}

		// Oops, And this one!!
		status |= clSetKernelArg(clwienerw,4,sizeof(int),&init);
		status |= clSetKernelArg(clwienerwminus,4,sizeof(int),&init);
		status |= clSetKernelArg(clwienerv,5,sizeof(int),&init);
		status |= clSetKernelArg(clwienert,5,sizeof(int),&init);
		status |= clSetKernelArg(clwieneru,5,sizeof(int),&init);


		// Define and index space of work items for execution
		// A workgroup size is not required but can be used
		size_t* globalWorkSize = new size_t[2];
		
		// There are 'elements' work items
		globalWorkSize[0] = sizeX;
		globalWorkSize[1] = sizeY;


		// Messed up when I tried this
		size_t* localWorkSize = new size_t[2];
		localWorkSize[0] = 32;
		localWorkSize[1] = 16;

		if(gotmtf)
		{
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfmtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfminusmtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		}
		else
		{
			status = clEnqueueNDRangeKernel(cmdQueue,clwtf,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
			status = clEnqueueNDRangeKernel(cmdQueue,clwtfminus,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		}
		
		status = clEnqueueNDRangeKernel(cmdQueue,clwienerw,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienerwminus,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienerv,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwienert,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
		status = clEnqueueNDRangeKernel(cmdQueue,clwieneru,2,NULL,globalWorkSize,NULL,0,NULL,NULL);
	}

	clFinish(cmdQueue);
	driftCorrectedLocker3.~ImageDataLocker();
	/* Did this earlier
	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[2];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	*/

	if(gotmtf&&gotnps)
	{
		status = clEnqueueNDRangeKernel(cmdQueue,clmakerestoredMTFNPS,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	}
	else
	{
		status = clEnqueueNDRangeKernel(cmdQueue,clmakerestored,2,NULL,globalWorkSize, NULL,0,NULL,NULL);
	}
	status = clAmdFftEnqueueTransform( fftplan, CLFFT_BACKWARD, 1, &cmdQueue, 0, NULL, NULL, 
			&clRestored, &dataBuffer2, clMedBuffer ); 

	clFinish(cmdQueue);

	// Actually want FFT of this :)

	std::vector< std::complex< double > > HostRestored2( sizeX*sizeY );
	clEnqueueReadBuffer( cmdQueue, dataBuffer2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<double>) , &HostRestored2[ 0 ], 
				0, NULL, NULL );


	/*
	if(gotmtf&&gotnps)
	{
		MakeRestoredWavefunctionmtfnps(W,Wminus,V,T,U,fft2,sizeX,sizeY,mtfdata,npsdata,mtflength,npslength,scaleMTFx,scaleMTFy);
	}
	else
	{
		MakeRestoredWavefunction(W,Wminus,V,T,U,fft2,sizeX,sizeY);
	}
	
	fftw_execute(fftplanInv2);
	*/
	
	DM::Image RestoredCorrected = DM::ComplexImage("RestoredWave - Corrected",16,sizeX,sizeY);

	Gatan::PlugIn::ImageDataLocker RestoredLocker2(RestoredCorrected);
	complex128* RestoredData2 = (complex128*) RestoredLocker2.get();

	for(int j = 0 ; j < sizeY ; j++)
		for(int i = 0 ; i < sizeX ; i++)
			{
				RestoredData2[i+j*sizeX].real() = HostRestored2[i+j*sizeX].real()+1;
				RestoredData2[i+j*sizeX].imag() = HostRestored2[i+j*sizeX].imag();
			}

	
	RestoredLocker2.~ImageDataLocker();
	RestoredCorrected.SetDimensionCalibration(0,0,pixelscale,"nm",0);
	RestoredCorrected.SetDimensionCalibration(1,0,pixelscale,"nm",0);
	RestoredCorrected.GetOrCreateImageDocument().Show();
	RestoredCorrected.GetImageDisplay(0).SetComplexMode(5);

	// More Cleanup
	HostRestored2.clear();


	/*
	// Process each Image
	for(int image = 0 ; image < numberOfImages ; image++)
	{
		// Images need to be normalised before FFT for FSR
		float sum = GetSum(driftCorrecteddata,sizeX*sizeY,image*sizeX*sizeY);
		
		for(int j = 0 ; j < sizeY ; j++)
			for(int i = 0 ; i < sizeX ; i++)
			{
				fft1[i+j*sizeX][0] = (sizeX*sizeY*driftCorrecteddata[image*sizeX*sizeY + j*sizeX + i]/(sum)) - 1;
				fft1[i+j*sizeX][1] = 0; 
			}

		fftw_execute(fftplan2);

		// Also normalise after FFT to rescale properly 1/sqrt(N)

		for(int j = 0 ; j < sizeY ; j++)
			for(int i = 0 ; i < sizeX ; i++)
			{
				fft2[i+j*sizeX][0] /= normalise;
				fft2[i+j*sizeX][1] /= normalise; 
			}

		// Now calculate w,wminus and add to W,W-,T,U,V

		// Because of image order first is most underfocus


		if(gotmtf&&gotnps)
		{
			wavetransferfunctionmtf(w,xFrequencies2,yFrequencies2,sizeX,sizeY,wavelength,beta,delta,A1rGuess2,A1iGuess2,0,0,DefocusGuess2+image*fStep,C3,objap,mtfdata,mtflength,scaleMTFx,scaleMTFy);
			wavetransferfunctionminusmtf(wminus,xFrequencies2,yFrequencies2,sizeX,sizeY,wavelength,beta,delta,A1rGuess2,A1iGuess2,0,0,DefocusGuess2+image*fStep,C3,objap,mtfdata,mtflength,scaleMTFx,scaleMTFy);
		}
		else
		{
			wavetransferfunction(w,xFrequencies2,yFrequencies2,sizeX,sizeY,wavelength,beta,delta,A1rGuess2,A1iGuess2,0,0,DefocusGuess2+image*fStep,C3,objap);
			wavetransferfunctionminus(wminus,xFrequencies2,yFrequencies2,sizeX,sizeY,wavelength,beta,delta,A1rGuess2,A1iGuess2,0,0,DefocusGuess2+image*fStep,C3,objap);
		}

		WienerW(W,w,sizeX,sizeY,image==0);
		WienerW(Wminus,wminus,sizeX,sizeY,image==0);
		WienerT(T,fft2,w,sizeX,sizeY,image==0);
		WienerU(U,fft2,wminus,sizeX,sizeY,image==0);
		WienerV(V,w,wminus,sizeX,sizeY,image==0);
	}

	
	if(gotmtf&&gotnps)
	{
		MakeRestoredWavefunctionmtfnps(W,Wminus,V,T,U,fft2,sizeX,sizeY,mtfdata,npsdata,mtflength,npslength,scaleMTFx,scaleMTFy);
	}
	else
	{
		MakeRestoredWavefunction(W,Wminus,V,T,U,fft2,sizeX,sizeY);
	}
	
	fftw_execute(fftplanInv2);

	// Now need to normalise then add 1 and stick in an image somewhere...
	//for(int j = 0 ; j < sizeY ; j++)
	//	for(int i = 0 ; i < sizeX ; i++)
	//		{
	//			fft1[i+j*sizeX][0] = fft1[i+j*sizeX][0] / normalise + 1;
	//			fft1[i+j*sizeX][1] = fft1[i+j*sizeX][1] / normalise + 1; 
	//		}
	
	DM::Image Restored2 = DM::ComplexImage("RestoredWave - Corrected",8,sizeX,sizeY);

	Gatan::PlugIn::ImageDataLocker RestoredLocker2(Restored2);
	complex64* RestoredData2 = (complex64*) RestoredLocker2.get();

	for(int j = 0 ; j < sizeY ; j++)
		for(int i = 0 ; i < sizeX ; i++)
			{
				RestoredData2[i+j*sizeX].real() = fft1[i+j*sizeX][0] / normalise + 1;
				RestoredData2[i+j*sizeX].imag() = fft1[i+j*sizeX][1] / normalise;
			}

	
	RestoredLocker2.~ImageDataLocker();
	Restored2.SetDimensionCalibration(0,0,pixelscale,"nm",0);
	Restored2.SetDimensionCalibration(1,0,pixelscale,"nm",0);
	Restored2.GetOrCreateImageDocument().Show();
	Restored2.GetImageDisplay(0).SetComplexMode(5);


	fftw_complex* fft1a = new fftw_complex[sizeX*sizeY];
	fftw_complex* fft2a = new fftw_complex[sizeX*sizeY];

	
	fftw_plan fftplanInv3;

	// These Plans can only operate on the two specified arrays, have to copy data to these locations to transform
	fftplanInv3 = fftw_plan_dft_2d(sizeY, sizeX, fft2a, fft1a, FFTW_BACKWARD, FFTW_ESTIMATE);

	// Maybe simulate the actual focal series again aswell

	DM::Image simulatedSeries = DM::RealImage("Simulated Stack",4,sizeX,sizeY,numberOfImages);

	Gatan::PlugIn::ImageDataLocker simulatedLocker(simulatedSeries);

	float* simulateddata = (float*) simulatedLocker.get();

	for(int image = 0; image < numberOfImages; image++ )
	{
		float focus = DefocusGuess2 + image*fStep;	
		makecorrected(xFrequencies2,yFrequencies2,sizeX,sizeY,wavelength,-A1rGuess2,-A1iGuess2,-focus,-C3,fft2,fft2a);
		fftw_execute(fftplanInv3);

		for(int j = 0; j < sizeY; j++)
			for(int i = 0; i < sizeX; i++)
			{
				simulateddata[image*sizeX*sizeY + j*sizeX + i] = sqrt(((fft1a[j*sizeX + i][0]/normalise)+1)* ((fft1a[j*sizeX + i][0]/normalise)+1) +  (fft1a[j*sizeX + i][1]/normalise) * (fft1a[j*sizeX + i][1]/normalise));
			}
	}

	simulatedLocker.~ImageDataLocker();
	simulatedSeries.SetDimensionCalibration(0,0,pixelscale,"nm",0);
	simulatedSeries.SetDimensionCalibration(1,0,pixelscale,"nm",0);
	simulatedSeries.GetOrCreateImageDocument().Show();


	// Final Cleanup

	fftw_destroy_plan(fftplan2);
	fftw_destroy_plan(fftplanInv2);
	fftw_destroy_plan(fftplanInv3);

	delete[] fft1a;
	delete[] fft2a;

	delete[] Fpci;
	delete[] Fpcic;
	delete[] Fpcis;
	delete[] Fpcim;
	delete[] Fpci2;
	delete[] xFrequencies2;
	delete[] yFrequencies2;
	delete[] W;
	delete[] Wminus;
	delete[] T;
	delete[] U;
	delete[] V;
	delete[] w;
	delete[] wminus;
	delete[] fft1;
	delete[] fft2;
	*/

// Now time to clean up everything...

	dataOne.clear();
	dataTwo.clear();

	clReleaseKernel(clPCPCF);
	clReleaseProgram(program);
	clReleaseProgram(wavetransferfunctionprogram);
	clReleaseKernel(clwtf);
	clReleaseProgram(wavetransferfunctionminusprogram);
	clReleaseKernel(clwtfminus);
	clReleaseProgram(wavetransferfunctionmtfprogram);
	clReleaseKernel(clwtfmtf);
	clReleaseProgram(wavetransferfunctionminusmtfprogram);
	clReleaseKernel(clwtfminusmtf);
	clReleaseProgram(wienerwprogram);
	clReleaseKernel(clwienerw);
	clReleaseProgram(wienerwminusprogram);
	clReleaseKernel(clwienerwminus);
	clReleaseProgram(wienervprogram);
	clReleaseKernel(clwienerv);
	clReleaseProgram(wienertprogram);
	clReleaseKernel(clwienert);
	clReleaseProgram(wieneruprogram);
	clReleaseKernel(clwieneru);
	clReleaseProgram(makerestoredprogram);
	clReleaseKernel(clmakerestored);
	clReleaseProgram(makerestoredMTFNPSprogram);
	clReleaseKernel(clmakerestoredMTFNPS);

	clReleaseProgram(getQprogram);
	clReleaseKernel(clgetQ);
	clReleaseProgram(minuswavefunctionprogram);
	clReleaseKernel(clminuswavefunction);
	clReleaseProgram(calculatePCIprogram);
	clReleaseKernel(clcalculatePCI);

	clReleaseMemObject(clW);
	clReleaseMemObject(clWminus);
	clReleaseMemObject(clw);
	clReleaseMemObject(clwminus);
	clReleaseMemObject(clT);
	clReleaseMemObject(clV);
	clReleaseMemObject(clU);
	clReleaseMemObject(clRestored);
	clReleaseMemObject(clRestoredMinus);
	clReleaseMemObject(clQ);
	clReleaseMemObject(clPCI);
	clReleaseMemObject(clPCIC);
	clReleaseMemObject(clPCIS);
	clReleaseMemObject(clPCIM);
	clReleaseMemObject(CLxFrequencies);
	clReleaseMemObject(CLyFrequencies);
	clReleaseMemObject(dataBuffer1);
	clReleaseMemObject(dataBuffer2);
	clReleaseMemObject(clMedBuffer);

	// Maybe I should just be destroying plan?
	clAmdFftTeardown();
}
