#include "stdafx.h"
#include "PCPCFFunction.h"
#include "vector"
#include "Utility.h"
#include "clKernels2.h"
#include "PCPCFLib.h"



void PCPCFWrapper(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice *cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb)
{
	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	cl_int status;
	// Construct PCPCF Kernel

	clKernel* PCPCF = new clKernel(pcpcfSource,context,cldev,"clPCPCF",clq);
	PCPCF->BuildKernel();

	// Setup required memory for PCPCF subroutine

	cl_mem clImage1 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clImage2 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage1 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage2 = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clPCPCFResult = clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);



	PCPCF->SetArgT(0,clFFTImage1);
	PCPCF->SetArgT(1,clFFTImage2);
	PCPCF->SetArgT(2,clPCPCFResult);
	PCPCF->SetArgT(3,clxFrequencies);
	PCPCF->SetArgT(4,clyFrequencies);
	PCPCF->SetArgT(5,sizeX);
	PCPCF->SetArgT(6,sizeY);
	PCPCF->SetArgT(8,wavelength);

	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[3];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	globalWorkSize[2] = 1;

	std::vector< std::complex< float > > dataOne( sizeX*sizeY );
	std::vector< std::complex< float > > dataTwo( sizeX*sizeY );

	// loop over all of the images
	for(int n = 0; n < numberOfImages - 1; n++)
	{
	
		//TODO: Could be clever and keep swapping the vectors to only load one image each time
		for(int j = 0; j< sizeY;j++)
			for(int i = 0; i< sizeX;i++)
			{
				dataOne[i+j*sizeX] = seriesdata[n*xDim*yDim + i + iLeft + (j+iTop)*xDim];
				dataTwo[i+j*sizeX] = seriesdata[(n+1)*xDim*yDim + i + iLeft + (j+iTop)*xDim];
			}

	
		clEnqueueWriteBuffer( clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
					0, NULL, NULL );
		clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
					0, NULL, NULL );

		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
		FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);


		// Work out what focus difference is expected between these 2 images.

		float expecteddifference = options.focalstep;

		int numberoftrials = 1;

		if(options.determinefocus)
		{
			numberoftrials = options.steps;
		}



		int* trialxshifts = new int[numberoftrials];
		int* trialyshifts = new int[numberoftrials];
		float* trialsubxshifts = new float[numberoftrials];
		float* trialsubyshifts = new float[numberoftrials];
		float* trialdefocus = new float[numberoftrials];
		float* peakheights = new float[numberoftrials];

		for(int trial = 0; trial < numberoftrials; trial++)
		{
			// go from expected - searchpercentage to expected = searchpercentage
			float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*expecteddifference 
						+ ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*expecteddifference; 	

			PCPCF->SetArgT(7,trialdef);
			PCPCF->Enqueue(globalWorkSize);

			// Now Inverse FFT
			// Using memory reserved for storing images as it is not needed at this point.
			FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);
			clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );

			// Could consider having this displayed to see roughly how well it is working...
			//Utility::PrintCLMemToImage(clFFTImage1,"FFT1",sizeX,sizeY,clFloat2,clq);
			//Utility::PrintCLMemToImage(clImage1,"PCPCF",sizeX,sizeY,clFloat2,clq);

			// Find shift from max height position
			int xShift;
			int yShift;
			float subXShift;
			float subYShift;
			float maxHeight1;

			// Translate linear array index into row and column.
			PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY);

			trialxshifts[trial] = xShift;
			trialyshifts[trial] = yShift;
			trialsubxshifts[trial] = xShift + subXShift;
			trialsubyshifts[trial] = yShift + subYShift;
			trialdefocus[trial] = trialdef;
			peakheights[trial] = maxHeight1;


		}

		// Remember the reference image has shifts of 0,0
		// Shifts now expressed in terms of difference to reference image.

		// Get best value and store in shifts.
		float bestheight = 0.0f;
		for(int trial = 0; trial < numberoftrials; trial++)
		{
			if(peakheights[trial] > bestheight)
			{
				bestheight = peakheights[trial];
				xShiftVals[n+1] = trialxshifts[trial];
				yShiftVals[n+1] = trialyshifts[trial];
				subXShifts[n+1] = trialsubxshifts[trial];
				subYShifts[n+1] = trialsubyshifts[trial];
				defocusshifts[n+1] = trialdefocus[trial];
			}
		}
	}

	dataOne.clear();
	dataTwo.clear();
	
};

// Sometimes getting funny results after a certain point in code, not sure which step is failing but it starts giving same defocus for every image which should be impossible...
void PCPCFWrapper2(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb)
{
	//Utility::SetResultWindow("Should Align = "+boost::lexical_cast<std::string>(options.determinefocus)+"\n");
	//Utility::SetResultWindow("Steps = "+boost::lexical_cast<std::string>(options.steps)+"\n");

	std::vector<int> imagelist;

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	cl_int status;
	// Construct PCPCF Kernel

	clKernel* PCPCF = new clKernel(pcpcfSource,context,cldev,"clPCPCF",clq);
	PCPCF->BuildKernel();

	// Setup required memory for PCPCF subroutine

	cl_mem clImage1			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clImage2			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage1		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage2		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clPCPCFResult	= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);

	cl_mem fullImage		= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem rotScaleImage	= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);

	cl_mem clW				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clWminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clw				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clwminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clT				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clU				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clV				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clRestored		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);

	int init = 1;
	int noinit = 0;

	PCPCF->SetArgT(0,clFFTImage1);
	PCPCF->SetArgT(1,clFFTImage2);
	PCPCF->SetArgT(2,clPCPCFResult);
	PCPCF->SetArgT(3,clxFrequencies);
	PCPCF->SetArgT(4,clyFrequencies);
	PCPCF->SetArgT(5,sizeX);
	PCPCF->SetArgT(6,sizeY);
	PCPCF->SetArgT(8,wavelength);
	PCPCF->SetArgT(9,options.pcpcfkmax);

	/* Make all the other kernels */
	clKernel* clPadCrop = new clKernel(PadCropSource,context,cldev,"clPadCrop",clq);
	clPadCrop->BuildKernel();

	clPadCrop->SetArgT(0,fullImage);
	clPadCrop->SetArgT(1,clImage1);
	clPadCrop->SetArgT(2,xDim); // long not int
	clPadCrop->SetArgT(3,yDim); // long not int
	clPadCrop->SetArgT(10,sizeX);
	clPadCrop->SetArgT(11,sizeY);
	clPadCrop->SetArgT(12,iTop);
	clPadCrop->SetArgT(13,iLeft);

	clKernel* clWaveTransferFunction = new clKernel(wavetransferfunctionsource,context,cldev,"clWaveTransferFunction",clq);
	clWaveTransferFunction->BuildKernel();

	clWaveTransferFunction->SetArgT(0,clw);
	clWaveTransferFunction->SetArgT(1,clxFrequencies);
	clWaveTransferFunction->SetArgT(2,clyFrequencies);
	clWaveTransferFunction->SetArgT(3,sizeX);
	clWaveTransferFunction->SetArgT(4,sizeY);
	clWaveTransferFunction->SetArgT(5,wavelength);
	clWaveTransferFunction->SetArgT(6,Abb.beta);
	clWaveTransferFunction->SetArgT(7,Abb.delta);
	clWaveTransferFunction->SetArgT(8,Abb.A1r);
	clWaveTransferFunction->SetArgT(9,Abb.A1i);
	clWaveTransferFunction->SetArgT(11,Abb.Cs);
	clWaveTransferFunction->SetArgT(12,Abb.kmax);

	clKernel* clWaveTransferFunctionMinus = new clKernel(wavetransferfunctionminussource,context,cldev,"clWaveTransferFunctionMinus",clq);
	clWaveTransferFunctionMinus->BuildKernel();

	clWaveTransferFunctionMinus->SetArgT(0,clwminus);
	clWaveTransferFunctionMinus->SetArgT(1,clxFrequencies);
	clWaveTransferFunctionMinus->SetArgT(2,clyFrequencies);
	clWaveTransferFunctionMinus->SetArgT(3,sizeX);
	clWaveTransferFunctionMinus->SetArgT(4,sizeY);
	clWaveTransferFunctionMinus->SetArgT(5,wavelength);
	clWaveTransferFunctionMinus->SetArgT(6,Abb.beta);
	clWaveTransferFunctionMinus->SetArgT(7,Abb.delta);
	clWaveTransferFunctionMinus->SetArgT(8,Abb.A1r);
	clWaveTransferFunctionMinus->SetArgT(9,Abb.A1i);
	clWaveTransferFunctionMinus->SetArgT(11,Abb.Cs);
	clWaveTransferFunctionMinus->SetArgT(12,Abb.kmax);

	clKernel* clWienerW = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerW->BuildKernel();

	clWienerW->SetArgT(0,clW);
	clWienerW->SetArgT(1,clw);
	clWienerW->SetArgT(2,sizeX);
	clWienerW->SetArgT(3,sizeY);
	clWienerW->SetArgT(4,init);

	clKernel* clWienerWMinus = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerWMinus->BuildKernel();

	clWienerWMinus->SetArgT(0,clWminus);
	clWienerWMinus->SetArgT(1,clwminus);
	clWienerWMinus->SetArgT(2,sizeX);
	clWienerWMinus->SetArgT(3,sizeY);
	clWienerWMinus->SetArgT(4,init);

	clKernel* clWienerV = new clKernel(wienervsource,context,cldev,"clWienerV",clq);
	clWienerV->BuildKernel();

	clWienerV->SetArgT(0,clV);
	clWienerV->SetArgT(1,clw);
	clWienerV->SetArgT(2,clwminus);
	clWienerV->SetArgT(3,sizeX);
	clWienerV->SetArgT(4,sizeY);
	clWienerV->SetArgT(5,init);

	clKernel* clWienerT = new clKernel(wienertsource,context,cldev,"clWienerT",clq);
	clWienerT->BuildKernel();

	clWienerT->SetArgT(0,clT);
	clWienerT->SetArgT(1,clw);
	clWienerT->SetArgT(2,clFFTImage1);
	clWienerT->SetArgT(3,sizeX);
	clWienerT->SetArgT(4,sizeY);
	clWienerT->SetArgT(5,init);

	clKernel* clWienerU = new clKernel(wienerusource,context,cldev,"clWienerU",clq);
	clWienerU->BuildKernel();

	clWienerU->SetArgT(0,clU);
	clWienerU->SetArgT(1,clwminus);
	clWienerU->SetArgT(2,clFFTImage1);
	clWienerU->SetArgT(3,sizeX);
	clWienerU->SetArgT(4,sizeY);
	clWienerU->SetArgT(5,init);
	
	clKernel* clMakeRestored = new clKernel(makerestoredsource,context,cldev,"clMakeRestored",clq);
	clMakeRestored->BuildKernel();

	clMakeRestored->SetArgT(0,clW);
	clMakeRestored->SetArgT(1,clWminus);
	clMakeRestored->SetArgT(2,clV);
	clMakeRestored->SetArgT(3,clT);
	clMakeRestored->SetArgT(4,clU);
	clMakeRestored->SetArgT(5,clFFTImage1);
	clMakeRestored->SetArgT(6,sizeX);
	clMakeRestored->SetArgT(7,sizeY);

	clKernel* clRotScale = new clKernel(RotScaleSource,context,cldev,"clRotScale",clq);
	clRotScale->BuildKernel();

	if(options.rotscale)
	{
		clRotScale->SetArgT(0,fullImage);
		clRotScale->SetArgT(1,rotScaleImage);
		clRotScale->SetArgT(2,xDim);
		clRotScale->SetArgT(3,yDim);
		clRotScale->SetArgT(4,options.magcal);
		clRotScale->SetArgT(5,options.rotcal);
	}


	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[3];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	globalWorkSize[2] = 1;

	size_t* fullWorkSize = new size_t[3];
	
	// There are 'elements' work items
	fullWorkSize[0] = xDim;
	fullWorkSize[1] = yDim;
	fullWorkSize[2] = 1;

	std::vector< std::complex< float > > dataOne( sizeX*sizeY );
	std::vector< std::complex< float > > dataTwo( sizeX*sizeY );

	int referenceimage = options.reference;
	int gonext = 1;
	int nextdir = -1;
	
	// First three have to be done via standard methods.

	// Initialise image numbers
	int imageone = referenceimage;
	int imagetwo = referenceimage;
	int currentimage;

	imagelist.push_back(referenceimage);

	std::vector<std::complex<float>> copyImage(xDim*yDim) ;

	// loop over all of the images
	for(int n = 0; n < numberOfImages - 1; n++)
	{
		
		 // For RECONPCPCF part

		// First 2 registrations cannot be done against reconstruction
		// have to be done against reference image...
		if(n < 2)
		{
			imageone = referenceimage;
			imagetwo += (gonext*nextdir);

			while(! (imagetwo < numberOfImages && imagetwo >=0)) // not an actual image 
			{
				gonext++;
				nextdir *=-1;
				imagetwo += (gonext*nextdir);
			}

			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = (imagetwo-imageone) * options.focalstep;

			// Rotate and Scale this image by expected defocus..
			if(options.rotscale)
			{
				clRotScale->SetArgT(6,expecteddifference);

				// Copy image into fullimage
				for(int j = 0; j < yDim;j++)
					for(int i = 0; i < xDim;i++)
					{
						copyImage[i+j*xDim] = seriesdata[imagetwo*xDim*yDim + i + (j)*xDim];
					}

					clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
					clRotScale->Enqueue(fullWorkSize);
					clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

					for(int j = 0; j < sizeY;j++)
						for(int i = 0; i < sizeX;i++)
						{
							dataOne[i+j*sizeX] = seriesdata[imageone*xDim*yDim + i + iLeft + (j+iTop)*xDim];
							dataTwo[i+j*sizeX] = copyImage[i + iLeft + (j+iTop)*xDim];
						}
			} 
			else
			{
				for(int j = 0; j < sizeY;j++)
					for(int i = 0; i < sizeX;i++)
					{
						dataOne[i+j*sizeX] = seriesdata[imageone*xDim*yDim + i + iLeft + (j+iTop)*xDim];
						dataTwo[i+j*sizeX] = seriesdata[imagetwo*xDim*yDim + i + iLeft + (j+iTop)*xDim];
					}
			}

	
			clEnqueueWriteBuffer( clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
						0, NULL, NULL );

			FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			for(int trial = 0; trial < numberoftrials; trial++)
			{

				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected + searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep + ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);

				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(imagetwo)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY);

				trialxshifts[trial] = xShift;
				trialyshifts[trial] = yShift;
				trialsubxshifts[trial] = xShift + subXShift;
				trialsubyshifts[trial] = yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;
			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[imagetwo] = trialxshifts[trial];
					yShiftVals[imagetwo] = trialyshifts[trial];
					subXShifts[imagetwo] = trialsubxshifts[trial];
					subYShifts[imagetwo] = trialsubyshifts[trial];
					defocusshifts[imagetwo] = trialdefocus[trial];
				}
			}

			imagelist.push_back(imagetwo);
		}

		if(n==2)
		{
			// To make sure it is initialised from end of first section
			currentimage = imagetwo;
		}
	
		if(n>=2)
		{
			// Now we can perform reconstruction and register to the reconstruction instead of standard pcpcf
			currentimage += (gonext*nextdir);

			while(! (currentimage < numberOfImages && currentimage >=0)) // not an actual image
			{
				gonext++;
				nextdir *=-1;
				currentimage += (gonext*nextdir);
			}

			// Need to reconstruct the previous images based on there focus and shifts relative to reference image...

			
			// Need to extract correct area from each image to perform reconstruction
			// 2 Possibilites - If an area of the correct size is not available then give an error...
			// or pad the area with zeroes for ALL images - need to first determine the extents of available area before we extract from any of the images.

			// Find maximum negative and positive shifts
			// All shifts have to be stored relative to the reference image.
			float maxnegx = 0;
			float maxnegy = 0;
			float maxposx = 0;
			float maxposy = 0;

			// Could just loop over this list of currently registered images instead...
			for(int i = 1; i <= imagelist.size(); i++)
			{
				if(subXShifts[imagelist[i-1]] < maxnegx)
					maxnegx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] < maxnegy)
					maxnegy = subYShifts[imagelist[i-1]];
				if(subXShifts[imagelist[i-1]] > maxposx)
					maxposx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] > maxposy)
					maxposy = subYShifts[imagelist[i-1]];
			}

			// Determine amount to pad on either direction
			int padTop = 0;
			int padLeft = 0;
			int padRight = 0;
			int padBottom = 0;

			if (abs(maxnegx)-iLeft > 0)
				padLeft = ceil(abs(maxnegx)-iLeft);

			if (maxposx-iRight > 0)
				padRight = ceil(maxposx-iRight);

			if (abs(maxnegy)-iTop > 0)
				padTop = ceil(abs(maxnegy)-iTop);

			if (maxposy-iBottom > 0)
				padBottom = ceil(maxposy-iBottom);	

			for(int i = 0 ; i < imagelist.size() ; i++)
			{
				if(i == 0) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,init);
					clWienerWMinus->SetArgT(4,init);
					clWienerV->SetArgT(5,init);
					clWienerT->SetArgT(5,init);
					clWienerU->SetArgT(5,init);
				}
				if(i == 1) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,noinit);
					clWienerWMinus->SetArgT(4,noinit);
					clWienerV->SetArgT(5,noinit);
					clWienerT->SetArgT(5,noinit);
					clWienerU->SetArgT(5,noinit);
				}
				
				// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
				clPadCrop->SetArgT(4,padLeft);
				clPadCrop->SetArgT(5,padRight);
				clPadCrop->SetArgT(6,padTop);
				clPadCrop->SetArgT(7,padBottom);
				clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
				clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
				// Copy correct image into full image..
				std::vector<cl_float2> image(xDim*yDim);
				/*
				for(int j = 0 ; j < xDim*yDim ; j++)
				{
					image[j].s[0] = seriesdata[imagelist[i]*xDim*yDim + j];
					image[j].s[1] = 0;
				}
				*/

				// Now add image to reconstruction - Make WTF - WTFminus			
				float defocus = defocusshifts[imagelist[i]];

				// Rotate and Scale this image by expected defocus..
				if(options.rotscale)
				{
					clRotScale->SetArgT(6,defocus);
					// Copy image into fullimage
					for(int j = 0 ; j < xDim*yDim ; j++)
					{
						copyImage[j] = seriesdata[imagelist[i]*xDim*yDim + j];
					}

					clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
					clRotScale->Enqueue(fullWorkSize);
					clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

					for(int j = 0; j < yDim;j++)
						for(int i = 0; i < xDim;i++)
						{
							image[i+j*xDim].s[0] = copyImage[i + (j)*xDim].real();
							image[i+j*xDim].s[1] = copyImage[i + (j)*xDim].imag();
						}
				}
				else
				{
					for(int j = 0 ; j < xDim*yDim ; j++)
					{
						image[j].s[0] = seriesdata[imagelist[i]*xDim*yDim + j];
						image[j].s[1] = 0;
					}
				}

				clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
				clPadCrop->Enqueue(globalWorkSize);

				image.clear();

				// Get FFT of this image
				FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

				clWaveTransferFunction->SetArgT(10,defocus);
				clWaveTransferFunctionMinus->SetArgT(10,defocus);

				clWaveTransferFunction->Enqueue(globalWorkSize);
				clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

				// add to W,W-,V,T,U 
				clWienerW->Enqueue(globalWorkSize);
				clWienerWMinus->Enqueue(globalWorkSize);
				clWienerV->Enqueue(globalWorkSize);
				clWienerT->Enqueue(globalWorkSize);
				clWienerU->Enqueue(globalWorkSize);
			}

			// make restored
			clMakeRestored->Enqueue(globalWorkSize);
			 
			// inverse fft for display
			FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);
			
			//Show Each reconstruction
			//Utility::PrintCLMemToImage(clRestored,"Restored",sizeX,sizeY,clFloat2,clq);
			
			// PCF with the reconstruction
			// Adjust reconstruction or do PCPCF?
			// TODO: at the moment reconstruction is in reference plane and images registered by PCPCF with this plane....
			// Also do i need to get modulus and FFT or not? - seems to work OK
			
			// Work out what focus difference is expected between these 2 images.
			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone


						
			// Rotate and Scale this image by expected defocus..
			if(options.rotscale)
			{
				// Might be bullshit
				float expecteddifference = (currentimage-referenceimage) * options.focalstep;

				clRotScale->SetArgT(6,expecteddifference);
				// Copy image into fullimage
				for(int j = 0 ; j < xDim*yDim ; j++)
				{
					copyImage[j] = seriesdata[currentimage*xDim*yDim + j];

				}

				clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
				clRotScale->Enqueue(fullWorkSize);
				clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

				//TODO: Could be clever and keep swapping the vectors to only load one image each time
				for(int j = 0; j < sizeY;j++)
					for(int i = 0; i < sizeX;i++)
					{
						dataTwo[i+j*sizeX] = copyImage[i + iLeft + (j+iTop)*xDim];
					}
			} 
			else 
			{	
				for(int j = 0; j < sizeY;j++)
					for(int i = 0; i < sizeX;i++)
					{
						dataTwo[i+j*sizeX] = seriesdata[currentimage*xDim*yDim + i + iLeft + (j+iTop)*xDim];
					}
			}

			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 0, NULL, NULL );

			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			float expecteddifference = (currentimage-referenceimage) * options.focalstep;

			// Get number of trial steps
			// Assume one if not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			for(int trial = 0; trial < numberoftrials; trial++)
			{
				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected = searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep + ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);
				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				// Could consider having this displayed to see roughly how well it is working...
				//Utility::PrintCLMemToImage(clFFTImage1,"FFT1",sizeX,sizeY,clFloat2,clq);
				//Utility::PrintCLMemToImage(clImage1,"PCPCF",sizeX,sizeY,clFloat2,clq);
				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(currentimage)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY);

				trialxshifts[trial] = xShift;
				trialyshifts[trial] = yShift;
				trialsubxshifts[trial] = xShift + subXShift;
				trialsubyshifts[trial] = yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;


			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[currentimage] = trialxshifts[trial];
					yShiftVals[currentimage] = trialyshifts[trial];
					subXShifts[currentimage] = trialsubxshifts[trial];
					subYShifts[currentimage] = trialsubyshifts[trial];
					defocusshifts[currentimage] = trialdefocus[trial];
				}
			}

			imagelist.push_back(currentimage);

			delete[] trialxshifts;
			delete[] trialyshifts;
			delete[] trialsubxshifts;
			delete[] trialsubyshifts;
			delete[] trialdefocus;
			delete[] peakheights;
		}
		
		// to get to next image..
		gonext++;
		nextdir *=-1;
		// image = previous +(gonext*nextdir);
		// checkbounds then if not available do again.
		
	}

	/* Now rebuild reconstruction including the final registered image */
	// All shifts have to be stored relative to the reference image.
	float maxnegx = 0;
	float maxnegy = 0;
	float maxposx = 0;
	float maxposy = 0;

	// Could just loop over this list of currently registered images instead...
	for(int i = 1; i <= numberOfImages; i++)
	{
		if(subXShifts[i-1] < maxnegx)
			maxnegx = subXShifts[i-1];
		if(subYShifts[i-1] < maxnegy)
			maxnegy = subYShifts[i-1];
		if(subXShifts[i-1] > maxposx)
			maxposx = subXShifts[i-1];
		if(subYShifts[i-1] > maxposy)
			maxposy = subYShifts[i-1];
	}

	// Determine amount to pad on either direction
	int padTop = 0;
	int padLeft = 0;
	int padRight = 0;
	int padBottom = 0;

	if (abs(maxnegx)-iLeft > 0)
		padLeft = ceil(abs(maxnegx)-iLeft);

	if (maxposx-iRight > 0)
		padRight = ceil(maxposx-iRight);

	if (abs(maxnegy)-iTop > 0)
		padTop = ceil(abs(maxnegy)-iTop);

	if (maxposy-iBottom > 0)
		padBottom = ceil(maxposy-iBottom);


	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]];
	
		// Rotate and Scale this image by expected defocus..
		if(options.rotscale)
		{
			clRotScale->SetArgT(6,defocus);
			// Copy image into fullimage
			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				copyImage[j] = seriesdata[imagelist[i]*xDim*yDim + j];
			}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			//TODO: Could be clever and keep swapping the vectors to only load one image each time
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					image[i+j*xDim].s[0] = copyImage[i  + (j)*xDim].real();
					image[i+j*xDim].s[1] = copyImage[i  + (j)*xDim].imag();
				}
		} 
		else 
		{
			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				image[j].s[0] = seriesdata[imagelist[i]*xDim*yDim + j];
				image[j].s[1] = 0;
			}
		}

		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();

		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);

		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);
			 
	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Restored EW",sizeX,sizeY,clFloat2,clq);

	// Now all should be registered - produce the final restored image and any other relevant data like Shift Graphs, and Drift Corrected Stack
	Gatan::DM::Image driftcorrected = DM::RealImage("Drift Corrected Stack",4,sizeX,sizeY,numberOfImages);
	Gatan::PlugIn::ImageDataLocker driftlocker = Gatan::PlugIn::ImageDataLocker(driftcorrected);
	float* driftstack = (float*) driftlocker.get();

	std::vector<cl_float2> cropimage(sizeX*sizeY);

	for(int i = 0 ; i < numberOfImages ; i++)
	{
			
		int padTop = 0;
		int padLeft = 0;
		int padRight = 0;
		int padBottom = 0;

		if (abs(maxnegx)-iLeft > 0)
			padLeft = ceil(abs(maxnegx)-iLeft);

		if (maxposx-iRight > 0)
			padRight = ceil(maxposx-iRight);

		if (abs(maxnegy)-iTop > 0)
			padTop = ceil(abs(maxnegy)-iTop);

		if (maxposy-iBottom > 0)
			padBottom = ceil(maxposy-iBottom);

		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[i]);
		clPadCrop->SetArgT(9,subYShifts[i]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);
		/*
		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = seriesdata[i * xDim * yDim + j];
			image[j].s[1] = 0;
		}
		*/

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[i];
				
		// Rotate and Scale this image by expected defocus..
		if(options.rotscale)
		{
			clRotScale->SetArgT(6,defocus);
			// Copy image into fullimage
			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				copyImage[j] = seriesdata[i*xDim*yDim + j];
			}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			//TODO: Could be clever and keep swapping the vectors to only load one image each time
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					image[i+j*xDim].s[0] = copyImage[i + (j)*xDim].real();
					image[i+j*xDim].s[1] = copyImage[i + (j)*xDim].imag();
				}

		} 
		else 
		{
			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				image[j].s[0] = seriesdata[i*xDim*yDim + j];
				image[j].s[1] = 0;
			}
		}


		// Copy full image to GPU and run crop kernel
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		// Copy image data back to host
		clEnqueueReadBuffer(clq->cmdQueue,clImage1,CL_TRUE,0,sizeX*sizeY*sizeof(cl_float2),&cropimage[0],0,NULL,NULL);
		
		// Copy to drift corrected stack
		for(int k = 0 ; k < sizeX*sizeY ; k++)
		{
			driftstack[i*sizeX*sizeY + k] = cropimage[k].s[0];
		}

		image.clear();
	}

	cropimage.clear();
	driftlocker.~ImageDataLocker();

	// Display drift corrected
	driftcorrected.GetOrCreateImageDocument().Show();

	dataOne.clear();
	dataTwo.clear();



	// Lots of these have been moved later on because I added extra steps that require them.
	
	//clReleaseMemObject(clImage1);
	//clReleaseMemObject(clFFTImage1);
	//clReleaseMemObject(fullImage);
	//clReleaseMemObject(clw);
	//clReleaseMemObject(clwminus);
	//clReleaseMemObject(clW);
	//clReleaseMemObject(clWminus);
	//clReleaseMemObject(clU);
	//clReleaseMemObject(clV);
	//clReleaseMemObject(clT);
	//clReleaseMemObject(clRestored);
	//clPadCrop->~clKernel();
	//clWienerW->~clKernel();
	//clWienerWMinus->~clKernel();
	//clWienerV->~clKernel();
	//clWienerT->~clKernel();
	//clWienerU->~clKernel();
	//clWaveTransferFunction->~clKernel();
	//clWaveTransferFunctionMinus->~clKernel();
	
	// Don't need these anymore
	PCPCF->~clKernel();
	clReleaseMemObject(clFFTImage2);
	clReleaseMemObject(clPCPCFResult);
	clReleaseMemObject(clImage2);

	// Now start working on finding the actual astigmatism and defocus from the reconstruction....

	// Make additional kernels
	clKernel* clGetQ				= new clKernel(getQsource, context, cldev, "clCalculateQ", clq);
	clKernel* clSumReduction		= new clKernel(sumReductionsource,context, cldev, "clSumReduction", clq);
	clKernel* clMinusWavefunction	= new clKernel(minuswavefunctionsource, context, cldev, "clMinusWavefunction", clq);
	clKernel* clCalculatePCI		= new clKernel(getPCIsource,context, cldev, "clCalculatePCI", clq);
	clGetQ->BuildKernel();
	clSumReduction->BuildKernel();
	clMinusWavefunction->BuildKernel();
	clCalculatePCI->BuildKernel();

	// New memory required
	cl_mem clRestoredMinus	= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clQ				= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCI			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIC			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIM			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIS			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);

	// Set Kernel Arguments
	clGetQ->SetArgT(0,clW);
	clGetQ->SetArgT(1,clWminus);
	clGetQ->SetArgT(2,clV);
	clGetQ->SetArgT(3,clQ);
	clGetQ->SetArgT(4,sizeX);
	clGetQ->SetArgT(5,sizeY);

	clMinusWavefunction->SetArgT(0,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clMinusWavefunction->SetArgT(1,clRestoredMinus);
	clMinusWavefunction->SetArgT(2,sizeX);
	clMinusWavefunction->SetArgT(3,sizeY);

	clGetQ->Enqueue(globalWorkSize);
	clMinusWavefunction->Enqueue(globalWorkSize);

	// Get sum of Q by my awesome reduction kernel ¬_¬ 
	int totalSize = sizeX*sizeY;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	cl_mem clSumOutput = clCreateBuffer(context,CL_MEM_READ_WRITE,nGroups*sizeof(std::complex<float>),0,&status);
	std::vector< std::complex< float > > sumQs( nGroups );

	clSumReduction->SetArgT(0,clQ);
	clSumReduction->SetArgT(1,clSumOutput);
	clSumReduction->SetArgT(2,totalSize);
	clSumReduction->SetArgLocalMemory(3,256,clFloat2);

	clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 0, NULL, NULL );

	
	// Find out which numbers to read back
	float sumQ2 = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sumQ2 += sumQs[i].real();
	}

	// Trial lots of C1 Values
	int trials2 = 50;

	float* Fpci = new float[trials2];
	float* Fpcic = new float[trials2];
	float* Fpcis = new float[trials2];
	float* Fpcim = new float[trials2];
	float* Fpci2 = new float[trials2];

	float A1r = 0.0f;
	float A1i = 0.0f;
	float zero2 = 0.0f;
	
	float StartDf = -75;
	float StepDf = 3;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	clCalculatePCI->SetArgT(0,clQ);
	clCalculatePCI->SetArgT(1,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clCalculatePCI->SetArgT(2,clRestoredMinus);
	clCalculatePCI->SetArgT(3,clxFrequencies);
	clCalculatePCI->SetArgT(4,clyFrequencies);
	clCalculatePCI->SetArgT(5,clPCI);
	clCalculatePCI->SetArgT(6,clPCIC);
	clCalculatePCI->SetArgT(7,clPCIM);
	clCalculatePCI->SetArgT(8,clPCIS);
	clCalculatePCI->SetArgT(9,sizeX);
	clCalculatePCI->SetArgT(10,sizeY);
	clCalculatePCI->SetArgT(12,zero2);
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);
	clCalculatePCI->SetArgT(15,wavelength);
	clCalculatePCI->SetArgT(16,Abb.kmax);
	

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
 
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);
	
		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.
		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpci[trial] = fpcitrial/sumQ2;
		Fpcic[trial] = fpcictrial/sumQ2;
		Fpcis[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcim[trial] = fpcimtrial/sumQ2;
	}

	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = StartDf + trial*StepDf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial <trials2 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess)
		);


	DM::Image pcigraph = DM::RealImage("PCI",4,trials2);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial <trials2 ; trial++)
	{
		pciData[trial] = Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	DM::UpdateImage(pcigraph);


	delete[] Fpci;
	delete[] Fpci2;
	delete[] Fpcic;
	delete[] Fpcim;
	delete[] Fpcis;

	// First PCI Stage Complete

	// Redo more accurate PCI in correct area 
	// Trial lots of C1 Values

	int trials3 = 50;
	
	float* Fpciref = new float[trials3];
	float* Fpcicref = new float[trials3];
	float* Fpcisref = new float[trials3];
	float* Fpcimref = new float[trials3];
	float* Fpci2ref = new float[trials3];

	A1r = A1rGuess;
	A1i = A1iGuess;
	
	StartDf = DefocusGuess-10;
	StepDf = 21.0f/trials2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	// Set Kernel Arguments - DONE ALREADY NOT CHANGED
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);


	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);

		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.

		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );
		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpciref[trial] = fpcitrial/sumQ2;
		Fpcicref[trial] = fpcictrial/sumQ2;
		Fpcisref[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcimref[trial] = fpcimtrial/sumQ2;

	}

	float maxheighta = 0;
	float C1aa = 0;
	int besttriala = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		if(Fpcimref[trial] > maxheighta)
		{
			maxheighta = Fpcimref[trial];
			C1aa = StartDf + trial*StepDf;
			besttriala = trial;
		}
	}

	float astigaxisa = atan2(Fpcisref[besttriala],Fpcicref[besttriala])/2;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		Fpci2ref[trial] = Fpciref[trial] - Fpcicref[trial] * cos(2*astigaxisa) - Fpcisref[trial]*sin(2*astigaxisa);
	}

	float maxheight2a = 0;
	float C1ba = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess2),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess2),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess2)
		);

	// Print out some results and graphs

	DM::Image pcigraph2 = DM::RealImage("PCI - Refined",4,trials3);
	Gatan::PlugIn::ImageDataLocker pciLocker2(pcigraph2);
	float* pciData2 = (float*) pciLocker2.get();

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		pciData2[trial]=Fpciref[trial];
	}

	pciLocker2.~ImageDataLocker();
	pcigraph2.GetOrCreateImageDocument().Show();
	pcigraph2.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	

	delete[] Fpciref;
	delete[] Fpcicref;
	delete[] Fpcisref;
	delete[] Fpcimref;
	delete[] Fpci2ref;

	// Second PCI Stage Complete


	// Then do restoration again at correct focal/astig value

	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..

		std::vector<cl_float2> image(xDim*yDim);

		// Want to rotate and scale correct this image -- do rotations from reference defocus not absolute defocus or it will (probably) affect shifts.

		float expecteddifference = (imagelist[i]-referenceimage) * options.focalstep;


		// Rotate and Scale this image by expected defocus..
		if(options.rotscale)
		{
			clRotScale->SetArgT(6,expecteddifference);
			// Copy image into fullimage
			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				copyImage[j] = seriesdata[imagelist[i]*xDim*yDim + j];

			}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			//TODO: Could be clever and keep swapping the vectors to only load one image each time
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					image[i+j*xDim].s[0] = copyImage[i  + (j)*xDim].real();
					image[i+j*xDim].s[1] = copyImage[i  + (j)*xDim].imag();
				}

		} else {

			for(int j = 0 ; j < xDim*yDim ; j++)
			{
				image[j].s[0] = seriesdata[imagelist[i]*xDim*yDim + j];
				image[j].s[1] = 0;
			}
		}








		/*
		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = seriesdata[imagelist[i] * xDim * yDim + j];
			image[j].s[1] = 0;
		}*/

		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();


		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]] + DefocusGuess2;

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);
		clWaveTransferFunction->SetArgT(8,A1rGuess2);
		clWaveTransferFunctionMinus->SetArgT(8,A1rGuess2);
		clWaveTransferFunction->SetArgT(9,A1iGuess2);
		clWaveTransferFunctionMinus->SetArgT(9,A1iGuess2);
		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);

	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Final Restored EW",sizeX,sizeY,clFloat2,clq);

	// New Image added for doing rotation and scale correcting
	copyImage.clear();

	clReleaseMemObject(clImage1);
	clReleaseMemObject(clFFTImage1);
	clReleaseMemObject(fullImage);
	clReleaseMemObject(clw);
	clReleaseMemObject(clwminus);
	clReleaseMemObject(clW);
	clReleaseMemObject(clWminus);
	clReleaseMemObject(clU);
	clReleaseMemObject(clV);
	clReleaseMemObject(clT);
	clReleaseMemObject(clRestored);
	clReleaseMemObject(clRestoredMinus);
	clReleaseMemObject(clQ);
	clReleaseMemObject(clPCI);
	clReleaseMemObject(clPCIC);
	clReleaseMemObject(clPCIM);
	clReleaseMemObject(clPCIS);
	clReleaseMemObject(clSumOutput);

	clPadCrop->~clKernel();
	clWienerW->~clKernel();
	clWienerWMinus->~clKernel();
	clWienerV->~clKernel();
	clWienerT->~clKernel();
	clWienerU->~clKernel();
	clWaveTransferFunction->~clKernel();
	clWaveTransferFunctionMinus->~clKernel();
	clMakeRestored->~clKernel();

	clGetQ->~clKernel();
	clSumReduction->~clKernel();
	clMinusWavefunction->~clKernel();
	clCalculatePCI->~clKernel();
	
};


// Maybe a third version that just prerotates and scales the entire stack then just runs the original reconstruction on that?
// Sometimes getting funny results after a certain point in code, not sure which step is failing but it starts giving same defocus for every image which should be impossible...
void PCPCFWrapper3(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb)
{
	//Utility::SetResultWindow("Should Align = "+boost::lexical_cast<std::string>(options.determinefocus)+"\n");
	//Utility::SetResultWindow("Steps = "+boost::lexical_cast<std::string>(options.steps)+"\n");

	std::vector<int> imagelist;

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	cl_int status;
	// Construct PCPCF Kernel

	clKernel* PCPCF = new clKernel(pcpcfSource,context,cldev,"clPCPCF",clq);
	PCPCF->BuildKernel();

	// Setup required memory for PCPCF subroutine

	cl_mem clImage1			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clImage2			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage1		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage2		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clPCPCFResult	= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);

	cl_mem fullImage		= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem rotScaleImage	= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);

	cl_mem clW				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clWminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clw				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clwminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clT				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clU				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clV				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clRestored		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);

	int init = 1;
	int noinit = 0;

	PCPCF->SetArgT(0,clFFTImage1);
	PCPCF->SetArgT(1,clFFTImage2);
	PCPCF->SetArgT(2,clPCPCFResult);
	PCPCF->SetArgT(3,clxFrequencies);
	PCPCF->SetArgT(4,clyFrequencies);
	PCPCF->SetArgT(5,sizeX);
	PCPCF->SetArgT(6,sizeY);
	PCPCF->SetArgT(8,wavelength);
	PCPCF->SetArgT(9,options.pcpcfkmax);

	/* Make all the other kernels */
	clKernel* clPadCrop = new clKernel(PadCropSource,context,cldev,"clPadCrop",clq);
	clPadCrop->BuildKernel();

	clPadCrop->SetArgT(0,fullImage);
	clPadCrop->SetArgT(1,clImage1);
	clPadCrop->SetArgT(2,xDim); // long not int
	clPadCrop->SetArgT(3,yDim); // long not int
	clPadCrop->SetArgT(10,sizeX);
	clPadCrop->SetArgT(11,sizeY);
	clPadCrop->SetArgT(12,iTop);
	clPadCrop->SetArgT(13,iLeft);

	clKernel* clWaveTransferFunction = new clKernel(wavetransferfunctionsource,context,cldev,"clWaveTransferFunction",clq);
	clWaveTransferFunction->BuildKernel();

	clWaveTransferFunction->SetArgT(0,clw);
	clWaveTransferFunction->SetArgT(1,clxFrequencies);
	clWaveTransferFunction->SetArgT(2,clyFrequencies);
	clWaveTransferFunction->SetArgT(3,sizeX);
	clWaveTransferFunction->SetArgT(4,sizeY);
	clWaveTransferFunction->SetArgT(5,wavelength);
	clWaveTransferFunction->SetArgT(6,Abb.beta);
	clWaveTransferFunction->SetArgT(7,Abb.delta);
	clWaveTransferFunction->SetArgT(8,Abb.A1r);
	clWaveTransferFunction->SetArgT(9,Abb.A1i);
	clWaveTransferFunction->SetArgT(11,Abb.Cs);
	clWaveTransferFunction->SetArgT(12,Abb.kmax);

	clKernel* clWaveTransferFunctionMinus = new clKernel(wavetransferfunctionminussource,context,cldev,"clWaveTransferFunctionMinus",clq);
	clWaveTransferFunctionMinus->BuildKernel();

	clWaveTransferFunctionMinus->SetArgT(0,clwminus);
	clWaveTransferFunctionMinus->SetArgT(1,clxFrequencies);
	clWaveTransferFunctionMinus->SetArgT(2,clyFrequencies);
	clWaveTransferFunctionMinus->SetArgT(3,sizeX);
	clWaveTransferFunctionMinus->SetArgT(4,sizeY);
	clWaveTransferFunctionMinus->SetArgT(5,wavelength);
	clWaveTransferFunctionMinus->SetArgT(6,Abb.beta);
	clWaveTransferFunctionMinus->SetArgT(7,Abb.delta);
	clWaveTransferFunctionMinus->SetArgT(8,Abb.A1r);
	clWaveTransferFunctionMinus->SetArgT(9,Abb.A1i);
	clWaveTransferFunctionMinus->SetArgT(11,Abb.Cs);
	clWaveTransferFunctionMinus->SetArgT(12,Abb.kmax);

	clKernel* clWienerW = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerW->BuildKernel();

	clWienerW->SetArgT(0,clW);
	clWienerW->SetArgT(1,clw);
	clWienerW->SetArgT(2,sizeX);
	clWienerW->SetArgT(3,sizeY);
	clWienerW->SetArgT(4,init);

	clKernel* clWienerWMinus = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerWMinus->BuildKernel();

	clWienerWMinus->SetArgT(0,clWminus);
	clWienerWMinus->SetArgT(1,clwminus);
	clWienerWMinus->SetArgT(2,sizeX);
	clWienerWMinus->SetArgT(3,sizeY);
	clWienerWMinus->SetArgT(4,init);

	clKernel* clWienerV = new clKernel(wienervsource,context,cldev,"clWienerV",clq);
	clWienerV->BuildKernel();

	clWienerV->SetArgT(0,clV);
	clWienerV->SetArgT(1,clw);
	clWienerV->SetArgT(2,clwminus);
	clWienerV->SetArgT(3,sizeX);
	clWienerV->SetArgT(4,sizeY);
	clWienerV->SetArgT(5,init);

	clKernel* clWienerT = new clKernel(wienertsource,context,cldev,"clWienerT",clq);
	clWienerT->BuildKernel();

	clWienerT->SetArgT(0,clT);
	clWienerT->SetArgT(1,clw);
	clWienerT->SetArgT(2,clFFTImage1);
	clWienerT->SetArgT(3,sizeX);
	clWienerT->SetArgT(4,sizeY);
	clWienerT->SetArgT(5,init);

	clKernel* clWienerU = new clKernel(wienerusource,context,cldev,"clWienerU",clq);
	clWienerU->BuildKernel();

	clWienerU->SetArgT(0,clU);
	clWienerU->SetArgT(1,clwminus);
	clWienerU->SetArgT(2,clFFTImage1);
	clWienerU->SetArgT(3,sizeX);
	clWienerU->SetArgT(4,sizeY);
	clWienerU->SetArgT(5,init);
	
	clKernel* clMakeRestored = new clKernel(makerestoredsource,context,cldev,"clMakeRestored",clq);
	clMakeRestored->BuildKernel();

	clMakeRestored->SetArgT(0,clW);
	clMakeRestored->SetArgT(1,clWminus);
	clMakeRestored->SetArgT(2,clV);
	clMakeRestored->SetArgT(3,clT);
	clMakeRestored->SetArgT(4,clU);
	clMakeRestored->SetArgT(5,clFFTImage1);
	clMakeRestored->SetArgT(6,sizeX);
	clMakeRestored->SetArgT(7,sizeY);

	clKernel* clRotScale = new clKernel(RotScaleSource,context,cldev,"clRotScale",clq);
	clRotScale->BuildKernel();

	if(options.rotscale)
	{
		clRotScale->SetArgT(0,fullImage);
		clRotScale->SetArgT(1,rotScaleImage);
		clRotScale->SetArgT(2,xDim);
		clRotScale->SetArgT(3,yDim);
		clRotScale->SetArgT(4,options.magcal);
		clRotScale->SetArgT(5,options.rotcal);
	}

	size_t* fullWorkSize = new size_t[3];
	
	// There are 'elements' work items
	fullWorkSize[0] = xDim;
	fullWorkSize[1] = yDim;
	fullWorkSize[2] = 1;

	std::vector<std::complex<float>> copyImage(xDim*yDim) ;
	std::vector<float> rotscaledseries(xDim*yDim*numberOfImages);

	// Make a new rotation aligned stack based on the expected defocus step.
	// To prevent edges going missing should start at end that should be smallest so all images are only scaled upwards....


	
	// Rotate and Scale this image by expected defocus..
	if(options.rotscale)
	{
		int startingimage = 0;

		if(options.magcal > 1)
		{
			startingimage = 0;
		}
		else
		{
			startingimage = numberOfImages -1;
		}


		for(int im = 0; im < numberOfImages ; im++ )
		{
			float expecteddifference = (im-startingimage)*options.focalstep;
			clRotScale->SetArgT(6,expecteddifference);

			// Copy image into fullimage
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					copyImage[i+j*xDim] = seriesdata[im*xDim*yDim + i + (j)*xDim];
				}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					rotscaledseries[im*xDim*yDim +i + j*xDim] = copyImage[i + (j)*xDim].real();
				}
		}		
	} 
	else
	{
		// Fill new series with old series data.
		for(int im = 0; im < numberOfImages ; im++ )
			for(int j = 0; j < yDim;j++)
					for(int i = 0; i < xDim;i++)
					{
						rotscaledseries[im*xDim*yDim +i + j*xDim] = seriesdata[im*xDim*yDim +i + j*xDim];
					}
	}

	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[3];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	globalWorkSize[2] = 1;



	std::vector< std::complex< float > > dataOne( sizeX*sizeY );
	std::vector< std::complex< float > > dataTwo( sizeX*sizeY );

	int referenceimage = options.reference;
	int gonext = 1;
	int nextdir = -1;
	
	// First three have to be done via standard methods.

	// Initialise image numbers
	int imageone = referenceimage;
	int imagetwo = referenceimage;
	int currentimage;

	imagelist.push_back(referenceimage);



	// loop over all of the images
	for(int n = 0; n < numberOfImages - 1; n++)
	{
		
		 // For RECONPCPCF part

		// First 2 registrations cannot be done against reconstruction
		// have to be done against reference image...
		if(n < 2)
		{
			imageone = referenceimage;
			imagetwo += (gonext*nextdir);

			while(! (imagetwo < numberOfImages && imagetwo >=0)) // not an actual image 
			{
				gonext++;
				nextdir *=-1;
				imagetwo += (gonext*nextdir);
			}

			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = (imagetwo-imageone) * options.focalstep;

			for(int j = 0; j < sizeY;j++)
				for(int i = 0; i < sizeX;i++)
				{
					dataOne[i+j*sizeX] = rotscaledseries[imageone*xDim*yDim + i + iLeft + (j+iTop)*xDim];
					dataTwo[i+j*sizeX] = rotscaledseries[imagetwo*xDim*yDim + i + iLeft + (j+iTop)*xDim];
				}

	
			clEnqueueWriteBuffer( clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
						0, NULL, NULL );

			FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			for(int trial = 0; trial < numberoftrials; trial++)
			{

				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected + searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep + ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);

				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(imagetwo)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY);

				trialxshifts[trial] = xShift;
				trialyshifts[trial] = yShift;
				trialsubxshifts[trial] = xShift + subXShift;
				trialsubyshifts[trial] = yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;
			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[imagetwo] = trialxshifts[trial];
					yShiftVals[imagetwo] = trialyshifts[trial];
					subXShifts[imagetwo] = trialsubxshifts[trial];
					subYShifts[imagetwo] = trialsubyshifts[trial];
					defocusshifts[imagetwo] = trialdefocus[trial];
				}
			}

			imagelist.push_back(imagetwo);
		}

		if(n==2)
		{
			// To make sure it is initialised from end of first section
			currentimage = imagetwo;
		}
	
		if(n>=2)
		{
			// Now we can perform reconstruction and register to the reconstruction instead of standard pcpcf
			currentimage += (gonext*nextdir);

			while(! (currentimage < numberOfImages && currentimage >=0)) // not an actual image
			{
				gonext++;
				nextdir *=-1;
				currentimage += (gonext*nextdir);
			}

			// Need to reconstruct the previous images based on there focus and shifts relative to reference image...

			
			// Need to extract correct area from each image to perform reconstruction
			// 2 Possibilites - If an area of the correct size is not available then give an error...
			// or pad the area with zeroes for ALL images - need to first determine the extents of available area before we extract from any of the images.

			// Find maximum negative and positive shifts
			// All shifts have to be stored relative to the reference image.
			float maxnegx = 0;
			float maxnegy = 0;
			float maxposx = 0;
			float maxposy = 0;

			// Could just loop over this list of currently registered images instead...
			for(int i = 1; i <= imagelist.size(); i++)
			{
				if(subXShifts[imagelist[i-1]] < maxnegx)
					maxnegx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] < maxnegy)
					maxnegy = subYShifts[imagelist[i-1]];
				if(subXShifts[imagelist[i-1]] > maxposx)
					maxposx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] > maxposy)
					maxposy = subYShifts[imagelist[i-1]];
			}

			// Determine amount to pad on either direction
			int padTop = 0;
			int padLeft = 0;
			int padRight = 0;
			int padBottom = 0;

			if (abs(maxnegx)-iLeft > 0)
				padLeft = ceil(abs(maxnegx)-iLeft);

			if (maxposx-iRight > 0)
				padRight = ceil(maxposx-iRight);

			if (abs(maxnegy)-iTop > 0)
				padTop = ceil(abs(maxnegy)-iTop);

			if (maxposy-iBottom > 0)
				padBottom = ceil(maxposy-iBottom);	

			for(int i = 0 ; i < imagelist.size() ; i++)
			{
				if(i == 0) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,init);
					clWienerWMinus->SetArgT(4,init);
					clWienerV->SetArgT(5,init);
					clWienerT->SetArgT(5,init);
					clWienerU->SetArgT(5,init);
				}
				if(i == 1) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,noinit);
					clWienerWMinus->SetArgT(4,noinit);
					clWienerV->SetArgT(5,noinit);
					clWienerT->SetArgT(5,noinit);
					clWienerU->SetArgT(5,noinit);
				}
				
				// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
				clPadCrop->SetArgT(4,padLeft);
				clPadCrop->SetArgT(5,padRight);
				clPadCrop->SetArgT(6,padTop);
				clPadCrop->SetArgT(7,padBottom);
				clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
				clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
				// Copy correct image into full image..
				std::vector<cl_float2> image(xDim*yDim);

				// Now add image to reconstruction - Make WTF - WTFminus			
				float defocus = defocusshifts[imagelist[i]];

				for(int j = 0 ; j < xDim*yDim ; j++)
				{
					image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
					image[j].s[1] = 0;
				}


				clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
				clPadCrop->Enqueue(globalWorkSize);

				image.clear();

				// Get FFT of this image
				FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

				clWaveTransferFunction->SetArgT(10,defocus);
				clWaveTransferFunctionMinus->SetArgT(10,defocus);

				clWaveTransferFunction->Enqueue(globalWorkSize);
				clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

				// add to W,W-,V,T,U 
				clWienerW->Enqueue(globalWorkSize);
				clWienerWMinus->Enqueue(globalWorkSize);
				clWienerV->Enqueue(globalWorkSize);
				clWienerT->Enqueue(globalWorkSize);
				clWienerU->Enqueue(globalWorkSize);
			}

			// make restored
			clMakeRestored->Enqueue(globalWorkSize);
			 
			// inverse fft for display
			FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);
			
			//Show Each reconstruction
			//Utility::PrintCLMemToImage(clRestored,"Restored",sizeX,sizeY,clFloat2,clq);
			
			// PCF with the reconstruction
			// Adjust reconstruction or do PCPCF?
			// TODO: at the moment reconstruction is in reference plane and images registered by PCPCF with this plane....
			// Also do i need to get modulus and FFT or not? - seems to work OK
			
			// Work out what focus difference is expected between these 2 images.
			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone


						
	
			for(int j = 0; j < sizeY;j++)
				for(int i = 0; i < sizeX;i++)
				{
					dataTwo[i+j*sizeX] = rotscaledseries[currentimage*xDim*yDim + i + iLeft + (j+iTop)*xDim];
				}

			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 0, NULL, NULL );

			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			float expecteddifference = (currentimage-referenceimage) * options.focalstep;

			// Get number of trial steps
			// Assume one if not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			for(int trial = 0; trial < numberoftrials; trial++)
			{
				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected = searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep + ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);
				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				// Could consider having this displayed to see roughly how well it is working...
				//Utility::PrintCLMemToImage(clFFTImage1,"FFT1",sizeX,sizeY,clFloat2,clq);
				//Utility::PrintCLMemToImage(clImage1,"PCPCF",sizeX,sizeY,clFloat2,clq);
				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(currentimage)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY);

				trialxshifts[trial] = xShift;
				trialyshifts[trial] = yShift;
				trialsubxshifts[trial] = xShift + subXShift;
				trialsubyshifts[trial] = yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;


			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[currentimage] = trialxshifts[trial];
					yShiftVals[currentimage] = trialyshifts[trial];
					subXShifts[currentimage] = trialsubxshifts[trial];
					subYShifts[currentimage] = trialsubyshifts[trial];
					defocusshifts[currentimage] = trialdefocus[trial];
				}
			}

			imagelist.push_back(currentimage);

			delete[] trialxshifts;
			delete[] trialyshifts;
			delete[] trialsubxshifts;
			delete[] trialsubyshifts;
			delete[] trialdefocus;
			delete[] peakheights;
		}
		
		// to get to next image..
		gonext++;
		nextdir *=-1;
		// image = previous +(gonext*nextdir);
		// checkbounds then if not available do again.
		
	}

	/* Now rebuild reconstruction including the final registered image */
	// All shifts have to be stored relative to the reference image.
	float maxnegx = 0;
	float maxnegy = 0;
	float maxposx = 0;
	float maxposy = 0;

	// Could just loop over this list of currently registered images instead...
	for(int i = 1; i <= numberOfImages; i++)
	{
		if(subXShifts[i-1] < maxnegx)
			maxnegx = subXShifts[i-1];
		if(subYShifts[i-1] < maxnegy)
			maxnegy = subYShifts[i-1];
		if(subXShifts[i-1] > maxposx)
			maxposx = subXShifts[i-1];
		if(subYShifts[i-1] > maxposy)
			maxposy = subYShifts[i-1];
	}

	// Determine amount to pad on either direction
	int padTop = 0;
	int padLeft = 0;
	int padRight = 0;
	int padBottom = 0;

	if (abs(maxnegx)-iLeft > 0)
		padLeft = ceil(abs(maxnegx)-iLeft);

	if (maxposx-iRight > 0)
		padRight = ceil(maxposx-iRight);

	if (abs(maxnegy)-iTop > 0)
		padTop = ceil(abs(maxnegy)-iTop);

	if (maxposy-iBottom > 0)
		padBottom = ceil(maxposy-iBottom);


	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]];

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();

		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);

		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);
			 
	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Restored EW",sizeX,sizeY,clFloat2,clq);

	// Now all should be registered - produce the final restored image and any other relevant data like Shift Graphs, and Drift Corrected Stack
	Gatan::DM::Image driftcorrected = DM::RealImage("Drift Corrected Stack",4,sizeX,sizeY,numberOfImages);
	Gatan::PlugIn::ImageDataLocker driftlocker = Gatan::PlugIn::ImageDataLocker(driftcorrected);
	float* driftstack = (float*) driftlocker.get();

	std::vector<cl_float2> cropimage(sizeX*sizeY);

	for(int i = 0 ; i < numberOfImages ; i++)
	{
			
		int padTop = 0;
		int padLeft = 0;
		int padRight = 0;
		int padBottom = 0;

		if (abs(maxnegx)-iLeft > 0)
			padLeft = ceil(abs(maxnegx)-iLeft);

		if (maxposx-iRight > 0)
			padRight = ceil(maxposx-iRight);

		if (abs(maxnegy)-iTop > 0)
			padTop = ceil(abs(maxnegy)-iTop);

		if (maxposy-iBottom > 0)
			padBottom = ceil(maxposy-iBottom);

		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[i]);
		clPadCrop->SetArgT(9,subYShifts[i]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[i];
			
		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[i*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		// Copy full image to GPU and run crop kernel
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		// Copy image data back to host
		clEnqueueReadBuffer(clq->cmdQueue,clImage1,CL_TRUE,0,sizeX*sizeY*sizeof(cl_float2),&cropimage[0],0,NULL,NULL);
		
		// Copy to drift corrected stack
		for(int k = 0 ; k < sizeX*sizeY ; k++)
		{
			driftstack[i*sizeX*sizeY + k] = cropimage[k].s[0];
		}

		image.clear();
	}

	cropimage.clear();
	driftlocker.~ImageDataLocker();

	// Display drift corrected
	driftcorrected.GetOrCreateImageDocument().Show();

	dataOne.clear();
	dataTwo.clear();



	// Lots of these have been moved later on because I added extra steps that require them.
	
	//clReleaseMemObject(clImage1);
	//clReleaseMemObject(clFFTImage1);
	//clReleaseMemObject(fullImage);
	//clReleaseMemObject(clw);
	//clReleaseMemObject(clwminus);
	//clReleaseMemObject(clW);
	//clReleaseMemObject(clWminus);
	//clReleaseMemObject(clU);
	//clReleaseMemObject(clV);
	//clReleaseMemObject(clT);
	//clReleaseMemObject(clRestored);
	//clPadCrop->~clKernel();
	//clWienerW->~clKernel();
	//clWienerWMinus->~clKernel();
	//clWienerV->~clKernel();
	//clWienerT->~clKernel();
	//clWienerU->~clKernel();
	//clWaveTransferFunction->~clKernel();
	//clWaveTransferFunctionMinus->~clKernel();
	
	// Don't need these anymore
	PCPCF->~clKernel();
	clReleaseMemObject(clFFTImage2);
	clReleaseMemObject(clPCPCFResult);
	clReleaseMemObject(clImage2);

	// Now start working on finding the actual astigmatism and defocus from the reconstruction....

	// Make additional kernels
	clKernel* clGetQ				= new clKernel(getQsource, context, cldev, "clCalculateQ", clq);
	clKernel* clSumReduction		= new clKernel(sumReductionsource,context, cldev, "clSumReduction", clq);
	clKernel* clMinusWavefunction	= new clKernel(minuswavefunctionsource, context, cldev, "clMinusWavefunction", clq);
	clKernel* clCalculatePCI		= new clKernel(getPCIsource,context, cldev, "clCalculatePCI", clq);
	clGetQ->BuildKernel();
	clSumReduction->BuildKernel();
	clMinusWavefunction->BuildKernel();
	clCalculatePCI->BuildKernel();

	// New memory required
	cl_mem clRestoredMinus	= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clQ				= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCI			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIC			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIM			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIS			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);

	// Set Kernel Arguments
	clGetQ->SetArgT(0,clW);
	clGetQ->SetArgT(1,clWminus);
	clGetQ->SetArgT(2,clV);
	clGetQ->SetArgT(3,clQ);
	clGetQ->SetArgT(4,sizeX);
	clGetQ->SetArgT(5,sizeY);

	clMinusWavefunction->SetArgT(0,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clMinusWavefunction->SetArgT(1,clRestoredMinus);
	clMinusWavefunction->SetArgT(2,sizeX);
	clMinusWavefunction->SetArgT(3,sizeY);

	clGetQ->Enqueue(globalWorkSize);
	clMinusWavefunction->Enqueue(globalWorkSize);

	// Get sum of Q by my awesome reduction kernel ¬_¬ 
	int totalSize = sizeX*sizeY;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	cl_mem clSumOutput = clCreateBuffer(context,CL_MEM_READ_WRITE,nGroups*sizeof(std::complex<float>),0,&status);
	std::vector< std::complex< float > > sumQs( nGroups );

	clSumReduction->SetArgT(0,clQ);
	clSumReduction->SetArgT(1,clSumOutput);
	clSumReduction->SetArgT(2,totalSize);
	clSumReduction->SetArgLocalMemory(3,256,clFloat2);

	clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 0, NULL, NULL );

	
	// Find out which numbers to read back
	float sumQ2 = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sumQ2 += sumQs[i].real();
	}

	// Trial lots of C1 Values
	int trials2 = 50;

	float* Fpci = new float[trials2];
	float* Fpcic = new float[trials2];
	float* Fpcis = new float[trials2];
	float* Fpcim = new float[trials2];
	float* Fpci2 = new float[trials2];

	float A1r = 0.0f;
	float A1i = 0.0f;
	float zero2 = 0.0f;
	
	float StartDf = -75;
	float StepDf = 3;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	clCalculatePCI->SetArgT(0,clQ);
	clCalculatePCI->SetArgT(1,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clCalculatePCI->SetArgT(2,clRestoredMinus);
	clCalculatePCI->SetArgT(3,clxFrequencies);
	clCalculatePCI->SetArgT(4,clyFrequencies);
	clCalculatePCI->SetArgT(5,clPCI);
	clCalculatePCI->SetArgT(6,clPCIC);
	clCalculatePCI->SetArgT(7,clPCIM);
	clCalculatePCI->SetArgT(8,clPCIS);
	clCalculatePCI->SetArgT(9,sizeX);
	clCalculatePCI->SetArgT(10,sizeY);
	clCalculatePCI->SetArgT(12,zero2);
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);
	clCalculatePCI->SetArgT(15,wavelength);
	clCalculatePCI->SetArgT(16,Abb.kmax);
	

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
 
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);
	
		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.
		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpci[trial] = fpcitrial/sumQ2;
		Fpcic[trial] = fpcictrial/sumQ2;
		Fpcis[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcim[trial] = fpcimtrial/sumQ2;
	}

	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = StartDf + trial*StepDf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial <trials2 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess)
		);


	DM::Image pcigraph = DM::RealImage("PCI",4,trials2);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial <trials2 ; trial++)
	{
		pciData[trial] = Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	DM::UpdateImage(pcigraph);


	delete[] Fpci;
	delete[] Fpci2;
	delete[] Fpcic;
	delete[] Fpcim;
	delete[] Fpcis;

	// First PCI Stage Complete

	// Redo more accurate PCI in correct area 
	// Trial lots of C1 Values

	int trials3 = 50;
	
	float* Fpciref = new float[trials3];
	float* Fpcicref = new float[trials3];
	float* Fpcisref = new float[trials3];
	float* Fpcimref = new float[trials3];
	float* Fpci2ref = new float[trials3];

	A1r = A1rGuess;
	A1i = A1iGuess;
	
	StartDf = DefocusGuess-10;
	StepDf = 21.0f/trials2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	// Set Kernel Arguments - DONE ALREADY NOT CHANGED
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);


	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);

		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.

		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );
		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpciref[trial] = fpcitrial/sumQ2;
		Fpcicref[trial] = fpcictrial/sumQ2;
		Fpcisref[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcimref[trial] = fpcimtrial/sumQ2;

	}

	float maxheighta = 0;
	float C1aa = 0;
	int besttriala = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		if(Fpcimref[trial] > maxheighta)
		{
			maxheighta = Fpcimref[trial];
			C1aa = StartDf + trial*StepDf;
			besttriala = trial;
		}
	}

	float astigaxisa = atan2(Fpcisref[besttriala],Fpcicref[besttriala])/2;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		Fpci2ref[trial] = Fpciref[trial] - Fpcicref[trial] * cos(2*astigaxisa) - Fpcisref[trial]*sin(2*astigaxisa);
	}

	float maxheight2a = 0;
	float C1ba = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess2),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess2),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess2)
		);

	// Print out some results and graphs

	DM::Image pcigraph2 = DM::RealImage("PCI - Refined",4,trials3);
	Gatan::PlugIn::ImageDataLocker pciLocker2(pcigraph2);
	float* pciData2 = (float*) pciLocker2.get();

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		pciData2[trial]=Fpciref[trial];
	}

	pciLocker2.~ImageDataLocker();
	pcigraph2.GetOrCreateImageDocument().Show();
	pcigraph2.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	

	delete[] Fpciref;
	delete[] Fpcicref;
	delete[] Fpcisref;
	delete[] Fpcimref;
	delete[] Fpci2ref;

	// Second PCI Stage Complete


	// Then do restoration again at correct focal/astig value

	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..

		std::vector<cl_float2> image(xDim*yDim);

		// Want to rotate and scale correct this image -- do rotations from reference defocus not absolute defocus or it will (probably) affect shifts.

		float expecteddifference = (imagelist[i]-referenceimage) * options.focalstep;

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}
	
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();


		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]] + DefocusGuess2;

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);
		clWaveTransferFunction->SetArgT(8,A1rGuess2);
		clWaveTransferFunctionMinus->SetArgT(8,A1rGuess2);
		clWaveTransferFunction->SetArgT(9,A1iGuess2);
		clWaveTransferFunctionMinus->SetArgT(9,A1iGuess2);
		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);

	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Final Restored EW",sizeX,sizeY,clFloat2,clq);

	// New Image added for doing rotation and scale correcting
	copyImage.clear();
	rotscaledseries.clear();

	clReleaseMemObject(clImage1);
	clReleaseMemObject(clFFTImage1);
	clReleaseMemObject(fullImage);
	clReleaseMemObject(clw);
	clReleaseMemObject(clwminus);
	clReleaseMemObject(clW);
	clReleaseMemObject(clWminus);
	clReleaseMemObject(clU);
	clReleaseMemObject(clV);
	clReleaseMemObject(clT);
	clReleaseMemObject(clRestored);
	clReleaseMemObject(clRestoredMinus);
	clReleaseMemObject(clQ);
	clReleaseMemObject(clPCI);
	clReleaseMemObject(clPCIC);
	clReleaseMemObject(clPCIM);
	clReleaseMemObject(clPCIS);
	clReleaseMemObject(clSumOutput);

	clPadCrop->~clKernel();
	clWienerW->~clKernel();
	clWienerWMinus->~clKernel();
	clWienerV->~clKernel();
	clWienerT->~clKernel();
	clWienerU->~clKernel();
	clWaveTransferFunction->~clKernel();
	clWaveTransferFunctionMinus->~clKernel();
	clMakeRestored->~clKernel();

	clGetQ->~clKernel();
	clSumReduction->~clKernel();
	clMinusWavefunction->~clKernel();
	clCalculatePCI->~clKernel();
	
};



// Maybe a third version that just prerotates and scales the entire stack then just runs the original reconstruction on that?
// Sometimes getting funny results after a certain point in code, not sure which step is failing but it starts giving same defocus for every image which should be impossible...
void PCPCFWrapper4(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb)
{
	//Utility::SetResultWindow("Should Align = "+boost::lexical_cast<std::string>(options.determinefocus)+"\n");
	//Utility::SetResultWindow("Steps = "+boost::lexical_cast<std::string>(options.steps)+"\n");

	std::vector<int> imagelist;

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	cl_int status;
	// Construct PCPCF Kernel

	clKernel* PCPCF = new clKernel(pcpcfSource,context,cldev,"clPCPCF",clq);
	PCPCF->BuildKernel();

	// Setup required memory for PCPCF subroutine

	cl_mem clImage1			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clImage2			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage1		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage2		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clPCPCFResult	= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);

	cl_mem fullImage		= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem rotScaleImage	= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);

	cl_mem clW				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clWminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clw				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clwminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clT				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clU				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clV				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clRestored		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);

	int init = 1;
	int noinit = 0;

	PCPCF->SetArgT(0,clFFTImage1);
	PCPCF->SetArgT(1,clFFTImage2);
	PCPCF->SetArgT(2,clPCPCFResult);
	PCPCF->SetArgT(3,clxFrequencies);
	PCPCF->SetArgT(4,clyFrequencies);
	PCPCF->SetArgT(5,sizeX);
	PCPCF->SetArgT(6,sizeY);
	PCPCF->SetArgT(8,wavelength);
	PCPCF->SetArgT(9,options.pcpcfkmax);

	/* Make all the other kernels */
	clKernel* clPadCrop = new clKernel(PadCropSource,context,cldev,"clPadCrop",clq);
	clPadCrop->BuildKernel();

	clPadCrop->SetArgT(0,fullImage);
	clPadCrop->SetArgT(1,clImage1);
	clPadCrop->SetArgT(2,xDim); // long not int
	clPadCrop->SetArgT(3,yDim); // long not int
	clPadCrop->SetArgT(10,sizeX);
	clPadCrop->SetArgT(11,sizeY);
	clPadCrop->SetArgT(12,iTop);
	clPadCrop->SetArgT(13,iLeft);

	clKernel* clWaveTransferFunction = new clKernel(wavetransferfunctionsource,context,cldev,"clWaveTransferFunction",clq);
	clWaveTransferFunction->BuildKernel();

	clWaveTransferFunction->SetArgT(0,clw);
	clWaveTransferFunction->SetArgT(1,clxFrequencies);
	clWaveTransferFunction->SetArgT(2,clyFrequencies);
	clWaveTransferFunction->SetArgT(3,sizeX);
	clWaveTransferFunction->SetArgT(4,sizeY);
	clWaveTransferFunction->SetArgT(5,wavelength);
	clWaveTransferFunction->SetArgT(6,Abb.beta);
	clWaveTransferFunction->SetArgT(7,Abb.delta);
	clWaveTransferFunction->SetArgT(8,Abb.A1r);
	clWaveTransferFunction->SetArgT(9,Abb.A1i);
	clWaveTransferFunction->SetArgT(11,Abb.Cs);
	clWaveTransferFunction->SetArgT(12,Abb.kmax);

	clKernel* clWaveTransferFunctionMinus = new clKernel(wavetransferfunctionminussource,context,cldev,"clWaveTransferFunctionMinus",clq);
	clWaveTransferFunctionMinus->BuildKernel();

	clWaveTransferFunctionMinus->SetArgT(0,clwminus);
	clWaveTransferFunctionMinus->SetArgT(1,clxFrequencies);
	clWaveTransferFunctionMinus->SetArgT(2,clyFrequencies);
	clWaveTransferFunctionMinus->SetArgT(3,sizeX);
	clWaveTransferFunctionMinus->SetArgT(4,sizeY);
	clWaveTransferFunctionMinus->SetArgT(5,wavelength);
	clWaveTransferFunctionMinus->SetArgT(6,Abb.beta);
	clWaveTransferFunctionMinus->SetArgT(7,Abb.delta);
	clWaveTransferFunctionMinus->SetArgT(8,Abb.A1r);
	clWaveTransferFunctionMinus->SetArgT(9,Abb.A1i);
	clWaveTransferFunctionMinus->SetArgT(11,Abb.Cs);
	clWaveTransferFunctionMinus->SetArgT(12,Abb.kmax);

	clKernel* clWienerW = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerW->BuildKernel();

	clWienerW->SetArgT(0,clW);
	clWienerW->SetArgT(1,clw);
	clWienerW->SetArgT(2,sizeX);
	clWienerW->SetArgT(3,sizeY);
	clWienerW->SetArgT(4,init);

	clKernel* clWienerWMinus = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerWMinus->BuildKernel();

	clWienerWMinus->SetArgT(0,clWminus);
	clWienerWMinus->SetArgT(1,clwminus);
	clWienerWMinus->SetArgT(2,sizeX);
	clWienerWMinus->SetArgT(3,sizeY);
	clWienerWMinus->SetArgT(4,init);

	clKernel* clWienerV = new clKernel(wienervsource,context,cldev,"clWienerV",clq);
	clWienerV->BuildKernel();

	clWienerV->SetArgT(0,clV);
	clWienerV->SetArgT(1,clw);
	clWienerV->SetArgT(2,clwminus);
	clWienerV->SetArgT(3,sizeX);
	clWienerV->SetArgT(4,sizeY);
	clWienerV->SetArgT(5,init);

	clKernel* clWienerT = new clKernel(wienertsource,context,cldev,"clWienerT",clq);
	clWienerT->BuildKernel();

	clWienerT->SetArgT(0,clT);
	clWienerT->SetArgT(1,clw);
	clWienerT->SetArgT(2,clFFTImage1);
	clWienerT->SetArgT(3,sizeX);
	clWienerT->SetArgT(4,sizeY);
	clWienerT->SetArgT(5,init);

	clKernel* clWienerU = new clKernel(wienerusource,context,cldev,"clWienerU",clq);
	clWienerU->BuildKernel();

	clWienerU->SetArgT(0,clU);
	clWienerU->SetArgT(1,clwminus);
	clWienerU->SetArgT(2,clFFTImage1);
	clWienerU->SetArgT(3,sizeX);
	clWienerU->SetArgT(4,sizeY);
	clWienerU->SetArgT(5,init);
	
	clKernel* clMakeRestored = new clKernel(makerestoredsource,context,cldev,"clMakeRestored",clq);
	clMakeRestored->BuildKernel();

	clMakeRestored->SetArgT(0,clW);
	clMakeRestored->SetArgT(1,clWminus);
	clMakeRestored->SetArgT(2,clV);
	clMakeRestored->SetArgT(3,clT);
	clMakeRestored->SetArgT(4,clU);
	clMakeRestored->SetArgT(5,clFFTImage1);
	clMakeRestored->SetArgT(6,sizeX);
	clMakeRestored->SetArgT(7,sizeY);

	clKernel* clRotScale = new clKernel(RotScaleSource,context,cldev,"clRotScale",clq);
	clRotScale->BuildKernel();

	if(options.rotscale)
	{
		clRotScale->SetArgT(0,fullImage);
		clRotScale->SetArgT(1,rotScaleImage);
		clRotScale->SetArgT(2,xDim);
		clRotScale->SetArgT(3,yDim);
		clRotScale->SetArgT(4,options.magcal);
		clRotScale->SetArgT(5,options.rotcal);
	}

	size_t* fullWorkSize = new size_t[3];
	
	// There are 'elements' work items
	fullWorkSize[0] = xDim;
	fullWorkSize[1] = yDim;
	fullWorkSize[2] = 1;

	std::vector<std::complex<float>> copyImage(xDim*yDim) ;
	std::vector<float> rotscaledseries(xDim*yDim*numberOfImages);

	// Make a new rotation aligned stack based on the expected defocus step.
	// To prevent edges going missing should start at end that should be smallest so all images are only scaled upwards....


	
	// Rotate and Scale this image by expected defocus..
	if(options.rotscale)
	{
		int startingimage = 0;

		if(options.magcal > 1)
		{
			startingimage = 0;
		}
		else
		{
			startingimage = numberOfImages -1;
		}


		for(int im = 0; im < numberOfImages ; im++ )
		{
			float expecteddifference = (im-startingimage)*options.focalstep;
			clRotScale->SetArgT(6,expecteddifference);

			// Copy image into fullimage
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					copyImage[i+j*xDim] = seriesdata[im*xDim*yDim + i + (j)*xDim];
				}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					rotscaledseries[im*xDim*yDim +i + j*xDim] = copyImage[i + (j)*xDim].real();
				}
		}		
	} 
	else
	{
		// Fill new series with old series data.
		for(int im = 0; im < numberOfImages ; im++ )
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					rotscaledseries[im*xDim*yDim +i + j*xDim] = seriesdata[im*xDim*yDim +i + j*xDim];
				}
	}

	// Define and index space of work items for execution
	// A workgroup size is not required but can be used
	size_t* globalWorkSize = new size_t[3];
	
	// There are 'elements' work items
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	globalWorkSize[2] = 1;

	std::vector< std::complex< float > > dataOne( sizeX*sizeY );
	std::vector< std::complex< float > > dataTwo( sizeX*sizeY );

	int referenceimage = options.reference;
	int gonext = 1;
	int nextdir = -1;
	
	// First three have to be done via standard methods.

	// Initialise image numbers
	int imageone = referenceimage;
	int imagetwo = referenceimage;
	int currentimage;

	imagelist.push_back(referenceimage);

	// loop over all of the images
	for(int n = 0; n < numberOfImages - 1; n++)
	{
		Debug("Start Image "+boost::lexical_cast<std::string>(n));

		 // For RECONPCPCF part

		// First 2 registrations cannot be done against reconstruction
		// have to be done against reference image...
		if(n < 2)
		{
			imageone = referenceimage;
			imagetwo += (gonext*nextdir);

			while(! (imagetwo < numberOfImages && imagetwo >=0)) // not an actual image 
			{
				gonext++;
				nextdir *=-1;
				imagetwo += (gonext*nextdir);
			}

			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = (imagetwo-imageone) * options.focalstep;

			for(int j = 0; j < sizeY;j++)
				for(int i = 0; i < sizeX;i++)
				{
					dataOne[i+j*sizeX] = rotscaledseries[imageone*xDim*yDim + i + iLeft + (j+iTop)*xDim];
					dataTwo[i+j*sizeX] = rotscaledseries[imagetwo*xDim*yDim + i + iLeft + (j+iTop)*xDim];
				}

	
			clEnqueueWriteBuffer( clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
						0, NULL, NULL );

			FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			for(int trial = 0; trial < numberoftrials; trial++)
			{

				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected + searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep 
							+ ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);

				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(imagetwo)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY,options.maxdrift);

				trialxshifts[trial] = xShift;
				trialyshifts[trial] = yShift;
				trialsubxshifts[trial] = xShift + subXShift;
				trialsubyshifts[trial] = yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;
			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[imagetwo] = trialxshifts[trial];
					yShiftVals[imagetwo] = trialyshifts[trial];
					subXShifts[imagetwo] = trialsubxshifts[trial];
					subYShifts[imagetwo] = trialsubyshifts[trial];
					defocusshifts[imagetwo] = trialdefocus[trial];
					Debug("Set Shifts for "+boost::lexical_cast<std::string>(imagetwo));
				}
			}

			imagelist.push_back(imagetwo);
		}

		if(n==2)
		{
			// To make sure it is initialised from end of first section
			currentimage = imagetwo;
		}
	
		if(n>=2)
		{
			// Now we can perform reconstruction and register to the reconstruction instead of standard pcpcf
			currentimage += (gonext*nextdir);

			while(! (currentimage < numberOfImages && currentimage >=0)) // not an actual image
			{
				gonext++;
				nextdir *=-1;
				currentimage += (gonext*nextdir);
			}

			// Need to reconstruct the previous images based on there focus and shifts relative to reference image...

			
			// Need to extract correct area from each image to perform reconstruction
			// 2 Possibilites - If an area of the correct size is not available then give an error...
			// or pad the area with zeroes for ALL images - need to first determine the extents of available area before we extract from any of the images.

			// Find maximum negative and positive shifts
			// All shifts have to be stored relative to the reference image.
			float maxnegx = 0;
			float maxnegy = 0;
			float maxposx = 0;
			float maxposy = 0;

			// Could just loop over this list of currently registered images instead...
			for(int i = 1; i <= imagelist.size(); i++)
			{
				if(subXShifts[imagelist[i-1]] < maxnegx)
					maxnegx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] < maxnegy)
					maxnegy = subYShifts[imagelist[i-1]];
				if(subXShifts[imagelist[i-1]] > maxposx)
					maxposx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] > maxposy)
					maxposy = subYShifts[imagelist[i-1]];
			}

			// Determine amount to pad on either direction
			int padTop = 0;
			int padLeft = 0;
			int padRight = 0;
			int padBottom = 0;

			if (abs(maxnegx)-iLeft > 0)
				padLeft = ceil(abs(maxnegx)-iLeft);

			if (maxposx-iRight > 0)
				padRight = ceil(maxposx-iRight);

			if (abs(maxnegy)-iTop > 0)
				padTop = ceil(abs(maxnegy)-iTop);

			if (maxposy-iBottom > 0)
				padBottom = ceil(maxposy-iBottom);	

			for(int i = 0 ; i < imagelist.size() ; i++)
			{
				if(i == 0) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,init);
					clWienerWMinus->SetArgT(4,init);
					clWienerV->SetArgT(5,init);
					clWienerT->SetArgT(5,init);
					clWienerU->SetArgT(5,init);
				}
				if(i == 1) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,noinit);
					clWienerWMinus->SetArgT(4,noinit);
					clWienerV->SetArgT(5,noinit);
					clWienerT->SetArgT(5,noinit);
					clWienerU->SetArgT(5,noinit);
				}
				
				// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
				clPadCrop->SetArgT(4,padLeft);
				clPadCrop->SetArgT(5,padRight);
				clPadCrop->SetArgT(6,padTop);
				clPadCrop->SetArgT(7,padBottom);
				clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
				clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
				// Copy correct image into full image..
				std::vector<cl_float2> image(xDim*yDim);

				// Now add image to reconstruction - Make WTF - WTFminus			
				float defocus = defocusshifts[imagelist[i]];

				for(int j = 0 ; j < xDim*yDim ; j++)
				{
					image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
					image[j].s[1] = 0;
				}


				clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
				clPadCrop->Enqueue(globalWorkSize);

				image.clear();

				// Get FFT of this image
				FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

				clWaveTransferFunction->SetArgT(10,defocus);
				clWaveTransferFunctionMinus->SetArgT(10,defocus);

				clWaveTransferFunction->Enqueue(globalWorkSize);
				clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

				// add to W,W-,V,T,U 
				clWienerW->Enqueue(globalWorkSize);
				clWienerWMinus->Enqueue(globalWorkSize);
				clWienerV->Enqueue(globalWorkSize);
				clWienerT->Enqueue(globalWorkSize);
				clWienerU->Enqueue(globalWorkSize);
			}

			// make restored
			clMakeRestored->Enqueue(globalWorkSize);
			 
			Debug("Made Restored");

			// inverse fft for display
			FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);
			
			//Show Each reconstruction
			//Utility::PrintCLMemToImage(clRestored,"Restored",sizeX,sizeY,clFloat2,clq);
			
			// PCF with the reconstruction
			// Adjust reconstruction or do PCPCF?
			// TODO: at the moment reconstruction is in reference plane and images registered by PCPCF with this plane....
			// Also do i need to get modulus and FFT or not? - seems to work OK
			
			// Work out what focus difference is expected between these 2 images.
			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone


			// Here take an image thats already shifted approximately to the right place... then add this back onto calculated shifts...
			// Need to check image above or below to see if its registered yet

			bool prevImageRegistered = false;
			bool nextImageRegistered = false;

			for(int i = 1; i <= imagelist.size(); i++)
			{
				Debug("imglist: "+boost::lexical_cast<std::string>(imagelist[i-1]));
								
				if(imagelist[i-1]==currentimage-1)
				{
					Debug("PREV");
					prevImageRegistered = true;
					break;
				}
				if(imagelist[i-1]==currentimage+1)
				{
					Debug("NEXT");
					nextImageRegistered = true;
					break;
				}
			}

			int preshiftx = 0;
			int preshifty = 0;
			
			
			if(prevImageRegistered)
			{
				preshiftx = xShiftVals[currentimage-1];
				preshifty = yShiftVals[currentimage-1];
			}

			if(nextImageRegistered)
			{
				preshiftx = xShiftVals[currentimage+1];
				preshifty = yShiftVals[currentimage+1];
				Debug("xshift: "+boost::lexical_cast<std::string>(preshiftx));
				Debug("yshift: "+boost::lexical_cast<std::string>(preshifty));
			}
			

			Debug("Calculated Preshifts");

			// Cant have it checking an image area that doesnt exist anymore
			if( iLeft + preshiftx < 0)
				preshiftx -=( preshiftx + iLeft );

			if( iLeft + preshiftx + sizeX >= xDim)
				preshiftx -=( iLeft + preshiftx + sizeX - xDim);

			if( iTop + preshifty < 0)
				preshifty -=( preshifty + iTop );

			if( iTop + preshifty + sizeY >= yDim)
				preshiftx -=( iTop + preshifty + sizeY - yDim);


			for(int j = 0; j < sizeY; j++)
				for(int i = 0; i < sizeX; i++)
				{
					dataTwo[i+j*sizeX] = rotscaledseries[currentimage*xDim*yDim + i + iLeft + preshiftx + (j+iTop + preshifty)*xDim];
				}

			Debug("Made Preshifted image");

			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 0, NULL, NULL );

			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			float expecteddifference = (currentimage-referenceimage) * options.focalstep;

			// Get number of trial steps
			// Assume one if not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			int* trialxshifts = new int[numberoftrials];
			int* trialyshifts = new int[numberoftrials];
			float* trialsubxshifts = new float[numberoftrials];
			float* trialsubyshifts = new float[numberoftrials];
			float* trialdefocus = new float[numberoftrials];
			float* peakheights = new float[numberoftrials];

			Debug("Starting Trials");

			for(int trial = 0; trial < numberoftrials; trial++)
			{
				Utility::SetResultWindow("Trial "+boost::lexical_cast<std::string>(trial+1)+" of "+boost::lexical_cast<std::string>(numberoftrials)+" \n");

				// go from expected - searchpercentage to expected = searchpercentage
				float trialdef = expecteddifference - ((float)options.searchpercentage/100.0f)*options.focalstep 
							+ ((float)trial/((float)numberoftrials-1.0f))*2.0f*((float)options.searchpercentage/100.0f)*options.focalstep; 

				// above formula wont work when number of trials is 1
				if(numberoftrials == 1)
				{
					trialdef = expecteddifference;
				}

				PCPCF->SetArgT(7,trialdef);
				PCPCF->Enqueue(globalWorkSize);

				// Now Inverse FFT
				// Using memory reserved for storing images as it is not needed at this point.
				FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);
				clEnqueueReadBuffer( clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
							0, NULL, NULL );

				// Could consider having this displayed to see roughly how well it is working...
				//Utility::PrintCLMemToImage(clFFTImage1,"FFT1",sizeX,sizeY,clFloat2,clq);
				//Utility::PrintCLMemToImage(clImage1,"PCPCF",sizeX,sizeY,clFloat2,clq);
				Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(currentimage)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

				// Find shift from max height position
				int xShift;
				int yShift;
				float subXShift;
				float subYShift;
				float maxHeight1;

				// Translate linear array index into row and column.
				PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY,options.maxdrift);

				trialxshifts[trial] = preshiftx + xShift;
				trialyshifts[trial] = preshifty + yShift;
				trialsubxshifts[trial] = preshiftx + xShift + subXShift;
				trialsubyshifts[trial] = preshifty + yShift + subYShift;
				trialdefocus[trial] = trialdef;
				peakheights[trial] = maxHeight1;

				Debug("Set Peak Height: "+boost::lexical_cast<std::string>(maxHeight1));

			}

			// Remember the reference image has shifts of 0,0
			// Shifts now expressed in terms of difference to reference image.

			// Get best value and store in shifts.
			float bestheight = 0.0f;
			for(int trial = 0; trial < numberoftrials; trial++)
			{
				if(peakheights[trial] > bestheight)
				{
					bestheight = peakheights[trial];
					xShiftVals[currentimage] = trialxshifts[trial];
					yShiftVals[currentimage] = trialyshifts[trial];
					subXShifts[currentimage] = trialsubxshifts[trial];
					subYShifts[currentimage] = trialsubyshifts[trial];
					defocusshifts[currentimage] = trialdefocus[trial];

					Debug("Set Shifts for "+boost::lexical_cast<std::string>(currentimage));
				}
			}

			imagelist.push_back(currentimage);

			delete[] trialxshifts;
			delete[] trialyshifts;
			delete[] trialsubxshifts;
			delete[] trialsubyshifts;
			delete[] trialdefocus;
			delete[] peakheights;
		}
		
		// to get to next image..
		gonext++;
		nextdir *=-1;
		// image = previous +(gonext*nextdir);
		// checkbounds then if not available do again.
		
	}

	Debug("Starting Final Reconstruction");


	/* Now rebuild reconstruction including the final registered image */
	// All shifts have to be stored relative to the reference image.
	float maxnegx = 0;
	float maxnegy = 0;
	float maxposx = 0;
	float maxposy = 0;

	// Could just loop over this list of currently registered images instead...
	for(int i = 1; i <= numberOfImages; i++)
	{
		if(subXShifts[i-1] < maxnegx)
			maxnegx = subXShifts[i-1];
		if(subYShifts[i-1] < maxnegy)
			maxnegy = subYShifts[i-1];
		if(subXShifts[i-1] > maxposx)
			maxposx = subXShifts[i-1];
		if(subYShifts[i-1] > maxposy)
			maxposy = subYShifts[i-1];
	}

	// Determine amount to pad on either direction
	int padTop = 0;
	int padLeft = 0;
	int padRight = 0;
	int padBottom = 0;

	if (abs(maxnegx)-iLeft > 0)
		padLeft = ceil(abs(maxnegx)-iLeft);

	if (maxposx-iRight > 0)
		padRight = ceil(maxposx-iRight);

	if (abs(maxnegy)-iTop > 0)
		padTop = ceil(abs(maxnegy)-iTop);

	if (maxposy-iBottom > 0)
		padBottom = ceil(maxposy-iBottom);


	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]];

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();

		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);

		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);
			 
	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Restored EW",sizeX,sizeY,clFloat2,clq);

	// Now all should be registered - produce the final restored image and any other relevant data like Shift Graphs, and Drift Corrected Stack
	Gatan::DM::Image driftcorrected = DM::RealImage("Drift Corrected Stack",4,sizeX,sizeY,numberOfImages);
	Gatan::PlugIn::ImageDataLocker driftlocker = Gatan::PlugIn::ImageDataLocker(driftcorrected);
	float* driftstack = (float*) driftlocker.get();

	std::vector<cl_float2> cropimage(sizeX*sizeY);

	for(int i = 0 ; i < numberOfImages ; i++)
	{
			
		int padTop = 0;
		int padLeft = 0;
		int padRight = 0;
		int padBottom = 0;

		if (abs(maxnegx)-iLeft > 0)
			padLeft = ceil(abs(maxnegx)-iLeft);

		if (maxposx-iRight > 0)
			padRight = ceil(maxposx-iRight);

		if (abs(maxnegy)-iTop > 0)
			padTop = ceil(abs(maxnegy)-iTop);

		if (maxposy-iBottom > 0)
			padBottom = ceil(maxposy-iBottom);

		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[i]);
		clPadCrop->SetArgT(9,subYShifts[i]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[i];
			
		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[i*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		// Copy full image to GPU and run crop kernel
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		// Copy image data back to host
		clEnqueueReadBuffer(clq->cmdQueue,clImage1,CL_TRUE,0,sizeX*sizeY*sizeof(cl_float2),&cropimage[0],0,NULL,NULL);
		
		// Copy to drift corrected stack
		for(int k = 0 ; k < sizeX*sizeY ; k++)
		{
			driftstack[i*sizeX*sizeY + k] = cropimage[k].s[0];
		}

		image.clear();
	}

	cropimage.clear();
	driftlocker.~ImageDataLocker();

	// Display drift corrected
	driftcorrected.GetOrCreateImageDocument().Show();

	dataOne.clear();
	dataTwo.clear();



	// Lots of these have been moved later on because I added extra steps that require them.
	
	//clReleaseMemObject(clImage1);
	//clReleaseMemObject(clFFTImage1);
	//clReleaseMemObject(fullImage);
	//clReleaseMemObject(clw);
	//clReleaseMemObject(clwminus);
	//clReleaseMemObject(clW);
	//clReleaseMemObject(clWminus);
	//clReleaseMemObject(clU);
	//clReleaseMemObject(clV);
	//clReleaseMemObject(clT);
	//clReleaseMemObject(clRestored);
	//clPadCrop->~clKernel();
	//clWienerW->~clKernel();
	//clWienerWMinus->~clKernel();
	//clWienerV->~clKernel();
	//clWienerT->~clKernel();
	//clWienerU->~clKernel();
	//clWaveTransferFunction->~clKernel();
	//clWaveTransferFunctionMinus->~clKernel();
	
	// Don't need these anymore
	PCPCF->~clKernel();
	clReleaseMemObject(clFFTImage2);
	clReleaseMemObject(clPCPCFResult);
	clReleaseMemObject(clImage2);

	// Now start working on finding the actual astigmatism and defocus from the reconstruction....

	// Make additional kernels
	clKernel* clGetQ				= new clKernel(getQsource, context, cldev, "clCalculateQ", clq);
	clKernel* clSumReduction		= new clKernel(sumReductionsource,context, cldev, "clSumReduction", clq);
	clKernel* clMinusWavefunction	= new clKernel(minuswavefunctionsource, context, cldev, "clMinusWavefunction", clq);
	clKernel* clCalculatePCI		= new clKernel(getPCIsource,context, cldev, "clCalculatePCI", clq);
	clGetQ->BuildKernel();
	clSumReduction->BuildKernel();
	clMinusWavefunction->BuildKernel();
	clCalculatePCI->BuildKernel();

	// New memory required
	cl_mem clRestoredMinus	= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clQ				= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCI			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIC			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIM			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIS			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);

	// Set Kernel Arguments
	clGetQ->SetArgT(0,clW);
	clGetQ->SetArgT(1,clWminus);
	clGetQ->SetArgT(2,clV);
	clGetQ->SetArgT(3,clQ);
	clGetQ->SetArgT(4,sizeX);
	clGetQ->SetArgT(5,sizeY);

	clMinusWavefunction->SetArgT(0,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clMinusWavefunction->SetArgT(1,clRestoredMinus);
	clMinusWavefunction->SetArgT(2,sizeX);
	clMinusWavefunction->SetArgT(3,sizeY);

	clGetQ->Enqueue(globalWorkSize);
	clMinusWavefunction->Enqueue(globalWorkSize);

	// Get sum of Q by my awesome reduction kernel ¬_¬ 
	int totalSize = sizeX*sizeY;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	cl_mem clSumOutput = clCreateBuffer(context,CL_MEM_READ_WRITE,nGroups*sizeof(std::complex<float>),0,&status);
	std::vector< std::complex< float > > sumQs( nGroups );

	clSumReduction->SetArgT(0,clQ);
	clSumReduction->SetArgT(1,clSumOutput);
	clSumReduction->SetArgT(2,totalSize);
	clSumReduction->SetArgLocalMemory(3,256,clFloat2);

	clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 0, NULL, NULL );

	
	// Find out which numbers to read back
	float sumQ2 = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sumQ2 += sumQs[i].real();
	}

	// Trial lots of C1 Values
	int trials2 = 50;

	float* Fpci = new float[trials2];
	float* Fpcic = new float[trials2];
	float* Fpcis = new float[trials2];
	float* Fpcim = new float[trials2];
	float* Fpci2 = new float[trials2];

	float A1r = 0.0f;
	float A1i = 0.0f;
	float zero2 = 0.0f;
	
	float StartDf = -50;
	float StepDf = 2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	clCalculatePCI->SetArgT(0,clQ);
	clCalculatePCI->SetArgT(1,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clCalculatePCI->SetArgT(2,clRestoredMinus);
	clCalculatePCI->SetArgT(3,clxFrequencies);
	clCalculatePCI->SetArgT(4,clyFrequencies);
	clCalculatePCI->SetArgT(5,clPCI);
	clCalculatePCI->SetArgT(6,clPCIC);
	clCalculatePCI->SetArgT(7,clPCIM);
	clCalculatePCI->SetArgT(8,clPCIS);
	clCalculatePCI->SetArgT(9,sizeX);
	clCalculatePCI->SetArgT(10,sizeY);
	clCalculatePCI->SetArgT(12,zero2);
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);
	clCalculatePCI->SetArgT(15,wavelength);
	clCalculatePCI->SetArgT(16,Abb.kmax);
	

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
 
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);
	
		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.
		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpci[trial] = fpcitrial/sumQ2;
		Fpcic[trial] = fpcictrial/sumQ2;
		Fpcis[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcim[trial] = fpcimtrial/sumQ2;
	}

	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = StartDf + trial*StepDf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial <trials2 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess)
		);


	DigitalMicrograph::Image pcigraph = DigitalMicrograph::RealImage("PCI",4,trials2);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial <trials2 ; trial++)
	{
		pciData[trial] = Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	DigitalMicrograph::UpdateImage(pcigraph);


	delete[] Fpci;
	delete[] Fpci2;
	delete[] Fpcic;
	delete[] Fpcim;
	delete[] Fpcis;

	// First PCI Stage Complete

	// Redo more accurate PCI in correct area 
	// Trial lots of C1 Values

	int trials3 = 50;
	
	float* Fpciref = new float[trials3];
	float* Fpcicref = new float[trials3];
	float* Fpcisref = new float[trials3];
	float* Fpcimref = new float[trials3];
	float* Fpci2ref = new float[trials3];

	A1r = A1rGuess;
	A1i = A1iGuess;
	
	StartDf = DefocusGuess-10;
	StepDf = 21.0f/trials2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	// Set Kernel Arguments - DONE ALREADY NOT CHANGED
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);


	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);

		//Gatan::PlugIn::ImageDataLocker pciLocker2(pcidisplay);
		//float* pciData2 = (float*) pciLocker2.get();

		// DO PCI Function and get sums for pcic,s,m etc here.

		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );
		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpciref[trial] = fpcitrial/sumQ2;
		Fpcicref[trial] = fpcictrial/sumQ2;
		Fpcisref[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcimref[trial] = fpcimtrial/sumQ2;

	}

	float maxheighta = 0;
	float C1aa = 0;
	int besttriala = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		if(Fpcimref[trial] > maxheighta)
		{
			maxheighta = Fpcimref[trial];
			C1aa = StartDf + trial*StepDf;
			besttriala = trial;
		}
	}

	float astigaxisa = atan2(Fpcisref[besttriala],Fpcicref[besttriala])/2;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		Fpci2ref[trial] = Fpciref[trial] - Fpcicref[trial] * cos(2*astigaxisa) - Fpcisref[trial]*sin(2*astigaxisa);
	}

	float maxheight2a = 0;
	float C1ba = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess2),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess2),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess2)
		);

	// Print out some results and graphs

	DigitalMicrograph::Image pcigraph2 = DigitalMicrograph::RealImage("PCI - Refined",4,trials3);
	Gatan::PlugIn::ImageDataLocker pciLocker2(pcigraph2);
	float* pciData2 = (float*) pciLocker2.get();

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		pciData2[trial]=Fpciref[trial];
	}

	pciLocker2.~ImageDataLocker();
	pcigraph2.GetOrCreateImageDocument().Show();
	pcigraph2.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	

	delete[] Fpciref;
	delete[] Fpcicref;
	delete[] Fpcisref;
	delete[] Fpcimref;
	delete[] Fpci2ref;

	// Second PCI Stage Complete


	// Then do restoration again at correct focal/astig value

	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..

		std::vector<cl_float2> image(xDim*yDim);

		// Want to rotate and scale correct this image -- do rotations from reference defocus not absolute defocus or it will (probably) affect shifts.

		float expecteddifference = (imagelist[i]-referenceimage) * options.focalstep;

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}
	
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();


		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]] + DefocusGuess2;

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);
		clWaveTransferFunction->SetArgT(8,A1rGuess2);
		clWaveTransferFunctionMinus->SetArgT(8,A1rGuess2);
		clWaveTransferFunction->SetArgT(9,A1iGuess2);
		clWaveTransferFunctionMinus->SetArgT(9,A1iGuess2);
		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);

	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Final Restored EW",sizeX,sizeY,clFloat2,clq);

	// New Image added for doing rotation and scale correcting
	copyImage.clear();
	rotscaledseries.clear();

	clReleaseMemObject(clImage1);
	clReleaseMemObject(clFFTImage1);
	clReleaseMemObject(fullImage);
	clReleaseMemObject(clw);
	clReleaseMemObject(clwminus);
	clReleaseMemObject(clW);
	clReleaseMemObject(clWminus);
	clReleaseMemObject(clU);
	clReleaseMemObject(clV);
	clReleaseMemObject(clT);
	clReleaseMemObject(clRestored);
	clReleaseMemObject(clRestoredMinus);
	clReleaseMemObject(clQ);
	clReleaseMemObject(clPCI);
	clReleaseMemObject(clPCIC);
	clReleaseMemObject(clPCIM);
	clReleaseMemObject(clPCIS);
	clReleaseMemObject(clSumOutput);

	clPadCrop->~clKernel();
	clWienerW->~clKernel();
	clWienerWMinus->~clKernel();
	clWienerV->~clKernel();
	clWienerT->~clKernel();
	clWienerU->~clKernel();
	clWaveTransferFunction->~clKernel();
	clWaveTransferFunctionMinus->~clKernel();
	clMakeRestored->~clKernel();

	clGetQ->~clKernel();
	clSumReduction->~clKernel();
	clMinusWavefunction->~clKernel();
	clCalculatePCI->~clKernel();
	
};





/*
void PCPCFWrapperReWrite(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb)
{
	//Utility::SetResultWindow("Should Align = "+boost::lexical_cast<std::string>(options.determinefocus)+"\n");
	//Utility::SetResultWindow("Steps = "+boost::lexical_cast<std::string>(options.steps)+"\n");

	std::vector<int> imagelist;

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	cl_int status;
	// Construct PCPCF Kernel

	clKernel* PCPCF = new clKernel(pcpcfSource,context,cldev,"clPCPCF",clq);
	PCPCF->BuildKernel();


	// Define and index space of work items for execution
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = sizeX;
	globalWorkSize[1] = sizeY;
	globalWorkSize[2] = 1;

	// Setup PCPCF Library

	PCPCFLib PCPCFLibrary;
	PCPCFLibrary.FFT = FFT;
	PCPCFLibrary.PCPCF = PCPCF;
	PCPCFLibrary.globalWorkSize = globalWorkSize;

	// Setup required memory for PCPCF subroutine

	cl_mem clImage1			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clImage2			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage1		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clFFTImage2		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);
	cl_mem clPCPCFResult	= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX * sizeY * sizeof(cl_float2), 0, &status);

	cl_mem fullImage		= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem rotScaleImage	= clCreateBuffer( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);

	cl_mem clW				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clWminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( float ), 0, &status);
	cl_mem clw				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clwminus			= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clT				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clU				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clV				= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);
	cl_mem clRestored		= clCreateBuffer ( context, CL_MEM_READ_WRITE, sizeX *sizeY * sizeof( std::complex< float > ), 0, &status);

	int init = 1;
	int noinit = 0;

	PCPCF->SetArgT(0,clFFTImage1);
	PCPCF->SetArgT(1,clFFTImage2);
	PCPCF->SetArgT(2,clPCPCFResult);
	PCPCF->SetArgT(3,clxFrequencies);
	PCPCF->SetArgT(4,clyFrequencies);
	PCPCF->SetArgT(5,sizeX);
	PCPCF->SetArgT(6,sizeY);
	PCPCF->SetArgT(8,wavelength);
	PCPCF->SetArgT(9,options.pcpcfkmax);

	// Make all the other kernels 
	clKernel* clPadCrop = new clKernel(PadCropSource,context,cldev,"clPadCrop",clq);
	clPadCrop->BuildKernel();

	clPadCrop->SetArgT(0,fullImage);
	clPadCrop->SetArgT(1,clImage1);
	clPadCrop->SetArgT(2,xDim); // long not int
	clPadCrop->SetArgT(3,yDim); // long not int
	clPadCrop->SetArgT(10,sizeX);
	clPadCrop->SetArgT(11,sizeY);
	clPadCrop->SetArgT(12,iTop);
	clPadCrop->SetArgT(13,iLeft);

	clKernel* clWaveTransferFunction = new clKernel(wavetransferfunctionsource,context,cldev,"clWaveTransferFunction",clq);
	clWaveTransferFunction->BuildKernel();

	clWaveTransferFunction->SetArgT(0,clw);
	clWaveTransferFunction->SetArgT(1,clxFrequencies);
	clWaveTransferFunction->SetArgT(2,clyFrequencies);
	clWaveTransferFunction->SetArgT(3,sizeX);
	clWaveTransferFunction->SetArgT(4,sizeY);
	clWaveTransferFunction->SetArgT(5,wavelength);
	clWaveTransferFunction->SetArgT(6,Abb.beta);
	clWaveTransferFunction->SetArgT(7,Abb.delta);
	clWaveTransferFunction->SetArgT(8,Abb.A1r);
	clWaveTransferFunction->SetArgT(9,Abb.A1i);
	clWaveTransferFunction->SetArgT(11,Abb.Cs);
	clWaveTransferFunction->SetArgT(12,Abb.kmax);

	clKernel* clWaveTransferFunctionMinus = new clKernel(wavetransferfunctionminussource,context,cldev,"clWaveTransferFunctionMinus",clq);
	clWaveTransferFunctionMinus->BuildKernel();

	clWaveTransferFunctionMinus->SetArgT(0,clwminus);
	clWaveTransferFunctionMinus->SetArgT(1,clxFrequencies);
	clWaveTransferFunctionMinus->SetArgT(2,clyFrequencies);
	clWaveTransferFunctionMinus->SetArgT(3,sizeX);
	clWaveTransferFunctionMinus->SetArgT(4,sizeY);
	clWaveTransferFunctionMinus->SetArgT(5,wavelength);
	clWaveTransferFunctionMinus->SetArgT(6,Abb.beta);
	clWaveTransferFunctionMinus->SetArgT(7,Abb.delta);
	clWaveTransferFunctionMinus->SetArgT(8,Abb.A1r);
	clWaveTransferFunctionMinus->SetArgT(9,Abb.A1i);
	clWaveTransferFunctionMinus->SetArgT(11,Abb.Cs);
	clWaveTransferFunctionMinus->SetArgT(12,Abb.kmax);

	clKernel* clWienerW = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerW->BuildKernel();

	clWienerW->SetArgT(0,clW);
	clWienerW->SetArgT(1,clw);
	clWienerW->SetArgT(2,sizeX);
	clWienerW->SetArgT(3,sizeY);
	clWienerW->SetArgT(4,init);

	clKernel* clWienerWMinus = new clKernel(wienerwsource,context,cldev,"clWienerW",clq);
	clWienerWMinus->BuildKernel();

	clWienerWMinus->SetArgT(0,clWminus);
	clWienerWMinus->SetArgT(1,clwminus);
	clWienerWMinus->SetArgT(2,sizeX);
	clWienerWMinus->SetArgT(3,sizeY);
	clWienerWMinus->SetArgT(4,init);

	clKernel* clWienerV = new clKernel(wienervsource,context,cldev,"clWienerV",clq);
	clWienerV->BuildKernel();

	clWienerV->SetArgT(0,clV);
	clWienerV->SetArgT(1,clw);
	clWienerV->SetArgT(2,clwminus);
	clWienerV->SetArgT(3,sizeX);
	clWienerV->SetArgT(4,sizeY);
	clWienerV->SetArgT(5,init);

	clKernel* clWienerT = new clKernel(wienertsource,context,cldev,"clWienerT",clq);
	clWienerT->BuildKernel();

	clWienerT->SetArgT(0,clT);
	clWienerT->SetArgT(1,clw);
	clWienerT->SetArgT(2,clFFTImage1);
	clWienerT->SetArgT(3,sizeX);
	clWienerT->SetArgT(4,sizeY);
	clWienerT->SetArgT(5,init);

	clKernel* clWienerU = new clKernel(wienerusource,context,cldev,"clWienerU",clq);
	clWienerU->BuildKernel();

	clWienerU->SetArgT(0,clU);
	clWienerU->SetArgT(1,clwminus);
	clWienerU->SetArgT(2,clFFTImage1);
	clWienerU->SetArgT(3,sizeX);
	clWienerU->SetArgT(4,sizeY);
	clWienerU->SetArgT(5,init);
	
	clKernel* clMakeRestored = new clKernel(makerestoredsource,context,cldev,"clMakeRestored",clq);
	clMakeRestored->BuildKernel();

	clMakeRestored->SetArgT(0,clW);
	clMakeRestored->SetArgT(1,clWminus);
	clMakeRestored->SetArgT(2,clV);
	clMakeRestored->SetArgT(3,clT);
	clMakeRestored->SetArgT(4,clU);
	clMakeRestored->SetArgT(5,clFFTImage1);
	clMakeRestored->SetArgT(6,sizeX);
	clMakeRestored->SetArgT(7,sizeY);

	clKernel* clRotScale = new clKernel(RotScaleSource,context,cldev,"clRotScale",clq);
	clRotScale->BuildKernel();

	if(options.rotscale)
	{
		clRotScale->SetArgT(0,fullImage);
		clRotScale->SetArgT(1,rotScaleImage);
		clRotScale->SetArgT(2,xDim);
		clRotScale->SetArgT(3,yDim);
		clRotScale->SetArgT(4,options.magcal);
		clRotScale->SetArgT(5,options.rotcal);
	}

	size_t* fullWorkSize = new size_t[3];
	
	// There are 'elements' work items
	fullWorkSize[0] = xDim;
	fullWorkSize[1] = yDim;
	fullWorkSize[2] = 1;

	std::vector<std::complex<float>> copyImage(xDim*yDim) ;
	std::vector<float> rotscaledseries(xDim*yDim*numberOfImages);

	// Make a new rotation aligned stack based on the expected defocus step.
	// To prevent edges going missing should start at end that should be smallest so all images are only scaled upwards....


	
	// Rotate and Scale this image by expected defocus..
	if(options.rotscale)
	{
		int startingimage = 0;

		if(options.magcal > 1)
		{
			startingimage = 0;
		}
		else
		{
			startingimage = numberOfImages -1;
		}


		for(int im = 0; im < numberOfImages ; im++ )
		{
			float expecteddifference = (im-startingimage)*options.focalstep;
			clRotScale->SetArgT(6,expecteddifference);

			// Copy image into fullimage
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					copyImage[i+j*xDim] = seriesdata[im*xDim*yDim + i + (j)*xDim];
				}

			clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);
			clRotScale->Enqueue(fullWorkSize);
			clEnqueueReadBuffer(clq->cmdQueue,rotScaleImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					rotscaledseries[im*xDim*yDim +i + j*xDim] = copyImage[i + (j)*xDim].real();
				}
		}		
	} 
	else
	{
		// Fill new series with old series data.
		for(int im = 0; im < numberOfImages ; im++ )
			for(int j = 0; j < yDim;j++)
				for(int i = 0; i < xDim;i++)
				{
					rotscaledseries[im*xDim*yDim +i + j*xDim] = seriesdata[im*xDim*yDim +i + j*xDim];
				}
	}



	std::vector< std::complex< float > > dataOne( sizeX*sizeY );
	std::vector< std::complex< float > > dataTwo( sizeX*sizeY );

	int referenceimage = options.reference;
	int gonext = 1;
	int nextdir = -1;
	
	// First three have to be done via standard methods.

	// Initialise image numbers
	int imageone = referenceimage;
	int imagetwo = referenceimage;
	int currentimage;

	imagelist.push_back(referenceimage);

	// loop over all of the images
	for(int n = 0; n < numberOfImages - 1; n++)
	{
		Debug("Start Image "+boost::lexical_cast<std::string>(n));

		 // For RECONPCPCF part

		// First 2 registrations cannot be done against reconstruction
		// have to be done against reference image...
		if(n < 2)
		{
			imageone = referenceimage;
			imagetwo += (gonext*nextdir);

			while(! (imagetwo < numberOfImages && imagetwo >=0)) // loop over adjacent images.. i.e    10,11,9,12,8,13,7,14,6etc..
			{
				gonext++;
				nextdir *=-1;
				imagetwo += (gonext*nextdir);
			}

			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = (imagetwo-imageone) * options.focalstep;

			for(int j = 0; j < sizeY;j++)
				for(int i = 0; i < sizeX;i++)
				{
					dataOne[i+j*sizeX] = rotscaledseries[imageone*xDim*yDim + i + iLeft + (j+iTop)*xDim];
					dataTwo[i+j*sizeX] = rotscaledseries[imagetwo*xDim*yDim + i + iLeft + (j+iTop)*xDim];
				}

	
			clEnqueueWriteBuffer( clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
						0, NULL, NULL );

			FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			PCPCFLibrary.PhaseCompensatedPCF(numberoftrials,expecteddifference,options,clPCPCFResult,clImage1,sizeX,sizeY,
						imagetwo,dataOne,xShiftVals,yShiftVals,subXShifts,subYShifts,defocusshifts,0,0);

			imagelist.push_back(imagetwo);
		}

		if(n==2)
		{
			// To make sure it is initialised from end of first section
			currentimage = imagetwo;
		}
	
		if(n>=2)
		{
			// Now we can perform reconstruction and register to the reconstruction instead of standard pcpcf
			currentimage += (gonext*nextdir);

			while(! (currentimage < numberOfImages && currentimage >=0)) // not an actual image
			{
				gonext++;
				nextdir *=-1;
				currentimage += (gonext*nextdir);
			}

			// Need to reconstruct the previous images based on there focus and shifts relative to reference image...

			
			// Need to extract correct area from each image to perform reconstruction
			// 2 Possibilites - If an area of the correct size is not available then give an error...
			// or pad the area with zeroes for ALL images - need to first determine the extents of available area before we extract from any of the images.

			// Find maximum negative and positive shifts
			// All shifts have to be stored relative to the reference image.
			float maxnegx = 0;
			float maxnegy = 0;
			float maxposx = 0;
			float maxposy = 0;

			// Could just loop over this list of currently registered images instead...
			for(int i = 1; i <= imagelist.size(); i++)
			{
				if(subXShifts[imagelist[i-1]] < maxnegx)
					maxnegx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] < maxnegy)
					maxnegy = subYShifts[imagelist[i-1]];
				if(subXShifts[imagelist[i-1]] > maxposx)
					maxposx = subXShifts[imagelist[i-1]];
				if(subYShifts[imagelist[i-1]] > maxposy)
					maxposy = subYShifts[imagelist[i-1]];
			}

			// Determine amount to pad on either direction
			int padTop = 0;
			int padLeft = 0;
			int padRight = 0;
			int padBottom = 0;

			if (abs(maxnegx)-iLeft > 0)
				padLeft = ceil(abs(maxnegx)-iLeft);

			if (maxposx-iRight > 0)
				padRight = ceil(maxposx-iRight);

			if (abs(maxnegy)-iTop > 0)
				padTop = ceil(abs(maxnegy)-iTop);

			if (maxposy-iBottom > 0)
				padBottom = ceil(maxposy-iBottom);	

			for(int i = 0 ; i < imagelist.size() ; i++)
			{
				if(i == 0) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,init);
					clWienerWMinus->SetArgT(4,init);
					clWienerV->SetArgT(5,init);
					clWienerT->SetArgT(5,init);
					clWienerU->SetArgT(5,init);
				}
				if(i == 1) // Gets set after first image (i=0) is completed.
				{
					clWienerW->SetArgT(4,noinit);
					clWienerWMinus->SetArgT(4,noinit);
					clWienerV->SetArgT(5,noinit);
					clWienerT->SetArgT(5,noinit);
					clWienerU->SetArgT(5,noinit);
				}
				
				// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
				clPadCrop->SetArgT(4,padLeft);
				clPadCrop->SetArgT(5,padRight);
				clPadCrop->SetArgT(6,padTop);
				clPadCrop->SetArgT(7,padBottom);
				clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
				clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
				// Copy correct image into full image..
				std::vector<cl_float2> image(xDim*yDim);

				// Now add image to reconstruction - Make WTF - WTFminus			
				float defocus = defocusshifts[imagelist[i]];

				for(int j = 0 ; j < xDim*yDim ; j++)
				{
					image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
					image[j].s[1] = 0;
				}

				clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
				clPadCrop->Enqueue(globalWorkSize);

				image.clear();

				// Get FFT of this image
				FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

				clWaveTransferFunction->SetArgT(10,defocus);
				clWaveTransferFunctionMinus->SetArgT(10,defocus);

				clWaveTransferFunction->Enqueue(globalWorkSize);
				clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

				// add to W,W-,V,T,U 
				clWienerW->Enqueue(globalWorkSize);
				clWienerWMinus->Enqueue(globalWorkSize);
				clWienerV->Enqueue(globalWorkSize);
				clWienerT->Enqueue(globalWorkSize);
				clWienerU->Enqueue(globalWorkSize);
			}

			// make restored
			clMakeRestored->Enqueue(globalWorkSize);
			 
			Debug("Made Restored");

			// inverse fft for display
			FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);
			
			//Show Each reconstruction
			//Utility::PrintCLMemToImage(clRestored,"Restored",sizeX,sizeY,clFloat2,clq);
			
			// PCF with the reconstruction
			// Adjust reconstruction or do PCPCF?
			// TODO: at the moment reconstruction is in reference plane and images registered by PCPCF with this plane....
			// Also do i need to get modulus and FFT or not? - seems to work OK
			
			// Work out what focus difference is expected between these 2 images.
			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone

			// Here take an image thats already shifted approximately to the right place... then add this back onto calculated shifts...
			// Need to check image above or below to see if its registered yet

			bool prevImageRegistered = false;
			bool nextImageRegistered = false;

			for(int i = 1; i <= imagelist.size(); i++)
			{
				Debug("imglist: "+boost::lexical_cast<std::string>(imagelist[i-1]));
								
				if(imagelist[i-1]==currentimage-1)
				{
					Debug("PREV");
					prevImageRegistered = true;
					break;
				}
				if(imagelist[i-1]==currentimage+1)
				{
					Debug("NEXT");
					nextImageRegistered = true;
					break;
				}
			}

			int preshiftx = 0;
			int preshifty = 0;
			
			if(prevImageRegistered)
			{
				preshiftx = xShiftVals[currentimage-1];
				preshifty = yShiftVals[currentimage-1];
			}

			if(nextImageRegistered)
			{
				preshiftx = xShiftVals[currentimage+1];
				preshifty = yShiftVals[currentimage+1];
				Debug("xshift: "+boost::lexical_cast<std::string>(preshiftx));
				Debug("yshift: "+boost::lexical_cast<std::string>(preshifty));
			}
			

			Debug("Calculated Preshifts");

			// Cant have it checking an image area that doesnt exist anymore
			if( iLeft + preshiftx < 0)
				preshiftx -=( preshiftx + iLeft );

			if( iLeft + preshiftx + sizeX >= xDim)
				preshiftx -=( iLeft + preshiftx + sizeX - xDim);

			if( iTop + preshifty < 0)
				preshifty -=( preshifty + iTop );

			if( iTop + preshifty + sizeY >= yDim)
				preshiftx -=( iTop + preshifty + sizeY - yDim);


			for(int j = 0; j < sizeY; j++)
				for(int i = 0; i < sizeX; i++)
				{
					dataTwo[i+j*sizeX] = rotscaledseries[currentimage*xDim*yDim + i + iLeft + preshiftx + (j+iTop + preshifty)*xDim];
				}

			Debug("Made Preshifted image");

			clEnqueueWriteBuffer( clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataTwo[ 0 ], 0, NULL, NULL );

			FFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

			float expecteddifference = (currentimage-referenceimage) * options.focalstep;

			// Get number of trial steps
			// Assume one if not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			PCPCFLibrary.PhaseCompensatedPCF(numberoftrials,expecteddifference,options,clPCPCFResult,clImage1,sizeX,sizeY,
						currentimage,dataOne,xShiftVals,yShiftVals,subXShifts,subYShifts,defocusshifts, preshiftx, preshifty);

			imagelist.push_back(currentimage);

		}
		
		// to get to next image..
		gonext++;
		nextdir *=-1;
		// image = previous +(gonext*nextdir);
		// checkbounds then if not available do again.
		
	}

	Debug("Starting Final Reconstruction");


	// Now rebuild reconstruction including the final registered image 
	// All shifts have to be stored relative to the reference image.
	float maxnegx = 0;
	float maxnegy = 0;
	float maxposx = 0;
	float maxposy = 0;

	// Could just loop over this list of currently registered images instead...
	for(int i = 1; i <= numberOfImages; i++)
	{
		if(subXShifts[i-1] < maxnegx)
			maxnegx = subXShifts[i-1];
		if(subYShifts[i-1] < maxnegy)
			maxnegy = subYShifts[i-1];
		if(subXShifts[i-1] > maxposx)
			maxposx = subXShifts[i-1];
		if(subYShifts[i-1] > maxposy)
			maxposy = subYShifts[i-1];
	}

	// Determine amount to pad on either direction
	int padTop = 0;
	int padLeft = 0;
	int padRight = 0;
	int padBottom = 0;

	if (abs(maxnegx)-iLeft > 0)
		padLeft = ceil(abs(maxnegx)-iLeft);

	if (maxposx-iRight > 0)
		padRight = ceil(maxposx-iRight);

	if (abs(maxnegy)-iTop > 0)
		padTop = ceil(abs(maxnegy)-iTop);

	if (maxposy-iBottom > 0)
		padBottom = ceil(maxposy-iBottom);


	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]];

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();

		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);

		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);
			 
	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Restored EW",sizeX,sizeY,clFloat2,clq);

	// Now all should be registered - produce the final restored image and any other relevant data like Shift Graphs, and Drift Corrected Stack
	DigitalMicrograph::Image driftcorrected = DigitalMicrograph::RealImage("Drift Corrected Stack",4,sizeX,sizeY,numberOfImages);
	Gatan::PlugIn::ImageDataLocker driftlocker = Gatan::PlugIn::ImageDataLocker(driftcorrected);
	float* driftstack = (float*) driftlocker.get();

	std::vector<cl_float2> cropimage(sizeX*sizeY);

	for(int i = 0 ; i < numberOfImages ; i++)
	{
		int padTop = 0;
		int padLeft = 0;
		int padRight = 0;
		int padBottom = 0;

		if (abs(maxnegx)-iLeft > 0)
			padLeft = ceil(abs(maxnegx)-iLeft);

		if (maxposx-iRight > 0)
			padRight = ceil(maxposx-iRight);

		if (abs(maxnegy)-iTop > 0)
			padTop = ceil(abs(maxnegy)-iTop);

		if (maxposy-iBottom > 0)
			padBottom = ceil(maxposy-iBottom);

		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[i]);
		clPadCrop->SetArgT(9,subYShifts[i]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(xDim*yDim);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[i];
			
		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[i*xDim*yDim + j];
			image[j].s[1] = 0;
		}

		// Copy full image to GPU and run crop kernel
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		// Copy image data back to host
		clEnqueueReadBuffer(clq->cmdQueue,clImage1,CL_TRUE,0,sizeX*sizeY*sizeof(cl_float2),&cropimage[0],0,NULL,NULL);
		
		// Copy to drift corrected stack
		for(int k = 0 ; k < sizeX*sizeY ; k++)
		{
			driftstack[i*sizeX*sizeY + k] = cropimage[k].s[0];
		}

		image.clear();
	}

	cropimage.clear();
	driftlocker.~ImageDataLocker();

	// Display drift corrected
	driftcorrected.GetOrCreateImageDocument().Show();

	dataOne.clear();
	dataTwo.clear();
	
	// Don't need these anymore
	PCPCF->~clKernel();
	clReleaseMemObject(clFFTImage2);
	clReleaseMemObject(clPCPCFResult);
	clReleaseMemObject(clImage2);

	// Now start working on finding the actual astigmatism and defocus from the reconstruction....

	// Make additional kernels
	clKernel* clGetQ				= new clKernel(getQsource, context, cldev, "clCalculateQ", clq);
	clKernel* clSumReduction		= new clKernel(sumReductionsource,context, cldev, "clSumReduction", clq);
	clKernel* clMinusWavefunction	= new clKernel(minuswavefunctionsource, context, cldev, "clMinusWavefunction", clq);
	clKernel* clCalculatePCI		= new clKernel(getPCIsource,context, cldev, "clCalculatePCI", clq);
	clGetQ->BuildKernel();
	clSumReduction->BuildKernel();
	clMinusWavefunction->BuildKernel();
	clCalculatePCI->BuildKernel();

	// New memory required
	cl_mem clRestoredMinus	= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clQ				= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCI			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIC			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIM			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);
	cl_mem clPCIS			= clCreateBuffer(context,CL_MEM_READ_WRITE,sizeX*sizeY*sizeof(std::complex<float>),0,&status);

	// Set Kernel Arguments
	clGetQ->SetArgT(0,clW);
	clGetQ->SetArgT(1,clWminus);
	clGetQ->SetArgT(2,clV);
	clGetQ->SetArgT(3,clQ);
	clGetQ->SetArgT(4,sizeX);
	clGetQ->SetArgT(5,sizeY);

	clMinusWavefunction->SetArgT(0,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clMinusWavefunction->SetArgT(1,clRestoredMinus);
	clMinusWavefunction->SetArgT(2,sizeX);
	clMinusWavefunction->SetArgT(3,sizeY);

	clGetQ->Enqueue(globalWorkSize);
	clMinusWavefunction->Enqueue(globalWorkSize);

	// Get sum of Q by my awesome reduction kernel ¬_¬ 
	int totalSize = sizeX*sizeY;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	cl_mem clSumOutput = clCreateBuffer(context,CL_MEM_READ_WRITE,nGroups*sizeof(std::complex<float>),0,&status);
	std::vector< std::complex< float > > sumQs( nGroups );

	clSumReduction->SetArgT(0,clQ);
	clSumReduction->SetArgT(1,clSumOutput);
	clSumReduction->SetArgT(2,totalSize);
	clSumReduction->SetArgLocalMemory(3,256,clFloat2);

	clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 0, NULL, NULL );


	// Find out which numbers to read back
	float sumQ2 = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sumQ2 += sumQs[i].real();
	}

	// Trial lots of C1 Values
	int trials2 = 50;

	float* Fpci = new float[trials2];
	float* Fpcic = new float[trials2];
	float* Fpcis = new float[trials2];
	float* Fpcim = new float[trials2];
	float* Fpci2 = new float[trials2];

	float A1r = 0.0f;
	float A1i = 0.0f;
	float zero2 = 0.0f;
	
	float StartDf = -50;
	float StepDf = 2;

	//DM::Image pcidisplay = DM::RealImage("PCI",4,sizeX,sizeY);
	//pcidisplay.GetOrCreateImageDocument().Show();
	//pcidisplay.GetImageDisplay(0).SetOutlierTrimLimits(0,0);
	//pcidisplay.GetImageDisplay(0).SetSurveyTechnique(1);

	clCalculatePCI->SetArgT(0,clQ);
	clCalculatePCI->SetArgT(1,clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	clCalculatePCI->SetArgT(2,clRestoredMinus);
	clCalculatePCI->SetArgT(3,clxFrequencies);
	clCalculatePCI->SetArgT(4,clyFrequencies);
	clCalculatePCI->SetArgT(5,clPCI);
	clCalculatePCI->SetArgT(6,clPCIC);
	clCalculatePCI->SetArgT(7,clPCIM);
	clCalculatePCI->SetArgT(8,clPCIS);
	clCalculatePCI->SetArgT(9,sizeX);
	clCalculatePCI->SetArgT(10,sizeY);
	clCalculatePCI->SetArgT(12,zero2);
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);
	clCalculatePCI->SetArgT(15,wavelength);
	clCalculatePCI->SetArgT(16,Abb.kmax);
	

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
 
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);

		// DO PCI Function and get sums for pcic,s,m etc here.
		clCalculatePCI->Enqueue(globalWorkSize);
		
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpci[trial] = fpcitrial/sumQ2;
		Fpcic[trial] = fpcictrial/sumQ2;
		Fpcis[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcim[trial] = fpcimtrial/sumQ2;
	}

	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = StartDf + trial*StepDf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < trials2 ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial <trials2 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess)
		);


	DigitalMicrograph::Image pcigraph = DigitalMicrograph::RealImage("PCI",4,trials2);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial <trials2 ; trial++)
	{
		pciData[trial] = Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	DigitalMicrograph::UpdateImage(pcigraph);


	delete[] Fpci;
	delete[] Fpci2;
	delete[] Fpcic;
	delete[] Fpcim;
	delete[] Fpcis;

	// First PCI Stage Complete

	// Redo more accurate PCI in correct area 
	// Trial lots of C1 Values

	int trials3 = 50;
	
	float* Fpciref = new float[trials3];
	float* Fpcicref = new float[trials3];
	float* Fpcisref = new float[trials3];
	float* Fpcimref = new float[trials3];
	float* Fpci2ref = new float[trials3];

	A1r = A1rGuess;
	A1i = A1iGuess;
	
	StartDf = DefocusGuess-10;
	StepDf = 21.0f/trials2;

	// Set Kernel Arguments - DONE ALREADY NOT CHANGED
	clCalculatePCI->SetArgT(13,A1r);
	clCalculatePCI->SetArgT(14,A1i);

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		float C1 = StartDf + trial*StepDf;

		clCalculatePCI->SetArgT(11,C1);

		// DO PCI Function and get sums for pcic,s,m etc here.
		clCalculatePCI->Enqueue(globalWorkSize);
		clSumReduction->SetArgT(0,clPCI);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );
		
		// Find out which numbers to read back
		float fpcitrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcitrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIC);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcictrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcictrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIS);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcistrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcistrial += sumQs[i].real();
		}

		clSumReduction->SetArgT(0,clPCIM);
		clSumReduction->Enqueue3D(globalSizeSum,localSizeSum);
		
		// Now copy back 
		clEnqueueReadBuffer( clq->cmdQueue, clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sumQs[0], 
					0, NULL, NULL );

		
		// Find out which numbers to read back
		float fpcimtrial = 0;

		for(int i = 0 ; i < nGroups; i++)
		{
			fpcimtrial += sumQs[i].real();
		}

		Fpciref[trial] = fpcitrial/sumQ2;
		Fpcicref[trial] = fpcictrial/sumQ2;
		Fpcisref[trial] = fpcistrial/sumQ2;
		// Not sure if this is ok, should do sumQ inside loop.
		Fpcimref[trial] = fpcimtrial/sumQ2;
	}

	float maxheighta = 0;
	float C1aa = 0;
	int besttriala = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		if(Fpcimref[trial] > maxheighta)
		{
			maxheighta = Fpcimref[trial];
			C1aa = StartDf + trial*StepDf;
			besttriala = trial;
		}
	}

	float astigaxisa = atan2(Fpcisref[besttriala],Fpcicref[besttriala])/2;

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		Fpci2ref[trial] = Fpciref[trial] - Fpcicref[trial] * cos(2*astigaxisa) - Fpcisref[trial]*sin(2*astigaxisa);
	}

	float maxheight2a = 0;
	float C1ba = 0;

	for(int trial = 0 ; trial < trials3 ; trial++)
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

	Utility::SetProgressWindow(
		"defocus: "+boost::lexical_cast<std::string>(DefocusGuess2),
		"A1 real: "+boost::lexical_cast<std::string>(A1rGuess2),
		"A1 imag: "+boost::lexical_cast<std::string>(A1iGuess2)
		);

	// Print out some results and graphs

	DigitalMicrograph::Image pcigraph2 = DigitalMicrograph::RealImage("PCI - Refined",4,trials3);
	Gatan::PlugIn::ImageDataLocker pciLocker2(pcigraph2);
	float* pciData2 = (float*) pciLocker2.get();

	for(int trial = 0 ; trial < trials3 ; trial++)
	{
		pciData2[trial]=Fpciref[trial];
	}

	pciLocker2.~ImageDataLocker();
	pcigraph2.GetOrCreateImageDocument().Show();
	pcigraph2.SetDimensionCalibration(0,StartDf,StepDf,"nm",0);
	

	delete[] Fpciref;
	delete[] Fpcicref;
	delete[] Fpcisref;
	delete[] Fpcimref;
	delete[] Fpci2ref;

	// Second PCI Stage Complete

	// Then do restoration again at correct focal/astig value

	for(int i = 0 ; i < imagelist.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,init);
			clWienerWMinus->SetArgT(4,init);
			clWienerV->SetArgT(5,init);
			clWienerT->SetArgT(5,init);
			clWienerU->SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			clWienerW->SetArgT(4,noinit);
			clWienerWMinus->SetArgT(4,noinit);
			clWienerV->SetArgT(5,noinit);
			clWienerT->SetArgT(5,noinit);
			clWienerU->SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		clPadCrop->SetArgT(4,padLeft);
		clPadCrop->SetArgT(5,padRight);
		clPadCrop->SetArgT(6,padTop);
		clPadCrop->SetArgT(7,padBottom);
		clPadCrop->SetArgT(8,subXShifts[imagelist[i]]);
		clPadCrop->SetArgT(9,subYShifts[imagelist[i]]);
				
		// Copy correct image into full image..

		std::vector<cl_float2> image(xDim*yDim);

		// Want to rotate and scale correct this image -- do rotations from reference defocus not absolute defocus or it will (probably) affect shifts.

		float expecteddifference = (imagelist[i]-referenceimage) * options.focalstep;

		for(int j = 0 ; j < xDim*yDim ; j++)
		{
			image[j].s[0] = rotscaledseries[imagelist[i]*xDim*yDim + j];
			image[j].s[1] = 0;
		}
	
		clEnqueueWriteBuffer(clq->cmdQueue,fullImage,CL_TRUE,0,xDim*yDim*sizeof(cl_float2),&image[0],0,NULL,NULL);
		clPadCrop->Enqueue(globalWorkSize);

		image.clear();


		// Get FFT of this image
		FFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[imagelist[i]] + DefocusGuess2;

		Utility::SetResultWindow("Adding Image "+boost::lexical_cast<std::string>(imagelist[i])+" at defocus "+boost::lexical_cast<std::string>(defocus)+"\n");

		clWaveTransferFunction->SetArgT(10,defocus);
		clWaveTransferFunctionMinus->SetArgT(10,defocus);
		clWaveTransferFunction->SetArgT(8,A1rGuess2);
		clWaveTransferFunctionMinus->SetArgT(8,A1rGuess2);
		clWaveTransferFunction->SetArgT(9,A1iGuess2);
		clWaveTransferFunctionMinus->SetArgT(9,A1iGuess2);
		clWaveTransferFunction->Enqueue(globalWorkSize);
		clWaveTransferFunctionMinus->Enqueue(globalWorkSize);

		// add to W,W-,V,T,U 
		clWienerW->Enqueue(globalWorkSize);
		clWienerWMinus->Enqueue(globalWorkSize);
		clWienerV->Enqueue(globalWorkSize);
		clWienerT->Enqueue(globalWorkSize);
		clWienerU->Enqueue(globalWorkSize);
	}

	// make restored
	clMakeRestored->Enqueue(globalWorkSize);

	// inverse fft for display
	FFT->Enqueue(clFFTImage1,clRestored,CLFFT_BACKWARD);

	Utility::PrintCLMemToImage(clRestored,"Final Restored EW",sizeX,sizeY,clFloat2,clq);

	// New Image added for doing rotation and scale correcting
	copyImage.clear();
	rotscaledseries.clear();

	clReleaseMemObject(clImage1);
	clReleaseMemObject(clFFTImage1);
	clReleaseMemObject(fullImage);
	clReleaseMemObject(clw);
	clReleaseMemObject(clwminus);
	clReleaseMemObject(clW);
	clReleaseMemObject(clWminus);
	clReleaseMemObject(clU);
	clReleaseMemObject(clV);
	clReleaseMemObject(clT);
	clReleaseMemObject(clRestored);
	clReleaseMemObject(clRestoredMinus);
	clReleaseMemObject(clQ);
	clReleaseMemObject(clPCI);
	clReleaseMemObject(clPCIC);
	clReleaseMemObject(clPCIM);
	clReleaseMemObject(clPCIS);
	clReleaseMemObject(clSumOutput);

	clPadCrop->~clKernel();
	clWienerW->~clKernel();
	clWienerWMinus->~clKernel();
	clWienerV->~clKernel();
	clWienerT->~clKernel();
	clWienerU->~clKernel();
	clWaveTransferFunction->~clKernel();
	clWaveTransferFunctionMinus->~clKernel();
	clMakeRestored->~clKernel();

	clGetQ->~clKernel();
	clSumReduction->~clKernel();
	clMinusWavefunction->~clKernel();
	clCalculatePCI->~clKernel();
};
*/