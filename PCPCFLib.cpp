#include "StdAfx.h"
#include "PCPCFLib.h"


void WrappedImageOffsets(int pos, int sizeX, int sizeY, int &t, int &b, int &l , int &r)
{
	// Translate linear array index into row and column.
	int indexr = floor(float((pos)/(sizeX)))+1;
	int indexc = pos + 1 - sizeX*floor(float((pos)/(sizeX)));
	
	/*
	t = indexr - 1 + sizeX*(indexc-2);
	b = indexr - 1 + sizeX*(indexc);
	l = indexr - 2 + sizeX*(indexc-1);
	r = indexr + sizeX*(indexc-1);
	*/

	t = indexc - 1 + sizeX*(indexr-2);
	b = indexc - 1 + sizeX*(indexr);
	l = indexc - 2 + sizeX*(indexr-1);
	r = indexc + sizeX*(indexr-1);

	// in top row so top wraps to bottom
	if(indexr == 1)
	{
		t = pos + (sizeX*sizeY) - sizeX;
	}

	// in bottom row so bottom wraps to top
	if(indexr == sizeY)
	{
		b= pos - (sizeX*sizeY) + sizeX;
	}

	// in far left so left wraps to right
	if(indexc == 1)
	{
		l += sizeX;
	}

	// in far right so right wraps to left
	if(indexc == sizeX)
	{
		r -= sizeX;
	}
}

PCPCFLib::PCPCFLib(void)
{
}


PCPCFLib::~PCPCFLib(void)
{
}

/* Doesn't account for possibility that max is at edge of image
void PCPCFLib::FindVertexParabola(float xl, float centre, float xr, float yt, float yb, float &xoffset, float &yoffset, float &peakheight)
{
	float x1 = -1;
	float x2 = 0;
	float x3 = 1;
	float y1 = xl;
	float y2 = centre;
	float y3 = xr;

	double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	double B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	xoffset = -B / (2*A);
	float peak1 = C - B*B / (4*A);

	y1 = yt;
	y3 = yb;

	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	yoffset = -B / (2*A);
	float peak2 = C - B*B / (4*A);

	peakheight = (peak1 + peak2) / 2;
}
*/


// BUG: Giving negative peak heights making comparison go wrong.... (possibly due to abs, or i used abs to try fix this?)
void PCPCFLib::FindVertexParabola(int maxposition1, int sizeX, int sizeY, std::vector<std::complex<float>> &data, float &xoffset, float &yoffset, float &peakheight)
{
	int t;
	int l;
	int b;
	int r;

	WrappedImageOffsets(maxposition1,sizeX,sizeY,t,b,l,r);

	// Think this could go badly if it was right at the edges of images, will need to add wrapping to positions
	float x1 = -1;
	float x2 = 0;
	float x3 = 1;
	float y1 = data[l].real();
	//float y1 = hypot(data[l].real(),data[l].imag());
	float y2 = data[maxposition1].real();
	//float y2 = hypot(data[maxposition1].real(),data[maxposition1].imag());
	float y3 = data[r].real();
	//float y3 = hypot(data[r].real(),data[r].imag());


	float denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	float A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	float B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	float C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	xoffset = max(min(-B / (2*A),1.0f),-1.0f); // Bad fits cause subshifts that are way bigger than +-1.0
	float peak1 = C - B*B / (4*A);

	y1 = data[t].real();
	//y1 = hypot(data[t].real(),data[t].imag());
	y3 = data[b].real();
	//y3 = hypot(data[b].real(),data[b].imag());


	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	yoffset = max(min(-B / (2*A),1.0f),-1.0f);
	float peak2 = C - B*B / (4*A);

	peakheight = (peak1 + peak2) / (2.0f);

	//Utility::SetResultWindow(boost::lexical_cast<std::string>(xoffset)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(yoffset)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(peakheight)+"\n");

}

void PCPCFLib::FindVertexParabolaMI(int maxposition1, int sizeX, int sizeY, float* data, float &xoffset, float &yoffset, float &peakheight)
{
	int t;
	int l;
	int b;
	int r;

	WrappedImageOffsets(maxposition1,sizeX,sizeY,t,b,l,r);

	// Think this could go badly if it was right at the edges of images, will need to add wrapping to positions
	float x1 = -1;
	float x2 = 0;
	float x3 = 1;
	float y1 = data[l];
	float y2 = data[maxposition1];
	float y3 = data[r];

	float denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	float A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	float B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	float C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	xoffset = max(min(-B / (2*A),1.0f),-1.0f); 
	float peak1 = C - B*B / (4*A);

	y1 = data[t];
	y3 = data[b];

	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	yoffset = max(min(-B / (2*A),1.0f),-1.0f); 
	float peak2 = C - B*B / (4*A);

	peakheight = (peak1 + peak2) / (2.0f);

	//Utility::SetResultWindow(boost::lexical_cast<std::string>(xoffset)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(yoffset)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(peakheight)+"\n");

}

void PCPCFLib::GetPeak(float* peak,int maxindexc,int maxindexr,int sizeX,int sizeY, std::vector<std::complex<float>> &data, float &xpos, float &ypos)
{
	//NOTE SHOULD DO ABS NOT REAL
	// NOTE NOT SURE THIS IS TRUE AT ALL...
	float sum = 0;
	float xweight = 0;
	float yweight = 0;

	for(int i = -2 ; i<=2 ; i++)
		for(int j = -2 ; j<=2 ; j++)
		{
			if(maxindexc+i < 0)
				maxindexc+=sizeX;
			if(maxindexc+i > sizeX-1)
				maxindexc-=sizeX;
			if(maxindexr+j < 0)
				maxindexr+=sizeY;
			if(maxindexr+j > sizeY - 1)
				maxindexr-=sizeY;

			maxindexc+=i;
			maxindexr+=j;

			sum+= sqrt(data[maxindexc + maxindexr+sizeX].real() *  data[maxindexc + maxindexr+sizeX].real() +  data[maxindexc + maxindexr+sizeX].imag() *  data[maxindexc + maxindexr+sizeX].imag());
			xweight += i * sqrt(data[maxindexc + maxindexr+sizeX].real() *  data[maxindexc + maxindexr+sizeX].real() +  data[maxindexc + maxindexr+sizeX].imag() *  data[maxindexc + maxindexr+sizeX].imag());
			yweight += j * sqrt(data[maxindexc + maxindexr+sizeX].real() *  data[maxindexc + maxindexr+sizeX].real() +  data[maxindexc + maxindexr+sizeX].imag() *  data[maxindexc + maxindexr+sizeX].imag());
		}

	xpos = xweight / sum;
	ypos = yweight / sum;

}

void PCPCFLib::GetShifts(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight, std::vector<std::complex<float>> &data, int sizeX, int sizeY)
{
	// Sometimes the highest value within the range set by maxshift is still negative...
	maxheight = -FLT_MAX;

	int maxPosition1 = 0;

	for(int j = 0; j< sizeX*sizeY;j++)
	{	
		//float val =  sqrt(data[j].real()*data[j].real() + data[j].imag()*data[j].imag());
		float val = data[j].real();
	//	pcpcfdata[j] = val;

		if(val > maxheight)
		{
			maxheight = val;
			maxPosition1 = j;
		}
	}

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
				
	// Construct 5*5 region around max to do peakfit (sortof). Probably just CoM. Can do vertex thing also
	//float* peak = new float[25];
	//GetPeak(peak,maxindexc,maxindexr,sizeX,sizeY, data, subXShift, subYShift);
	//delete[] peak;

	//Utility::SetResultWindow(boost::lexical_cast<std::string>(xShift)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(yShift)+"\n");

	// Parabola Vertex Mode
	PCPCFLib::FindVertexParabola(maxPosition1,sizeX,sizeY, data, subXShift, subYShift, maxheight);


}

void PCPCFLib::GetShifts(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight, 
						 std::vector<std::complex<float>> &data, int sizeX, int sizeY, float maxshift)
{
	// Sometimes the highest value within the range set by maxshift is still negative...
	maxheight = -FLT_MAX;
	int maxPosition1 = 0;

	// Modify to only look for a max with a certain range.
	int testpixels = 0;

	for(int j = 0; j< sizeY;j++)
		for(int i = 0; i< sizeX;i++)
		{	
			if(maxshift == 0)
			{
				//float val =  sqrt(data[i+j*sizeX].real()*data[i+j*sizeX].real() + data[i+j*sizeX].imag()*data[i+j*sizeX].imag());
				float val =  data[i+j*sizeX].real();

				if(val > maxheight)
				{
					maxheight = val;
					maxPosition1 = i+j*sizeX;
				}
			}
			else
			{		
				// This selects a square not a circular range :(
				//if(i<maxshift || (i+maxshift)>=sizeX)
					//if(j<maxshift || (j+maxshift)>=sizeY)
				if(i*i + j*j <= maxshift*maxshift || (i-sizeX)*(i-sizeX)+ j*j <= maxshift*maxshift || (i-sizeX)*(i-sizeX) + (j-sizeY)*(j-sizeY) <= maxshift*maxshift || i*i + (j-sizeY)*(j-sizeY) <= maxshift*maxshift)
					{
						//float val =  sqrt(data[i+j*sizeX].real()*data[i+j*sizeX].real() + data[i+j*sizeX].imag()*data[i+j*sizeX].imag());
						testpixels++;
						float val =  data[i+j*sizeX].real();
						//	pcpcfdata[j] = val;

						if(val > maxheight)
						{
							Debug("found \n");
							maxheight = val;
							maxPosition1 = i+j*sizeX;
						}
					}
			}
		}

	Debug(Lex(testpixels)+"\n");

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
				
	// Construct 5*5 region around max to do peakfit (sortof). Probably just CoM. Can do vertex thing also
	//float* peak = new float[25];
	//GetPeak(peak,maxindexc,maxindexr,sizeX,sizeY, data, subXShift, subYShift);
	//delete[] peak;

	//Utility::SetResultWindow(boost::lexical_cast<std::string>(xShift)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(yShift)+"\n");

	// Parabola Vertex Mode
	PCPCFLib::FindVertexParabola(maxPosition1,sizeX,sizeY, data, subXShift, subYShift, maxheight);

	Debug("subx "+Lex(subXShift)+"\n");
}

void PCPCFLib::GetShiftsCL(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight, 
						 std::vector<std::complex<float>> &data, cl_mem &cldata, int sizeX, int sizeY, float maxshift, clKernel &clMaxReduction, cl_mem &result, cl_mem &position)
{
	// Sometimes the highest value within the range set by maxshift is still negative...
	maxheight = -FLT_MAX;
	int maxPosition1 = 0;

	// Need to set memory and workgroup sizes sensibly for this kernel...

	int length = sizeX*sizeY;

	int ReductionGroups =  256;

	size_t* WorkSize = new size_t[3];
	WorkSize[0] = ReductionGroups*256;
	WorkSize[1] = 1;
	WorkSize[2] = 1;
	
	size_t* LocalSize = new size_t[3];
	LocalSize[0] = 256;
	LocalSize[1] = 1;
	LocalSize[2] = 1;


	clMaxReduction.SetArgT(0,cldata);
	clMaxReduction.SetArgLocalMemory(1,256,clFloat);
	clMaxReduction.SetArgLocalMemory(2,256,clUInt);
	clMaxReduction.SetArgT(3,length);
	clMaxReduction.SetArgT(4,result);
	clMaxReduction.SetArgT(5,position);

	clMaxReduction.Enqueue3D(WorkSize,LocalSize);

	// Now read back 256 maxs and there offsets and find best one.
	std::vector<float> hostres(ReductionGroups);
	std::vector<cl_uint> hostpos(ReductionGroups);

	clEnqueueReadBuffer( clState::clq->cmdQueue, result, CL_FALSE, 0, ReductionGroups*sizeof(float) , &hostres[0], 0, NULL, NULL );
	clEnqueueReadBuffer( clState::clq->cmdQueue, position, CL_TRUE, 0, ReductionGroups*sizeof(unsigned int) , &hostpos[0], 0, NULL, NULL );
	
	// Find out which numbers to read back
	float res = -1000000.0f;
	cl_uint pos = 0;

	for(int i = 0 ; i < ReductionGroups; i++)
	{
		if(hostres[i] > res)
		{
			res = hostres[i];
			pos = hostpos[i];

			//Debug(Lex(res));
			//Debug(Lex(pos));
		}
	}

	maxPosition1 = pos;

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
				
	// Construct 5*5 region around max to do peakfit (sortof). Probably just CoM. Can do vertex thing also
	//float* peak = new float[25];
	//GetPeak(peak,maxindexc,maxindexr,sizeX,sizeY, data, subXShift, subYShift);
	//delete[] peak;

	//Utility::SetResultWindow(boost::lexical_cast<std::string>(xShift)+"\n");
	//Utility::SetResultWindow(boost::lexical_cast<std::string>(yShift)+"\n");

	// Parabola Vertex Mode
	PCPCFLib::FindVertexParabola(maxPosition1,sizeX,sizeY, data, subXShift, subYShift, maxheight);

	//Debug("subx "+Lex(subXShift)+"\n");
}

void PCPCFLib::GetShiftsMI(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight, 
						 float* data, int sizeX, int sizeY, float maxshift)
{
	
	maxheight = -FLT_MAX;
	int maxPosition1 = 0;

	// Modify to only ook for a max with a certain range.
	for(int j = 0; j< sizeY;j++)
		for(int i = 0; i< sizeX;i++)
		{	
			if(i!=sizeX/2 && j!=sizeY/2)
			{
				float val = data[i+j*sizeX];

				if(val > maxheight)
				{
					maxheight = val;
					maxPosition1 = i+j*sizeX;
				}
			}
		}

	// Translate linear array index into row and column.
	int maxindexr = floor(float((maxPosition1)/(sizeX)))+1;
	int maxindexc = maxPosition1 + 1 - sizeX*floor(float((maxPosition1)/(sizeX)));


	// The zero position is at INDEX width/2 and height/2 i.e 16th element if width is 30...


	// Shift is positive or negative depending on image quadrant it appears in.
	yShift = maxindexr - (1+sizeX/2); // 0 if max is at 31st element of the 60..
	xShift = maxindexc - (1+sizeY/2);
				

	// Parabola Vertex Mode
	PCPCFLib::FindVertexParabolaMI(maxPosition1,sizeX,sizeY, data, subXShift, subYShift, maxheight);
}

void PCPCFLib::PrepareImages(float* Data, int width, int height, int numberofimages)
{
	// Loop over each image to find the average value and divide through by it.
	// Also determine location of hot pixels and remove them by comparison with median filter

	DigitalMicrograph::Image source = DigitalMicrograph::RealImage("source",4,width,height);
	
	for(int z = 0 ; z < numberofimages ; z++)
	{

		Gatan::PlugIn::ImageDataLocker SourceLocker(source);
		float* sourcedata = (float*) SourceLocker.get();

		for(int j = 0 ; j < height ; j++)
			for(int i = 0 ; i < width ; i++)
			{
				sourcedata[i+j*width] = Data[z*width*height + i+j*width];
			}

		DigitalMicrograph::Image med = DigitalMicrograph::MedianFilter(source,3,1);
		
		Gatan::PlugIn::ImageDataLocker MedLocker(med);
		float* meddata = (float*) MedLocker.get();

		float imagesum = 0.0f;

		for(int j = 0 ; j < height ; j++)
			for(int i = 0 ; i < width ; i++)
			{
				if(sourcedata[i+j*width]/meddata[i+j*width] > 1.5f || sourcedata[i+j*width]/meddata[i+j*width] < 0.5f)
				{
					sourcedata[i+j*width] = meddata[i+j*width];
				}

				imagesum+=sourcedata[i+j*width];
			}

		SourceLocker.~ImageDataLocker();
		MedLocker.~ImageDataLocker();

		for(int j = 0 ; j < height ; j++)
			for(int i = 0 ; i < width ; i++)
			{
				// Normalise to average of 1, then -1.
				Data[z*width*height + i+j*width] = (Data[z*width*height + i+j*width]/(imagesum/(width*height)))-1.0f;
			}

		// Clear med image...
		//DigitalMicrograph::DeleteImage(med);
	}

	// Clear source image....
	DigitalMicrograph::DeleteImage(source);
}
/*
void PCPCFLib::PhaseCompensatedPCF(int numberoftrials, float expectedDF, PCPCFOptions options, cl_mem &clPCPCFResult, cl_mem &clImage1, int sizeX, int sizeY, int imagenumber, std::vector<std::complex<float>> &dataOne, int* xShiftVals, int* yShiftVals, float* subXShifts, float* subYShifts, float* defocusshifts, int preshiftx, int preshifty)
{
	// Create arrays to hold results of each individual trial run
	int* trialxshifts = new int[numberoftrials];
	int* trialyshifts = new int[numberoftrials];
	float* trialsubxshifts = new float[numberoftrials];
	float* trialsubyshifts = new float[numberoftrials];
	float* trialdefocus = new float[numberoftrials];
	float* peakheights = new float[numberoftrials];

	for(int trial = 0; trial < numberoftrials; trial++)
	{

		Utility::SetResultWindow("Trial " + Lex(trial+1) + " of " + Lex(numberoftrials) + "\n");

		// go from expected - searchpercentage to expected + searchpercentage
		float trialdef = expectedDF - (options.searchpercentage/100.0f)*options.focalstep 
					+ ((float)trial/((float)numberoftrials-1.0f))*2.0f*(options.searchpercentage/100.0f)*options.focalstep; 

		// above formula wont work when number of trials is 1
		if(numberoftrials == 1)
		{
			trialdef = expectedDF;
		}

		PCPCF->SetArgT(7,trialdef);
		PCPCF->Enqueue(globalWorkSize);

		// Now Inverse FFT
		// Using memory reserved for storing images as it is not needed at this point.
		FFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);

		clEnqueueReadBuffer( clState::clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &dataOne[ 0 ], 
					0, NULL, NULL );

		Utility::SetResultWindow("Registering Image "+boost::lexical_cast<std::string>(imagenumber)+" at defocus "+boost::lexical_cast<std::string>(trialdef)+"\n");

		// Find shift from max height position
		int xShift;
		int yShift;
		float subXShift;
		float subYShift;
		float maxHeight1;

		// Translate linear array index into row and column.
		PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,sizeX,sizeY,options.maxdrift);

		trialxshifts[trial] = preshiftx + xShift;
		trialyshifts[trial] = preshifty+ yShift;
		trialsubxshifts[trial] =  preshiftx + xShift + subXShift;
		trialsubyshifts[trial] = preshifty + yShift + subYShift;
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
			xShiftVals[imagenumber] = trialxshifts[trial];
			yShiftVals[imagenumber] = trialyshifts[trial];
			subXShifts[imagenumber] = trialsubxshifts[trial];
			subYShifts[imagenumber] = trialsubyshifts[trial];
			defocusshifts[imagenumber] = trialdefocus[trial];
			Debug("Set Shifts for "+boost::lexical_cast<std::string>(imagenumber));
		}
	}
}
*/