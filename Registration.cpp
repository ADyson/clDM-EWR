#include "StdAfx.h"
#include "clKernels3.h"
#include "Registration.h"
#include "clState.h" // gives access to global clContext etc..
#include "PCPCFLib.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


const char* code_clJointHistogram2 = 
"__kernel void clJointHistogram(__global const float* restrict ImageData1, __global const float* restrict ImageData2, __global uint* JointHistogram, int sizeX, int sizeY, float max, float min, float max2, float min2, int xs, int ys) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{ \n"
"		int xp = xid + xs; \n"
"		int yp = yid + ys; \n"
"		if(xp>=sizeX) \n"
"			xp-=sizeX; \n"
"		if(yp>=sizeY) \n"
"			yp-=sizeY; \n"
"		if(xp<0) \n"
"			xp+=sizeX; \n"
"		if(yp<0) \n"
"			yp+=sizeY; \n"
"		float val = ImageData1[xid + sizeX * yid]; \n"
"		float val2 = ImageData2[xp + sizeX * yp]; \n"
"		int bin = floor((val - min)/(max-min) * 255.0f); \n"
"		int bin2 = floor((val2 - min2)/(max2-min2) * 255.0f); \n"
"		if(bin>255) \n"
"			bin=255; \n"
"		if(bin2>255) \n"
"			bin2=255; \n"
"		if(bin<0) \n"
"			bin=0; \n"
"		if(bin2<0) \n"
"			bin2=0; \n"
"		atomic_inc(&JointHistogram[bin+256*bin2]); \n"
"	} \n"
" } \n"
;

// Set upper and lower limit bins to zero before summing the number of none outlier bins to get reduced size....
const char* code_clRSize = 
"__kernel void clRSize(__global uint* JH, int histo, int numhisto, int llimit, int ulimit) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<256&&yid<256) \n"
"	{ \n"
"		for(int i = histo; i <= numhisto; i++) { \n"
"			if(xid<llimit||yid<llimit) { \n"
"				if(xid>ulimit||yid>ulimit) { \n"
"					JH[i*256*256 + xid + 256*yid] = 0; \n"
"				} \n"
"			} \n"
"		} \n"
"	} \n"
"} \n"
;


const char* code_clEntropy = 
"__kernel void clEntropy(__global const uint* JH, __global float* Entropy, int histo, float size) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<256&&yid<256) \n"
"	{ \n"
"		if(JH[histo*256*256 + xid+256*yid]!=0) {\n"
"			Entropy[xid+256*yid]= -((float)JH[histo*256*256 + xid+256*yid]/size) * log(((float)JH[histo*256*256 + xid+256*yid]/size)); \n"
"		} else { \n"
"			Entropy[xid+256*yid]= 0; \n"
"		} \n"
"	} \n"
"} \n"
;


// 400 Joint histograms at once like a badass kernel...
// So Parallel, Much Histogram, Very GPU, WOW
const char* code_clJointHistogramMULTIFix =
"__kernel void clJointHistogramMULTI(__global const float* restrict ImageData1, __global const float* restrict ImageData2, __global int* restrict JointHistograms, int sizeX, int sizeY, float max, float min, float max2, float min2, int xs, int ys,  __local int* restrict tmp_histogram2) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
" // Read bin values into appropriate part of shared memory in parallel \n"
"	for( int lidx =  get_local_id(0); lidx < (get_local_size(0)+19); lidx += get_local_size(0)) \n"
"		for( int lidy =  get_local_id(1); lidy < (get_local_size(1)+19); lidy += get_local_size(1)) \n"
"		{ \n"
"			int xp = lidx + get_group_id(0)*get_local_size(0) + xs; \n"
"			int yp = lidy + get_group_id(1)*get_local_size(1) + ys; \n"
"			//int x2 = lidx + get_group_id(0)*get_local_size(0); \n"
"			//int y2 = lidy + get_group_id(1)*get_local_size(1); \n"
"			if(xp>=sizeX) \n"
"				xp-=sizeX; \n"
"			if(yp>=sizeY) \n"
"				yp-=sizeY; \n"
"			if(xp<0) \n"
"				xp+=sizeX; \n"
"			if(yp<0) \n"
"				yp+=sizeY; \n"
"			//if(x2>=sizeX) // Technically I shouldnt have to wrap these because I dont need to calculate the values of B1 outside normal area but this is easier? \n"
"			//	x2-=sizeX; \n"
"			//if(y2>=sizeY) \n"
"			//	y2-=sizeY; \n"
"			//Want normal B1s from first image, but the second image should also take into account the xs and ys so we can the full range in multiple kernel launches \n"
"			//tmp_histogram[lidx + (get_local_size(0)+19)*lidy] = floor((ImageData1[x2 + sizeX * y2] - min)/(max-min) * 255.0f); \n"
"			tmp_histogram2[lidx + (get_local_size(0)+19)*lidy] = floor((ImageData2[xp + sizeX * yp] - min2)/(max2-min2) * 255.0f); \n"
"		} \n"
"	barrier(CLK_LOCAL_MEM_FENCE);\n"
"	int b1 = floor((ImageData1[xid + sizeX * yid] - min)/(max-min) * 255.0f); \n"
"	if(b1>255) \n"
"		b1=255; \n"
"	if(b1<0) \n"
"		b1=0; \n"
" 	// At this point have a shared memory full of B1 and B2 values to make a 20*20 set of Joint Histograms for every pixel in this workgroup \n"
"	for( int xx = 0; xx < 20; xx++) \n"
"		for( int yy = 0; yy < 20; yy++) \n"
"		{ \n"
"			int histo = xx + 20*yy; \n"
"			int b2 = tmp_histogram2[get_local_id(0) + xx + (get_local_id(1) + yy)*(get_local_size(0)+19)]; \n"
"			if(b2>255) \n"
"				b2=255; \n"
"			if(b2<0) \n"
"				b2=0; \n"
"			atomic_inc(&JointHistograms[(histo*256*256) + b1 + (256 * b2)]); \n"
"		} \n"
" } \n"
;

// 400 Joint histograms at once like a badass kernel...
// So Parallel, Much Histogram, Very GPU, WOW
const char* code_clJointHistogramMULTI =
"__kernel void clJointHistogramMULTI(__global const float* ImageData1, __global const float* ImageData2, __global int* JointHistograms, int sizeX, int sizeY, float max, float min, float max2, float min2, int xs, int ys, __local int* tmp_histogram, __local int* tmp_histogram2) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
" // Read bin values into appropriate part of shared memory in parallel \n"
"	for( int lidx =  get_local_id(0); lidx < (get_local_size(0)+19); lidx += get_local_size(0)) \n"
"		for( int lidy =  get_local_id(1); lidy < (get_local_size(1)+19); lidy += get_local_size(1)) \n"
"		{ \n"
"			int xp = lidx + get_group_id(0)*get_local_size(0) + xs; \n"
"			int yp = lidy + get_group_id(1)*get_local_size(1) + ys; \n"
"			int x2 = lidx + get_group_id(0)*get_local_size(0); \n"
"			int y2 = lidy + get_group_id(1)*get_local_size(1); \n"
"			if(xp>=sizeX) \n"
"				xp-=sizeX; \n"
"			if(yp>=sizeY) \n"
"				yp-=sizeY; \n"
"			if(xp<0) \n"
"				xp+=sizeX; \n"
"			if(yp<0) \n"
"				yp+=sizeY; \n"
"			if(x2>=sizeX) // Technically I shouldnt have to wrap these because I dont need to calculate the values of B1 outside normal area but this is easier? \n"
"				x2-=sizeX; \n"
"			if(y2>=sizeY) \n"
"				y2-=sizeY; \n"
"			//Want normal B1s from first image, but the second image should also take into account the xs and ys so we can the full range in multiple kernel launches \n"
"			tmp_histogram[lidx + (get_local_size(0)+19)*lidy] = floor((ImageData1[x2 + sizeX * y2] - min)/(max-min) * 255.0f); \n"
"			tmp_histogram2[lidx + (get_local_size(0)+19)*lidy] = floor((ImageData2[xp + sizeX * yp] - min2)/(max2-min2) * 255.0f); \n"
"			//tmp_histogram[lidx + (get_local_size(0)+19)*lidy] = 1; \n"
"			//tmp_histogram2[lidx + (get_local_size(0)+19)*lidy] = 1; \n"
"		} \n"
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
" // At this point have a shared memory full of B1 and B2 values to make a 20*20 set of Joint Histograms for every pixel in this workgroup \n"
"	for( int xx = 0; xx < 20; xx++) \n"
"		for( int yy = 0; yy < 20; yy++) \n"
"		{ \n"
"			int histo = xx + 20*yy; \n"
"			int b1 = tmp_histogram[get_local_id(0) + get_local_id(1) * (get_local_size(0)+19)]; \n"
"			int b2 = tmp_histogram2[get_local_id(0) + xx + (get_local_id(1) + yy)*(get_local_size(0)+19)]; \n"
"			if(b1>255) \n"
"				b1=255; \n"
"			if(b2>255) \n"
"				b2=255; \n"
"			if(b1<0) \n"
"				b1=0; \n"
"			if(b2<0) \n"
"				b2=0; \n"
"			atomic_inc(&JointHistograms[(histo*256*256) + b1 + (256 * b2)]); \n"
"		} \n"
" } \n"
;

const char* code_clJointHistogramMULTITest =
"__kernel void clJointHistogramMULTI(__global const float* ImageData1, __global const float* ImageData2, __global uint* JointHistograms, int sizeX, int sizeY, float max, float min, float max2, float min2, int xs, int ys) \n"
"{ \n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"   __local uint  tmp_histogram[(32+9)*(8+9)];\n"
"   __local uint  tmp_histogram2[(32+9)*(8+9)];\n"
" // Read bin values into appropriate part of shared memory in parallel \n"
"	for( int lidx =  get_local_id(0); lidx < (32+9); lidx += get_local_size(0)) \n"
"		for( int lidy =  get_local_id(1); lidy < (8+9); lidy += get_local_size(1)) \n"
"		{ \n"
"			int xp = lidx + get_group_id(0)*get_local_size(0) + xs; \n"
"			int yp = lidy + get_group_id(1)*get_local_size(1) + ys; \n"
"			int x2 = lidx + get_group_id(0)*get_local_size(0); \n"
"			int y2 = lidy + get_group_id(1)*get_local_size(1); \n"
"			if(xp>=sizeX) \n"
"				xp-=sizeX; \n"
"			if(yp>=sizeY) \n"
"				yp-=sizeY; \n"
"			if(xp<0) \n"
"				xp+=sizeX; \n"
"			if(yp<0) \n"
"				yp+=sizeY; \n"
"			if(x2>=sizeX) // Technically I shouldnt have to wrap these because I dont need to calculate the values of B1 outside normal area but this is easier? \n"
"				x2-=sizeX; \n"
"			if(y2>=sizeY) \n"
"				y2-=sizeY; \n"
"			//Want normal B1s from first image, but the second image should also take into account the xs and ys so we can the full range in multiple kernel launches \n"
"			tmp_histogram[lidx + (32+9)*lidy] = floor((ImageData1[x2 + sizeX * y2] - min)/(max-min) * 255.0f); \n"
"			tmp_histogram2[lidx + (32+9)*lidy] = floor((ImageData2[xp + sizeX * yp] - min2)/(max2-min2) * 255.0f); \n"
"			//tmp_histogram[lidx + (32+9)*lidy] = 1; \n"
"			//tmp_histogram2[lidx + (32+9)*lidy] = 1; \n"
"		} \n"
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
" // At this point have a shared memory full of B1 and B2 values to make a 20*20 set of Joint Histograms for every pixel in this workgroup \n"
"	for( int xx = 0; xx < 10; xx++) \n"
"		for( int yy = 0; yy < 10; yy++) \n"
"		{ \n"
"			int histo = xx + 10*yy; \n"
"			uint b1 = tmp_histogram[get_local_id(0) + get_local_id(1) * (32+9)]; \n"
"			uint b2 = tmp_histogram2[get_local_id(0) + xx + (get_local_id(1) + yy)*(32+9)]; \n"
"			if(b1>255) \n"
"				b1=255; \n"
"			if(b2>255) \n"
"				b2=255; \n"
"			if(b1<0) \n"
"				b1=0; \n"
"			if(b2<0) \n"
"				b2=0; \n"
"			atomic_inc(&JointHistograms[(histo*256*256) + b1 + (256 * b2)]); \n"
"		} \n"
" } \n"
;



const char* code_clPartialHistogram2 = 
"__kernel void clPartialHistogram(__global const float* img, __global uint *histogram, int sizeX, int sizeY) \n"
"{ \n"
"	\n"
"    int     local_size = (int)get_local_size(0) *\n"
"                         (int)get_local_size(1);\n"
"    int     image_width = sizeX;\n"
"    int     image_height = sizeY;\n"
"    int     group_indx = (get_group_id(1) * get_num_groups(0)\n"
"                                + get_group_id(0)) * 256;\n"
"    int     x = get_global_id(0);\n"
"    int     y = get_global_id(1);\n"
"   local uint  tmp_histogram[256];\n"
"\n"
"    int     tid = get_local_id(1) * get_local_size(0)\n"
"                                  + get_local_id(0);\n"
"    int     j = 256;\n"
 "   int     indx = 0;\n"
" \n"
"   // clear the local buffer that will generate the partial\n"
"   // histogram\n"
"   do\n"
"    {\n"
"        if (tid < j)\n"
"           tmp_histogram[indx+tid] = 0;\n"
" \n"
"       j -= local_size;\n"
"       indx += local_size;\n"
"    } while (j > 0);\n"
" \n"
"    barrier(CLK_LOCAL_MEM_FENCE);\n"
" \n"
"    if ((x < image_width) && (y < image_height))\n"
"    {\n"
"       float clr = img[x + width*y]; \n"
"      uchar   indx_x;\n"
"      indx_x = convert_uchar_sat(clr * 255.0f);\n"
"     atomic_inc(&tmp_histogram[indx_x]);\n"
"   }\n"
" \n"
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
" \n"
"  // copy the partial histogram to appropriate location in\n"
"  // histogram given by group_indx\n"
" if (local_size >= (256))\n"
"  {\n"
"      if (tid < (256))\n"
"           histogram[group_indx + tid] = tmp_histogram[tid];\n"
"   }\n"
"   else\n"
"   {\n"
"       j = 256;\n"
"       indx = 0;\n"
"       do\n"
"       {\n"
"           if (tid < j)\n"
"               histogram[group_indx + indx + tid] = \n"
"               tmp_histogram[indx + tid];\n"
" \n"
"         j -= local_size;\n"
"        indx += local_size;\n"
"      } while (j > 0);\n"
"  } \n"
" } \n"
;


Registration::Registration(void)
{
	gonext = 1;
	nextdir = -1;

	// Default to false unless changed in load MTFNPS
	gotMTF = false;
	gotNPS = false;
}


Registration::~Registration(void)
{
}

void Registration::Setup(clFourier* FFT,int ref, int num, int wid, int hei, int fullwid, int fullhei, int* xs, int* ys, float*sxs, float* sys, float* dfs, float wavel, float pixscale, float volts, float min, float max)
{
	// Clear values if it has been run once already.
	ImageList.clear();

	OpenCLFFT = FFT;
	referenceimage = ref;
	NumberOfImages = num;
	width = wid;
	height = hei;
	fullwidth = fullwid;
	fullheight = fullhei;
	xShiftVals = xs;
	yShiftVals = ys;
	subXShifts = sxs;
	subYShifts = sys;
	defocusshifts = dfs;
	wavelength = wavel;
	pixelscale = pixscale;
	voltage = volts;
	this->min = min;
	this->max = max;

	gonext = 1;
	nextdir = -1;

	scaleMTFx = 1024.0f/float(wid);
	scaleMTFy = 1024.0f/float(hei);

	//DM::Result(t_to_string(scaleMTFx));
}

void Registration::BuildKernels()
{
	cl_k_Predicted.SetCodeAndName(predictedsource,"clPredicted");
	cl_k_Predicted.BuildKernel();

	cl_k_PCPCF.SetCodeAndName(pcpcfsource2,"clPCPCF");
	cl_k_PCPCF.BuildKernel();

	cl_k_PadCrop.SetCodeAndName(PadCropsource2,"clPadCrop");
	cl_k_PadCrop.BuildKernel();

	cl_k_WaveTransferFunction.SetCodeAndName(wavetransferfunctionsource2,"clWaveTransferFunction");
	cl_k_WaveTransferFunction.BuildKernel();

	cl_k_WaveTransferFunctionMinus.SetCodeAndName(wavetransferfunctionminussource2,"clWaveTransferFunctionMinus");
	cl_k_WaveTransferFunctionMinus.BuildKernel();

	cl_k_WienerW.SetCodeAndName(wienerwsource2,"clWienerW");
	cl_k_WienerW.BuildKernel();

	cl_k_WienerWMinus.SetCodeAndName(wienerwsource2,"clWienerW");
	cl_k_WienerWMinus.BuildKernel();

	cl_k_WienerV.SetCodeAndName(wienervsource2,"clWienerV");
	cl_k_WienerV.BuildKernel();

	cl_k_WienerT.SetCodeAndName(wienertsource2,"clWienerT");
	cl_k_WienerT.BuildKernel();

	cl_k_WienerU.SetCodeAndName(wienerusource2,"clWienerU");
	cl_k_WienerU.BuildKernel();

	cl_k_MakeRestored.SetCodeAndName(makerestoredsource2,"clMakeRestored");
	cl_k_MakeRestored.BuildKernel();

	cl_k_RotScale.SetCodeAndName(RotScalesource2,"clRotScale");
	cl_k_RotScale.BuildKernel();

	cl_k_Q.SetCodeAndName(getQsource2,"clCalculateQ");
	cl_k_Q.BuildKernel();




	cl_k_SumReduction.SetCodeAndName(sumReductionsource2,"clSumReduction");
	cl_k_SumReduction.BuildKernel();

	cl_k_SumReductionFloat.SetCodeAndName(sumReductionsource3,"clSumReductionFloat");
	cl_k_SumReductionFloat.BuildKernel();

	cl_k_SumReductionUint.SetCodeAndName(sumReductionsource4,"clSumReductionUint");
	cl_k_SumReductionUint.BuildKernel();

	cl_k_MinusWavefunction.SetCodeAndName(minuswavefunctionsource2,"clMinusWavefunction");
	cl_k_MinusWavefunction.BuildKernel();

	cl_k_CalculatePCI.SetCodeAndName(getPCIsource2,"clCalculatePCI");
	cl_k_CalculatePCI.BuildKernel();

	cl_k_Abs.SetCodeAndName(abssource2,"clAbs");
	cl_k_Abs.BuildKernel();

	cl_k_HanningWindow.SetCodeAndName(HanningWindowSource	,"clHanningWindow");
	cl_k_HanningWindow.BuildKernel();

	if(gotMTF&&gotNPS)
	{
		cl_k_WaveTransferFunctionMTF.SetCodeAndName(wavetransferfunctionmtfsource2,"clWaveTransferFunctionMTF");
		cl_k_WaveTransferFunctionMTF.BuildKernel();

		cl_k_WaveTransferFunctionMinusMTF.SetCodeAndName(wavetransferfunctionminusmtfsource2,"clWaveTransferFunctionMinusMTF");
		cl_k_WaveTransferFunctionMinusMTF.BuildKernel();

		cl_k_MakeRestoredMTFNPS.SetCodeAndName(makerestoredMTFNPSsource2,"clMakeRestoredMTFNPS");
		cl_k_MakeRestoredMTFNPS.BuildKernel();

		cl_k_Q2.SetCodeAndName(getQsource2NPS,"clCalculateQ");
		cl_k_Q2.BuildKernel();
	}

}

void Registration::RegisterSeries(float* seriesdata, ROIPositions ROIpos, cl_mem &clxFrequencies, cl_mem &clyFrequencies)
{
	// Build all the kernels required for series registration
	BuildKernels();

	int width = ROIpos.Width();
	int height = ROIpos.Height();

	// Calculate work group sizes
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = width;
	globalWorkSize[1] = height;
	globalWorkSize[2] = 1;

	size_t* fullWorkSize = new size_t[3];
	fullWorkSize[0] = fullwidth;
	fullWorkSize[1] = fullheight;
	fullWorkSize[2] = 1;

	// Allocate Group1 Memories
	clMem.SetupGroupOne(width,height,fullwidth,fullheight);

	// Set All Fixed Kernel Arguments.
	SetFixedArgs(ROIpos,clxFrequencies,clyFrequencies,wavelength);

	// Make a new rotation aligned stack based on the expected defocus step.
	std::vector<float> rotscaledseries(fullwidth*fullheight*NumberOfImages);

	// Make new series of rotated and scaled data or just copy data if no rotation magnification is required.
	RotationScale(seriesdata,fullWorkSize,rotscaledseries);

	// Store each of the two images in correlation
	std::vector< std::complex< float > > dataOne( width*height );
	std::vector< std::complex< float > > dataTwo( width*height );

	// Check if reference image is allowed...
	if(referenceimage >=NumberOfImages )
	{
		referenceimage = NumberOfImages - 1;
	}


	referenceimage = options.reference;
	int imageone = referenceimage;
	currentimage = referenceimage;

	// Reference image already registered, it has fixed position.
	ImageList.push_back(referenceimage);

	// Loop over all images to perform registration
	for(int n = 0; n < NumberOfImages - 1; n++)
	{
		if(n < 2)
		{
			// Set Current Image to next image to register.
			IterateImageNumber();

			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = CalculateExpectedDF(imageone,currentimage);

			CopyImageData(imageone,ROIpos.iLeft,ROIpos.iTop,dataOne,rotscaledseries,0,0);
			CopyImageData(currentimage,ROIpos.iLeft,ROIpos.iTop,dataTwo,rotscaledseries,0,0);


			// Insert Magnification Finding Routine
			//MagnificationPCF(30,expecteddifference,dataOne,0,0,globalWorkSize,rotscaledseries,ROIpos.iLeft,ROIpos.iTop,currentimage);

			clEnqueueWriteBuffer( clState::clq->cmdQueue , clMem.clImage1, CL_FALSE, 0, width*height*sizeof(std::complex<float>) , &dataOne[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clState::clq->cmdQueue, clMem.clImage2, CL_TRUE, 0, width*height*sizeof(std::complex<float>) , &dataTwo[ 0 ], 
						0, NULL, NULL );

			Window(clMem.clImage1,width,height);
			Window(clMem.clImage2,width,height);


			OpenCLFFT->Enqueue(clMem.clImage1,clMem.clFFTImage1,CLFFT_FORWARD);
			OpenCLFFT->Enqueue(clMem.clImage2,clMem.clFFTImage2,CLFFT_FORWARD);

			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			
			PhaseCompensatedPCF(numberoftrials,expecteddifference,dataOne,0,0,globalWorkSize);

			// Add to list of registered images.
			ImageList.push_back(currentimage);

		} 
		else 
		{
			// Set Current Image to next image to register.
			IterateImageNumber();
			
			// Find out if it is necessary to pad the images.
			DeterminePadding(ROIpos);
			
			// Add currently registered images to a reconstruction
			AddToReconstruction(rotscaledseries,globalWorkSize,0,0,0);


			cl_k_Predicted.Enqueue(globalWorkSize);
			OpenCLFFT->Enqueue(clMem.clFFTImage2,clMem.clRestored,CLFFT_BACKWARD);
			
			// Need to turn reconstructed ew into fft of image intensity(probably).
			//OpenCLFFT->Enqueue(clMem.clFFTImage1,clMem.clRestored,CLFFT_BACKWARD);
			//cl_k_Abs.Enqueue(globalWorkSize);

			// Hanning Window the first image.
			Window(clMem.clRestored,width,height);

			OpenCLFFT->Enqueue(clMem.clRestored,clMem.clFFTImage1,CLFFT_FORWARD);


			// Determine the amount to preshift the next image by so it is already closer to alignment (i.e start from alignment of next closest image)
			bool prevImageRegistered = false;
			bool nextImageRegistered = false;

			for(int i = 1; i <= ImageList.size(); i++)
			{
								
				if(ImageList[i-1]==currentimage-1)
				{
					prevImageRegistered = true;
					break;
				}
				if(ImageList[i-1]==currentimage+1)
				{
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
			}

			// Cant have it checking an image area that doesnt exist anymore
			if( ROIpos.iLeft + preshiftx < 0)
				preshiftx -=( preshiftx + ROIpos.iLeft );

			if( ROIpos.iLeft + preshiftx + width >= fullwidth)
				preshiftx -=( ROIpos.iLeft + preshiftx + width - fullwidth);

			if( ROIpos.iTop + preshifty < 0)
				preshifty -=( preshifty + ROIpos.iTop );

			if( ROIpos.iTop + preshifty + height >= fullheight)
				preshiftx -=( ROIpos.iTop + preshifty + height - fullheight);


			CopyImageData(currentimage,ROIpos.iLeft,ROIpos.iTop,dataTwo,rotscaledseries,preshiftx,preshifty);

			clEnqueueWriteBuffer( clState::clq->cmdQueue, clMem.clImage2, CL_TRUE, 0, width*height*sizeof(std::complex<float>) , &dataTwo[ 0 ], 0, NULL, NULL );

			// Hanning Window the second image.
			Window(clMem.clImage2,width,height);

			OpenCLFFT->Enqueue(clMem.clImage2,clMem.clFFTImage2,CLFFT_FORWARD);

			float expecteddifference = CalculateExpectedDF(referenceimage,currentimage);

			int numberoftrials = 1;

			if(options.determinefocus)
			{
				numberoftrials = options.steps;
			}

			
			PhaseCompensatedPCF(numberoftrials,expecteddifference,dataOne,preshiftx, preshifty,globalWorkSize);

			ImageList.push_back(currentimage);
		}

	}

	// All images now registered
	Debug("Starting Full Reconstruction");
	
	// Determine final padding requirements.
	DeterminePadding(ROIpos);

	// Add all images to reconstruction
	AddToReconstruction(rotscaledseries,globalWorkSize,0,0,0);

	
	// Output this reconstruction
	OpenCLFFT->Enqueue(clMem.clFFTImage1,clMem.clRestored,CLFFT_BACKWARD);
	DigitalMicrograph::Image PrelimRecon = Utility::PrintCLMemToImagePlusOne(clMem.clRestored,"Restored EW",width,height,clFloat2,clState::clq);

	
	DigitalMicrograph::TagGroup PrelimReconTags = PrelimRecon.GetTagGroup();
	PrelimReconTags.SetTagAsFloat("Reconstruction:Kmax",Aberrations.kmax);
	PrelimReconTags.SetTagAsFloat("Reconstruction:Wavelength",wavelength);
	PrelimRecon.SetDimensionScale(0,pixelscale);
	PrelimRecon.SetDimensionScale(1,pixelscale);

	// Produce a drift corrected series
	MakeDriftCorrectedSeries(rotscaledseries,globalWorkSize);

	// Bit of cleanup
	dataOne.clear();
	dataTwo.clear();
	
	// Don't need these anymore
	cl_k_PCPCF.~clKernel();
	clReleaseMemObject(clMem.clFFTImage2);
	clReleaseMemObject(clMem.clPCPCFResult);
	clReleaseMemObject(clMem.clImage2);

	// Setup work group sizes for reduction kernel and allocate memory.
	int totalSize = width*height;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	// Setup new required memory
	clMem.SetupGroupTwo(width,height,nGroups);

	// Start on PCI to determine focus and astigmatism accurately
	SetupPCI(globalWorkSize);
	
	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	float sumQ = SumReduction(clMem.clQ,globalSizeSum,localSizeSum,nGroups,totalSize);

	// Now trial lots of defocus values
	float DfGuess;
	float A1rGuess;
	float A1iGuess;
	float DfGuess2;
	float A1rGuess2;
	float A1iGuess2;
	float DfGuess3;
	float A1rGuess3;
	float A1iGuess3;

	//Start with non zero to prevent errors....
	float step = (options.enddf-options.startdf)/50;

	PCITrials(sumQ,options.startdf,step,2,2,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess,A1rGuess,A1iGuess,true);
	PCITrials(sumQ,DfGuess-25,1,A1rGuess,A1iGuess,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess2,A1rGuess2,A1iGuess2,false);
	PCITrials(sumQ,DfGuess2-25,1,A1rGuess2,A1iGuess2,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess3,A1rGuess3,A1iGuess3,false);

	// Second PCI Stage Complete
	// Now do restoration again at correct focal/astig value

	AddToReconstruction(rotscaledseries,globalWorkSize,DfGuess2,A1rGuess2,A1iGuess2);

	// inverse fft for display
	OpenCLFFT->Enqueue(clMem.clFFTImage1,clMem.clRestored,CLFFT_BACKWARD);
	DigitalMicrograph::Image Reconstruction = Utility::PrintCLMemToImagePlusOne(clMem.clRestored,"Final Restored EW",width,height,clFloat2,clState::clq);

	DigitalMicrograph::TagGroup ReconTags = Reconstruction.GetTagGroup();
	ReconTags.SetTagAsFloat("Reconstruction:Kmax",Aberrations.kmax);
	ReconTags.SetTagAsFloat("Reconstruction:Wavelength",wavelength);
	Reconstruction.SetDimensionScale(0,pixelscale);
	Reconstruction.SetDimensionScale(1,pixelscale);
	
	// Cleanup
	clMem.CleanUp(gotMTF&&gotNPS);
	rotscaledseries.clear();
	KernelCleanUp();
}

void Registration::MIRegisterSeries(float* seriesdata, ROIPositions ROIpos, cl_mem &clxFrequencies, cl_mem &clyFrequencies)
{
	// Build all the kernels required for series registration
	BuildKernels();

	int width = ROIpos.Width();
	int height = ROIpos.Height();

	// Calculate work group sizes
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = width;
	globalWorkSize[1] = height;
	globalWorkSize[2] = 1;

	size_t* fullWorkSize = new size_t[3];
	fullWorkSize[0] = fullwidth;
	fullWorkSize[1] = fullheight;
	fullWorkSize[2] = 1;


	// Create stack image to store the mutial information maps
	int miSize = 80;
	DigitalMicrograph::Image MIMap = DigitalMicrograph::RealImage("MI Map",4,miSize,miSize,NumberOfImages-1);

	// Display the empty image
	DigitalMicrograph::ImageDocument MIMap_doc = MIMap.GetOrCreateImageDocument();
	MIMap_doc.Show();

	// Allocate Group1 Memories
	clMem.SetupGroupOne(width,height,fullwidth,fullheight);

	// Set All Fixed Kernel Arguments.
	SetFixedArgs(ROIpos,clxFrequencies,clyFrequencies,wavelength);

	// Make a new rotation aligned stack based on the expected defocus step.
	std::vector<float> rotscaledseries(fullwidth*fullheight*NumberOfImages);

	// Make new series of rotated and scaled data or just copy data if no rotation magnification is required.
	RotationScale(seriesdata,fullWorkSize,rotscaledseries);

	// Store each of the two images in correlation
	std::vector< std::complex< float > > dataOne( width*height );
	std::vector< std::complex< float > > dataTwo( width*height );

	// Check if reference image is allowed...
	if(referenceimage >=NumberOfImages )
	{
		referenceimage = NumberOfImages - 1;
	}


	referenceimage = options.reference;
	int imageone = referenceimage;
	currentimage = referenceimage;

	// Reference image already registered, it has fixed position.
	ImageList.push_back(referenceimage);

	// Loop over all images to perform registration
	for(int n = 0; n < NumberOfImages - 1; n++)
	{
		if(n < 2)
		{
			// Set Current Image to next image to register.
			IterateImageNumber();
			// Work out what focus difference is expected between these 2 images.

			// Images should start at most underfocus
			// Positive focus difference if imagetwo > imageone
			float expecteddifference = CalculateExpectedDF(imageone,currentimage);

			if(options.knownfocus)
			{
				if(sgn(currentimage-imageone) > 0)
					expecteddifference = sgn(currentimage-imageone) * options.focuslist[currentimage-1];
				else
					expecteddifference = sgn(currentimage-imageone) * options.focuslist[currentimage];
			}

			CopyImageData(imageone,ROIpos.iLeft,ROIpos.iTop,dataOne,rotscaledseries,0,0);
			CopyImageData(currentimage,ROIpos.iLeft,ROIpos.iTop,dataTwo,rotscaledseries,0,0);
		
			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;


			// with regards to passing the image number to specify what location to put the map. there are N-1 different maps, first image is 0, last image is 19..
			// If image number is less than reference pass image number, if its greater than reference pass image number -1....
			MutualInformationFast(numberoftrials,expecteddifference,dataOne,dataTwo,0,0,globalWorkSize,miSize,MIMap,(currentimage > referenceimage) ? currentimage-1 : currentimage );

			// Add to list of registered images.
			ImageList.push_back(currentimage);

		} 
		else 
		{
			// Determine position we expect the peak to be at based on average shift size on either side of reference image....


			// Set Current Image to next image to register.
			IterateImageNumber();
			
			// Find out if it is necessary to pad the images.
			DeterminePadding(ROIpos);
			
			// Determine the amount to preshift the next image by so it is already closer to alignment (i.e start from alignment of next closest image)
			bool prevImageRegistered = false;
			bool nextImageRegistered = false;

			int otherimage;

			for(int i = 1; i <= ImageList.size(); i++)
			{
								
				if(ImageList[i-1]==currentimage-1)
				{
					prevImageRegistered = true;
					break;
				}
				if(ImageList[i-1]==currentimage+1)
				{
					nextImageRegistered = true;
					break;
				}
			}
			
			int preshiftx = 0;
			int preshifty = 0;
			
			if(prevImageRegistered)
			{
				otherimage = currentimage - 1;
				preshiftx = xShiftVals[currentimage-1];
				preshifty = yShiftVals[currentimage-1];
			}

			if(nextImageRegistered)
			{
				otherimage = currentimage+1;
				preshiftx = xShiftVals[currentimage+1];
				preshifty = yShiftVals[currentimage+1];
			}

			// Cant have it checking an image area that doesnt exist anymore
			if( ROIpos.iLeft + preshiftx < 0)
				preshiftx -=( preshiftx + ROIpos.iLeft );

			if( ROIpos.iLeft + preshiftx + width >= fullwidth)
				preshiftx -=( ROIpos.iLeft + preshiftx + width - fullwidth);

			if( ROIpos.iTop + preshifty < 0)
				preshifty -=( preshifty + ROIpos.iTop );

			if( ROIpos.iTop + preshifty + height >= fullheight)
				preshiftx -=( ROIpos.iTop + preshifty + height - fullheight);

			float expecteddifference = CalculateExpectedDF(otherimage,currentimage);

			// Should select the correct focal displacement from the list...
			if(options.knownfocus)
			{
				if(sgn(currentimage-otherimage) > 0)
					expecteddifference = sgn(currentimage-otherimage) * options.focuslist[currentimage-1];
				else
					expecteddifference = sgn(currentimage-otherimage) * options.focuslist[currentimage];
			}

			// Add currently registered images to a reconstruction

			// NOTE: Want to make reconstruction offset by a defocus to get in plane of next registration attempt...
			// Also want to put result into dataOne for the MI routines
			// Image is added with w at defocusshifts + dfguess... this should make it so the current image would be reconstructed at zero??
			AddToReconstruction(rotscaledseries,globalWorkSize,-(expecteddifference+defocusshifts[otherimage]),0,0);

			// Need to turn reconstructed ew into fft of image intensity(probably).

			cl_k_Predicted.Enqueue(globalWorkSize);

			OpenCLFFT->Enqueue(clMem.clFFTImage2,clMem.clRestored,CLFFT_BACKWARD);
			
			
			//cl_k_Abs.Enqueue(globalWorkSize); // Check this acts on clRestored?


			// Don't need to check reconstructed images..
			//Utility::PrintCLMemToImage(clMem.clRestored,"restored"+Lex(n),width,height,clFloat2,clState::clq);

			clEnqueueReadBuffer(clState::clq->cmdQueue,clMem.clRestored,CL_TRUE,0,width*height*sizeof(std::complex<float>),&dataOne[0],0,0,0);

			//CopyImageData(otherimage,ROIpos.iLeft,ROIpos.iTop,dataOne,rotscaledseries,preshiftx,preshifty);
			CopyImageData(currentimage,ROIpos.iLeft,ROIpos.iTop,dataTwo,rotscaledseries,preshiftx,preshifty);
		
			// Get number of trial steps
			// Assume one if it is not specified
			int numberoftrials = 1;

			Debug("About to MI with reconstructed");
			// with regards to passing the image number to specify what location to put the map. there are N-1 different maps, first image is 0, last image is 19..
			// If image number is less than reference pass image number, if its greater than reference pass image number -1....
			MutualInformationFast(numberoftrials,expecteddifference,dataOne,dataTwo,preshiftx,preshifty,globalWorkSize,miSize,MIMap,(currentimage > referenceimage) ? currentimage-1 : currentimage);

			Debug("MI finished");

			defocusshifts[currentimage]+=defocusshifts[otherimage];

			ImageList.push_back(currentimage);
		}

	}

	// All images now registered
	Debug("Starting Full Reconstruction");

	// Determine final padding requirements.
	DeterminePadding(ROIpos);

	// Add all images to reconstruction
	AddToReconstruction(rotscaledseries,globalWorkSize,0,0,0);
	
	// Output this reconstruction
	OpenCLFFT->Enqueue(clMem.clFFTImage1,clMem.clRestored,CLFFT_BACKWARD);
	DigitalMicrograph::Image PrelimRecon = Utility::PrintCLMemToImagePlusOne(clMem.clRestored,"Restored EW",width,height,clFloat2,clState::clq);

	
	DigitalMicrograph::TagGroup PrelimReconTags = PrelimRecon.GetTagGroup();
	PrelimReconTags.SetTagAsFloat("Reconstruction:Kmax",Aberrations.kmax);
	PrelimReconTags.SetTagAsFloat("Reconstruction:Wavelength",wavelength);
	PrelimRecon.SetDimensionScale(0,pixelscale);
	PrelimRecon.SetDimensionScale(1,pixelscale);

	// Produce a drift corrected series
	MakeDriftCorrectedSeries(rotscaledseries,globalWorkSize);

	// Bit of cleanup
	dataOne.clear();
	dataTwo.clear();
	
	// Don't need these anymore
	cl_k_PCPCF.~clKernel();
	clReleaseMemObject(clMem.clFFTImage2);
	clReleaseMemObject(clMem.clPCPCFResult);
	clReleaseMemObject(clMem.clImage2);

	// Setup work group sizes for reduction kernel and allocate memory.
	int totalSize = width*height;

	// Need to know number of workgroups (wont work for not power 2)
	int nGroups = totalSize / 256;

	// Setup new required memory
	clMem.SetupGroupTwo(width,height,nGroups);

	// Start on PCI to determine focus and astigmatism accurately
	SetupPCI(globalWorkSize);
	
	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = totalSize;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	float sumQ = SumReduction(clMem.clQ,globalSizeSum,localSizeSum,nGroups,totalSize);

	// Now trial lots of defocus values
	float DfGuess;
	float A1rGuess;
	float A1iGuess;
	float DfGuess2;
	float A1rGuess2;
	float A1iGuess2;
	float DfGuess3;
	float A1rGuess3;
	float A1iGuess3;

	//Start with non zero to prevent errors....
	float step = (options.enddf-options.startdf)/50;

	PCITrials(sumQ,options.startdf,step,2,2,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess,A1rGuess,A1iGuess,true);
	PCITrials(sumQ,DfGuess-25,1,A1rGuess,A1iGuess,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess2,A1rGuess2,A1iGuess2,false);
	PCITrials(sumQ,DfGuess2-25,1,A1rGuess2,A1iGuess2,clxFrequencies,clyFrequencies,globalWorkSize,globalSizeSum,localSizeSum,nGroups,totalSize,DfGuess3,A1rGuess3,A1iGuess3,false);

	// Second PCI Stage Complete
	// Now do restoration again at correct focal/astig value

	AddToReconstruction(rotscaledseries,globalWorkSize,DfGuess2,A1rGuess2,A1iGuess2);

	// inverse fft for display
	OpenCLFFT->Enqueue(clMem.clFFTImage1,clMem.clRestored,CLFFT_BACKWARD);
	DigitalMicrograph::Image Reconstruction = Utility::PrintCLMemToImagePlusOne(clMem.clRestored,"Final Restored EW",width,height,clFloat2,clState::clq);

	DigitalMicrograph::TagGroup ReconTags = Reconstruction.GetTagGroup();
	ReconTags.SetTagAsFloat("Reconstruction:Kmax",Aberrations.kmax);
	ReconTags.SetTagAsFloat("Reconstruction:Wavelength",wavelength);
	Reconstruction.SetDimensionScale(0,pixelscale);
	Reconstruction.SetDimensionScale(1,pixelscale);
	
	// Cleanup
	clMem.CleanUp(gotMTF&&gotNPS);
	rotscaledseries.clear();
	KernelCleanUp();
}

void Registration::AddToReconstruction(std::vector<float> &rotscaledseries, size_t* globalWorkSize, float DfGuess, float A1rGuess, float A1iGuess)
{

	int init = 1;
	int noinit = 0;


	// Check imglist
	for(int i = 0 ; i < ImageList.size() ; i++)
	{
		if(i == 0) // Gets set after first image (i=0) is completed.
		{
			cl_k_WienerW.SetArgT(4,init);
			cl_k_WienerWMinus.SetArgT(4,init);
			cl_k_WienerV.SetArgT(5,init);
			cl_k_WienerT.SetArgT(5,init);
			cl_k_WienerU.SetArgT(5,init);
		}
		if(i == 1) // Gets set after first image (i=0) is completed.
		{
			cl_k_WienerW.SetArgT(4,noinit);
			cl_k_WienerWMinus.SetArgT(4,noinit);
			cl_k_WienerV.SetArgT(5,noinit);
			cl_k_WienerT.SetArgT(5,noinit);
			cl_k_WienerU.SetArgT(5,noinit);
		}
				
		// Loop through the registered images extracting the correct area and padding otherwise then add to reconstruction
		cl_k_PadCrop.SetArgT(4,padLeft);
		cl_k_PadCrop.SetArgT(5,padRight);
		cl_k_PadCrop.SetArgT(6,padTop);
		cl_k_PadCrop.SetArgT(7,padBottom);
		cl_k_PadCrop.SetArgT(8,subXShifts[ImageList[i]]);
		cl_k_PadCrop.SetArgT(9,subYShifts[ImageList[i]]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(fullwidth*fullheight);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[ImageList[i]] + DfGuess;

		for(int j = 0 ; j < fullwidth*fullheight ; j++)
		{
			image[j].s[0] = rotscaledseries[ImageList[i]*fullwidth*fullheight + j];
			image[j].s[1] = 0;
		}

		clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.fullImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&image[0],0,NULL,NULL);
		cl_k_PadCrop.Enqueue(globalWorkSize);

		// Check if PadCrop did anything
		//Utility::PrintCLMemToImage(clMem.clImage1,"PadCropTest",width,height,clFloat2,clState::clq);


		// Get FFT of this image
		OpenCLFFT->Enqueue(clMem.clImage1,clMem.clFFTImage1,CLFFT_FORWARD);

		if(gotMTF&&gotNPS)
		{
			cl_k_WaveTransferFunctionMTF.SetArgT(8,A1rGuess);
			cl_k_WaveTransferFunctionMinusMTF.SetArgT(8,A1rGuess);
			cl_k_WaveTransferFunctionMTF.SetArgT(9,A1iGuess);
			cl_k_WaveTransferFunctionMinusMTF.SetArgT(9,A1iGuess);
			cl_k_WaveTransferFunctionMTF.SetArgT(10,defocus);
			cl_k_WaveTransferFunctionMinusMTF.SetArgT(10,defocus);

			cl_k_WaveTransferFunctionMTF.Enqueue(globalWorkSize);
			cl_k_WaveTransferFunctionMinusMTF.Enqueue(globalWorkSize);
		}
		else
		{
			cl_k_WaveTransferFunction.SetArgT(8,A1rGuess);
			cl_k_WaveTransferFunctionMinus.SetArgT(8,A1rGuess);
			cl_k_WaveTransferFunction.SetArgT(9,A1iGuess);
			cl_k_WaveTransferFunctionMinus.SetArgT(9,A1iGuess);
			cl_k_WaveTransferFunction.SetArgT(10,defocus);
			cl_k_WaveTransferFunctionMinus.SetArgT(10,defocus);

			cl_k_WaveTransferFunction.Enqueue(globalWorkSize);
			cl_k_WaveTransferFunctionMinus.Enqueue(globalWorkSize);
		}
	

		// add to W,W-,V,T,U 
		cl_k_WienerW.Enqueue(globalWorkSize);
		cl_k_WienerWMinus.Enqueue(globalWorkSize);
		cl_k_WienerV.Enqueue(globalWorkSize);
		cl_k_WienerT.Enqueue(globalWorkSize);
		cl_k_WienerU.Enqueue(globalWorkSize);
		clFinish(clState::clq->cmdQueue);

		Debug(Lex(i));
	}
	// make restored
	if(!(gotMTF&&gotNPS))
	{
		cl_k_MakeRestored.Enqueue(globalWorkSize);
	//	Debug("Made Restored");
	}
	else
	{
		cl_k_MakeRestoredMTFNPS.Enqueue(globalWorkSize);
		Debug("Made Restored MTFNPS");
		//cl_k_MakeRestored.Enqueue(globalWorkSize);
	//	Debug("Made Restored");
	}


}

void::Registration::MagnificationPCF(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, int preshiftx, int preshifty, size_t* globalWorkSize, std::vector<float> &seriesdata, int iLeft, int iTop, int im )
{
	float minscale = 0.9999f;
	float maxscale = 1.0005f;

	size_t* fullWorkSize = new size_t[3];
	fullWorkSize[0] = fullwidth;
	fullWorkSize[1] = fullheight;
	fullWorkSize[2] = 1;

	// Put into full, image, call rotscale, get from rotscale image.. arg 4 = magcal...
	std::vector<std::complex<float>> copyImage(fullwidth * fullheight) ;
	std::vector<std::complex<float>> returnImage(fullwidth * fullheight) ;

	std::vector<std::complex<float>> data(width * height) ;

	// Copy image into fullimage
	for(int j = 0; j < fullheight; j++)
		for(int i = 0; i < fullwidth; i++)
		{
			copyImage[i+j*fullwidth] = seriesdata[im*fullwidth*fullheight + i + (j)*fullwidth];
		}

	float bestheight = 0.0f;
	float bestscale = minscale;


	for(int trial = 0; trial < numberoftrials; trial++)
	{
		// Rescale one of the images between the min and max scales to test.
		float scale = minscale + trial * (maxscale-minscale)/numberoftrials;	
		// Retrieve and put back in copy image...

		clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.fullImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

		cl_k_RotScale.SetArgT(6,expectedDF);
		cl_k_RotScale.SetArgT(4,scale);
		cl_k_RotScale.Enqueue(fullWorkSize);

		clEnqueueReadBuffer(clState::clq->cmdQueue,clMem.rotScaleImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&returnImage[0],0,NULL,NULL);

		// Now get the roi region and determine shifts with it...

		for(int j = 0; j < height; j++)
			for(int i = 0; i < width; i++)
			{
				data[i+j*width] = returnImage[i + iLeft + preshiftx + (j+iTop + preshifty)*fullwidth];
			}

		clEnqueueWriteBuffer( clState::clq->cmdQueue , clMem.clImage1, CL_FALSE, 0, width*height*sizeof(std::complex<float>) , &dataOne[ 0 ], 
					0, NULL, NULL );
		clEnqueueWriteBuffer( clState::clq->cmdQueue, clMem.clImage2, CL_TRUE, 0, width*height*sizeof(std::complex<float>) , &data[ 0 ], 
					0, NULL, NULL );

			
		OpenCLFFT->Enqueue(clMem.clImage1,clMem.clFFTImage1,CLFFT_FORWARD);
		OpenCLFFT->Enqueue(clMem.clImage2,clMem.clFFTImage2,CLFFT_FORWARD);

		cl_k_PCPCF.SetArgT(7,expectedDF);
		cl_k_PCPCF.Enqueue(globalWorkSize);

		OpenCLFFT->Enqueue(clMem.clPCPCFResult,clMem.clImage1,CLFFT_BACKWARD);

		clEnqueueReadBuffer( clState::clq->cmdQueue, clMem.clImage1, CL_TRUE, 0, width*height*sizeof(std::complex<float>) , &data[ 0 ], 
				0, NULL, NULL );

		float maxheight = 0.0f;
		int maxPosition1 = 0;

		for(int j = 0; j< width * height;j++)
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
		int maxindexr = floor(float((maxPosition1)/(width)))+1;
		int maxindexc = maxPosition1 + 1 - width*floor(float((maxPosition1)/(width)));
	
		float xShift;
		float yShift;
		float subXShift;
		float subYShift;

		// Shift is positive or negative depending on image quadrant it appears in.
		if(maxindexr > floor(height/2 +.5))
			yShift = maxindexr - (height) -1;
		else
			yShift = maxindexr - 1;

		if(maxindexc > floor(width/2 +.5))
			xShift = maxindexc - (width) -1;
		else
			xShift = maxindexc - 1;
				
		// Construct 5*5 region around max to do peakfit (sortof). Probably just CoM. Can do vertex thing also
		//float* peak = new float[25];
		//GetPeak(peak,maxindexc,maxindexr,sizeX,sizeY, data, subXShift, subYShift);
		//delete[] peak;

		//Utility::SetResultWindow(boost::lexical_cast<std::string>(xShift)+"\n");
		//Utility::SetResultWindow(boost::lexical_cast<std::string>(yShift)+"\n");

		// Parabola Vertex Mode
		PCPCFLib::FindVertexParabola(maxPosition1,width,height, data, subXShift, subYShift, maxheight);

		if(maxheight > bestheight) {
				bestheight = maxheight;
				bestscale = scale;
		}

		Debug("Height = "+Lex(maxheight)+" for "+Lex(scale));
	}

	Debug("BEST SCALE IS "+Lex(bestscale));
}

void Registration::PhaseCompensatedPCF(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, int preshiftx, int preshifty, size_t* globalWorkSize)
{
	// Create arrays to hold results of each individual trial run
	int* trialxshifts = new int[numberoftrials];
	int* trialyshifts = new int[numberoftrials];
	float* trialsubxshifts = new float[numberoftrials];
	float* trialsubyshifts = new float[numberoftrials];
	float* trialdefocus = new float[numberoftrials];
	float* peakheights = new float[numberoftrials];
	

	Debug("Max Shift is " + boost::lexical_cast<string>(options.maxdrift)+"\n");

	for(int trial = 0; trial < numberoftrials; trial++)
	{

		//Utility::SetResultWindow("Trial " + Lex(trial+1) + " of " + Lex(numberoftrials));

		// go from expected - searchpercentage to expected + searchpercentage
		float trialdef = CalculateTrialDF(trial,numberoftrials,expectedDF);

		// above formula wont work when number of trials is 1
		if(numberoftrials == 1)
		{
			trialdef = expectedDF;
		}

		cl_k_PCPCF.SetArgT(7,trialdef);
		cl_k_PCPCF.Enqueue(globalWorkSize);

		// Now Inverse FFT
		// Using memory reserved for storing images as it is not needed at this point.
		OpenCLFFT->Enqueue(clMem.clPCPCFResult,clMem.clImage1,CLFFT_BACKWARD);

		clEnqueueReadBuffer( clState::clq->cmdQueue, clMem.clImage1, CL_TRUE, 0, width*height*sizeof(std::complex<float>) , &dataOne[ 0 ], 
					0, NULL, NULL );

		//Utility::SetResultWindow("Registering Image "+Lex(currentimage)+" at defocus "+Lex(trialdef));

		// Temporarily display results of PCPCF's
		//Utility::PrintCLMemToImage(clMem.clImage1,"PCPCF",width,height,clFloat2,clState::clq);

		// Find shift from max height position
		int xShift;
		int yShift;
		float subXShift;
		float subYShift;
		float maxHeight1;

		// Translate linear array index into row and column.
		PCPCFLib::GetShifts(xShift,yShift,subXShift,subYShift,maxHeight1,dataOne,width,height,options.maxdrift);

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
	float bestheight = -FLT_MAX;
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

	if(noalign)
	{
		xShiftVals[currentimage] = 0;
		yShiftVals[currentimage] = 0;
		subXShifts[currentimage] = 0;
		subYShifts[currentimage] = 0;
	}
}

void Registration::MutualInformation(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, std::vector<std::complex<float>> &dataTwo, int preshiftx, int preshifty, size_t* globalWorkSize, int miSize, DigitalMicrograph::Image &MIMap, int imagenumber)
{
	// Create arrays to hold results of each individual trial run
	unsigned int* zeroes = new unsigned int[256*256];
	for (int i = 0 ; i < 256*256; i++)
	{
		zeroes[i] = 0;
	}

	// Get pointer to map data at correct position
	Gatan::PlugIn::ImageDataLocker MILocker(MIMap);
	float* mapdata = (float*)MILocker.get();

	 // std::vector<float> mapdata(miSize*miSize);

	std::vector<unsigned int> hA(256);
	std::vector<unsigned int> hB(256);
	std::vector<unsigned int> hAB(256*256);

	float size = width*height;

	for(int i = 0 ; i < 256; i++)
	{
		hA[i]=0;
		hB[i]=0;

		for(int j = 0 ; j < 256; j++)
		{
			hAB[i+256*j]=0;
		}
	}

	// Build Histogram Kernels..
	clKernel cl_k_JointHistogram;
	cl_k_JointHistogram.SetCodeAndName(code_clJointHistogram2,"clJointHistogram");
	cl_k_JointHistogram.BuildKernel();


	// Calculate work group sizes
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = width;
	globalWorkSize[1] = height;
	globalWorkSize[2] = 1;

	size_t* LocalWorkSize = new size_t[3];
	LocalWorkSize[0] = 32;
	LocalWorkSize[1] = 8;
	LocalWorkSize[2] = 1;

	int nG = ceil((float)width*height/(32*8));

	// Create memory buffers
	cl_mem clImage1		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);
	cl_mem clImage2		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);
	cl_mem clHA		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 256 * nG * sizeof(unsigned int), 0, &clState::status);
	cl_mem clHB		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 256 * nG * sizeof(unsigned int), 0, &clState::status);
	cl_mem clJH		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 256 * 256 * sizeof(unsigned int), 0, &clState::status);


	// Upload the two images
	std::vector<float> im1(width*height);
	std::vector<float> im2(width*height);

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			im1[i+width*j] = dataOne[i+width*j].real();
			im2[i+width*j] = dataTwo[i+width*j].real();
		}


	clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage1,CL_TRUE,0,width*height*sizeof(float),&im1[0],0,NULL,NULL);
	clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage2,CL_TRUE,0,width*height*sizeof(float),&im2[0],0,NULL,NULL);
		
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			int bin1 = floor((im1[i+width*j]-min)/(max-min) * 255.0f);
			int bin2 = floor((im2[i+width*j]-min)/(max-min) * 255.0f);

			if(bin1>255)
				bin1=255;
			if(bin1<0)
				bin1=0;
				
			if(bin2>255)
				bin2=255;
			if(bin2<0)
				bin2=0;

			hA[bin1]++;
			hB[bin2]++;
		}

		// Set Kernel Arguments

		for(int k = 0; k < miSize; k++)
			for(int l = 0; l < miSize; l++)
			{
				int a2 = 0 - miSize/2 + k;
				int b2 = 0 - miSize/2 + l;

				if(!(a2==0&&b2==0))
				{

					//int a2 = preshiftx-miSize/2 + k;
					//int b2 = preshifty - miSize/2 + l;

					clEnqueueWriteBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,0,256*256*sizeof(unsigned int),zeroes,0,NULL,NULL);

					cl_k_JointHistogram.SetArgT(0,clImage1);
					cl_k_JointHistogram.SetArgT(1,clImage2);
					cl_k_JointHistogram.SetArgT(2,clJH);
					cl_k_JointHistogram.SetArgT(3,width);
					cl_k_JointHistogram.SetArgT(4,height);
					cl_k_JointHistogram.SetArgT(5,max);
					cl_k_JointHistogram.SetArgT(6,min);
					cl_k_JointHistogram.SetArgT(7,max);
					cl_k_JointHistogram.SetArgT(8,min);
					cl_k_JointHistogram.SetArgT(9,a2);
					cl_k_JointHistogram.SetArgT(10,b2);

					cl_k_JointHistogram.Enqueue(globalWorkSize);

					// Retrieve Joint Histogram;

					clEnqueueReadBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,0,256*256*sizeof(unsigned int),&hAB[0],0,NULL,NULL);

					//Utility::DisplayArray(&hA[0],"A",256,1);
					//Utility::DisplayArray(&hB[0],"B",256,1);
					//Utility::DisplayArray(&hAB[0],"AB",256,256);


					float eA(0);
					float eB(0);
					float eAB(0);


					for(int i = 0 ; i < 256; i++)
					{
						if(hA[i]!=0)
						{
							eA+= -(hA[i]/size) * log((hA[i]/size));
						}
						if(hB[i]!=0)
						{
							eB+= -(hB[i]/size) * log((hB[i]/size));
						}

						for(int j = 0 ; j < 256; j++)
						{
							if(hAB[i+256*j]!=0)
							{
								eAB+= -(hAB[i+256*j]/(size)) * log(hAB[i+256*j]/(size));
							}
						}
					}

					//if(k==0&&l==0)
					//{
					//	Utility::DisplayArray(&hA[0],"A",256,1);
					//	Utility::DisplayArray(&hB[0],"B",256,1);
					//	Utility::DisplayArray(&hAB[0],"AB",256,256);
					//}

					//DigitalMicrograph::Result(Lex(eA)+"\n");
					//DigitalMicrograph::Result(Lex(eB)+"\n");
				//	DigitalMicrograph::Result(Lex(eAB)+"\n");
					//DigitalMicrograph::Result(Lex(eA + eB - eAB)+"\n");

					float mi = eA + eB - eAB;

					mapdata[imagenumber*miSize*miSize + k+miSize*l] = mi;
				}

			}

		// Fill centre point with average of surrounding
		mapdata[imagenumber*miSize*miSize + miSize/2 + miSize*miSize/2] =
				0.25f * (mapdata[imagenumber*miSize*miSize +miSize/2 +1  + miSize*miSize/2] + mapdata[imagenumber*miSize*miSize +miSize/2 -1 + miSize*miSize/2] + mapdata[imagenumber*miSize*miSize +miSize/2 + miSize*(1+miSize/2)] + mapdata[imagenumber*miSize*miSize +miSize/2 + miSize*(miSize/2-1)]);

		//Utility::DisplayArray(&mapdata[0],"Map",miSize,miSize);

		clReleaseMemObject(clImage1);
		clReleaseMemObject(clImage2);
		clReleaseMemObject(clJH);
		clReleaseMemObject(clHA);
		clReleaseMemObject(clHB);

		int xShift;
		int yShift;
		float subXShift;
		float subYShift;
		float maxHeight1;

		// Translate linear array index into row and column.
		PCPCFLib::GetShiftsMI(xShift,yShift,subXShift,subYShift,maxHeight1,mapdata+imagenumber*miSize*miSize,miSize,miSize,options.maxdrift);

		xShiftVals[currentimage] =  preshiftx + xShift;
		yShiftVals[currentimage] = preshifty + yShift;
		subXShifts[currentimage] = preshiftx + xShift + subXShift;
		subYShifts[currentimage] = preshifty + yShift + subYShift;
		defocusshifts[currentimage] = expectedDF;

		MILocker.~ImageDataLocker();
		MIMap.GetImageDisplay(0).SetDisplayedLayers(imagenumber,imagenumber);

}

void Registration::MutualInformationFast(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, std::vector<std::complex<float>> &dataTwo, int preshiftx, int preshifty, size_t* globalWorkSize, int miSize, DigitalMicrograph::Image &MIMap, int imagenumber)
{

			float averagexshift =0;
			float averageyshift =0;

			for(int i = 1; i < ImageList.size(); i++)
			{
				float xshift = xShiftVals[ImageList[i]]/(ImageList[i]-referenceimage);
				float yshift = yShiftVals[ImageList[i]]/(ImageList[i]-referenceimage);

				averagexshift+=xshift;
				averageyshift+=yshift;
			}

			averagexshift/=(ImageList.size()-1);
			averageyshift/=(ImageList.size()-1);

			Debug("Average x shift is " + Lex(averagexshift));
			Debug("Average y shift is " + Lex(averageyshift));

			int xc = round(averagexshift)*sgn(currentimage-referenceimage) + miSize/2;
			int yc = round(averageyshift)*sgn(currentimage-referenceimage) + miSize/2;


	// Create arrays to hold results of each individual trial run
	int* zeroes = new int[256*256*20*20];
	for (int i = 0 ; i < 256*256*20*20; i++)
	{
		zeroes[i] = 0;
	}

	float* mapdata;
	Gatan::PlugIn::ImageDataLocker MILocker;

	try{
		// Get pointer to map data at correct position
		MILocker = Gatan::PlugIn::ImageDataLocker(MIMap);
		mapdata = (float*)MILocker.get();
	}
	catch(...)
	{
		DigitalMicrograph::Result("Problem getting locker to MI Map\n");
	}
	 // std::vector<float> mapdata(miSize*miSize);

	std::vector<int> hA(256);
	std::vector<int> hB(256);
	std::vector<int> hAB(256*256);

	float size = width*height;

	for(int i = 0 ; i < 256; i++)
	{
		hA[i]=0;
		hB[i]=0;

		for(int j = 0 ; j < 256; j++)
		{
			hAB[i+256*j]=0;
		}
	}

	clKernel cl_k_JointHistogram;
	clKernel cl_k_RSize;
	clKernel cl_k_Entropy;

	// Build Histogram Kernels..
	cl_k_JointHistogram.SetCodeAndName(code_clJointHistogramMULTIFix,"clJointHistogramMULTI");
	cl_k_JointHistogram.BuildKernel();

	cl_k_RSize.SetCodeAndName(code_clRSize,"clRSize");
	cl_k_RSize.BuildKernel();
	
	cl_k_Entropy.SetCodeAndName(code_clEntropy,"clEntropy");
	cl_k_Entropy.BuildKernel();


	// Calculate work group sizes
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = width;
	globalWorkSize[1] = height;
	globalWorkSize[2] = 1;

	size_t* LocalWorkSize = new size_t[3];
	LocalWorkSize[0] = 32;
	LocalWorkSize[1] = 8;
	LocalWorkSize[2] = 1;

	size_t* HistoWorkSize = new size_t[3];
	HistoWorkSize[0] = 256;
	HistoWorkSize[1] = 256;
	HistoWorkSize[2] = 1;

	int nGroups = 256;

	size_t* globalSizeSum = new size_t[3];
	size_t* localSizeSum = new size_t[3];

	globalSizeSum[0] = 256*256;
	globalSizeSum[1] = 1;
	globalSizeSum[2] = 1;
	localSizeSum[0] = 256;
	localSizeSum[1] = 1;
	localSizeSum[2] = 1;

	clMem.clSumOutputFloat= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,nGroups*sizeof(float),0,&clState::status);
	clMem.clSumOutputUint= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,nGroups*sizeof(unsigned int),0,&clState::status);


	int nG = ceil((float)width*height/(32*8));

	// Create memory buffers 
	// MAKE SURE WE USE THESE AND NOT THE CLMEM.clImage1
	cl_mem clImage1		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);
	cl_mem clImage2		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);
	cl_mem clJH		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 256 * 256 * 20 * 20 *  sizeof(int), 0, &clState::status);
	cl_mem clGPUEntropy = clCreateBuffer( clState::context, CL_MEM_READ_WRITE, 256 * 256 *  sizeof(float), 0, &clState::status);


	// Upload the two images
	std::vector<float> im1(width*height);
	std::vector<float> im2(width*height);

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			im1[i+width*j] = dataOne[i+width*j].real();
			im2[i+width*j] = dataTwo[i+width*j].real();
		}


	clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage1,CL_TRUE,0,width*height*sizeof(float),&im1[0],0,NULL,NULL);
	clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage2,CL_TRUE,0,width*height*sizeof(float),&im2[0],0,NULL,NULL);

	// Size not including outlier rejection.
	float reducedsize = 0;
	float reducedsize2 = 0;
	int llimit = 2;
	int ulimit  = 253;
		
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			int bin1 = floor((im1[i+width*j]-min)/(max-min) * 255.0f);
			int bin2 = floor((im2[i+width*j]-min)/(max-min) * 255.0f);

			if(bin1>255)
				bin1=255;
			if(bin1<0)
				bin1=0;
				
			if(bin2>255)
				bin2=255;
			if(bin2<0)
				bin2=0;

			hA[bin1]++;
			hB[bin2]++;

			
			if(!(bin1<llimit||bin1>ulimit))
			{
				reducedsize++;
			}
			if(!(bin2<llimit||bin2>ulimit))
			{
				reducedsize2++;
			}
		}
		
		float eA(0);
		float eB(0);

		for(int i = llimit; i < ulimit; i++)
		{
			if(hA[i]!=0)
			{
				eA+= -((float)hA[i]/reducedsize) * log(((float)hA[i]/reducedsize));
			}

			if(hB[i]!=0)
			{
				eB+= -((float)hB[i]/reducedsize2) * log(((float)hB[i]/reducedsize2));
			}
		}
		
		// Set Kernel Arguments

		for(int k = 0; k < miSize; k+=20)
			for(int l = 0; l < miSize; l+=20)
			{
				int a2 = 0 - miSize/2 + k; // is 0 at 30 or 31???
				int b2 = 0 - miSize/2 + l;

				//if(!(a2==0&&b2==0))
				//{

					//int a2 = preshiftx-miSize/2 + k;
					//int b2 = preshifty - miSize/2 + l;

					clEnqueueWriteBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,0,256*256*20*20*sizeof(int),zeroes,0,NULL,NULL);

					cl_k_JointHistogram.SetArgT(0,clImage1);
					cl_k_JointHistogram.SetArgT(1,clImage2);
					cl_k_JointHistogram.SetArgT(2,clJH);
					cl_k_JointHistogram.SetArgT(3,width);
					cl_k_JointHistogram.SetArgT(4,height);
					cl_k_JointHistogram.SetArgT(5,max);
					cl_k_JointHistogram.SetArgT(6,min);
					cl_k_JointHistogram.SetArgT(7,max);
					cl_k_JointHistogram.SetArgT(8,min);
					cl_k_JointHistogram.SetArgT(9,a2);
					cl_k_JointHistogram.SetArgT(10,b2);
					cl_k_JointHistogram.SetArgLocalMemory(11,51*27,clInt);
					//cl_k_JointHistogram.SetArgLocalMemory(12,51*27,clInt);

					cl_k_JointHistogram.Enqueue3D(globalWorkSize,LocalWorkSize);					

					// Retrieve Joint Histogram;
					//Utility::DisplayArray(&hA[0],"A",256,1);
					//Utility::DisplayArray(&hB[0],"B",256,1);
					//Utility::DisplayArray(&hAB[0],"AB",256,256);

					cl_k_RSize.SetArgT(0,clJH);
					cl_k_RSize.SetArgT(3,llimit);
					cl_k_RSize.SetArgT(4,ulimit);
					cl_k_Entropy.SetArgT(0,clJH);
					cl_k_Entropy.SetArgT(1,clGPUEntropy);

					for(int xx = 0; xx < 20; xx++)
						for(int yy = 0; yy < 20; yy++)
						{
							int histo =xx+20*yy;
							cl_k_RSize.SetArgT(1,histo);
							cl_k_RSize.SetArgT(2,histo);
							cl_k_RSize.Enqueue(HistoWorkSize);

							float rsize = SumReductionUint(clJH,globalSizeSum,localSizeSum,256,256*256,histo*256*256);

							//Debug(Lex(rsize));

							//clEnqueueReadBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,( xx + 20 * yy ) * (256*256*sizeof(int)),256*256*sizeof(int),&hAB[0],0,NULL,NULL);
							
							cl_k_Entropy.SetArgT(2,histo);
							cl_k_Entropy.SetArgT(3,rsize);
							cl_k_Entropy.Enqueue(HistoWorkSize);

							float mi = eA + eB - SumReductionFloat(clGPUEntropy,globalSizeSum,localSizeSum,256,256*256,0);

							mapdata[imagenumber*miSize*miSize + (k+xx)+miSize*(l+yy)] = mi * (80.0f/(80.0f + abs(k+xx-xc)))*(80.0f/(80.0f + abs(l+yy-yc)));

						}
			}




		clReleaseMemObject(clImage1);
		clReleaseMemObject(clImage2);
		clReleaseMemObject(clJH);
		clReleaseMemObject(clGPUEntropy);
		clReleaseMemObject(clMem.clSumOutputFloat);
		clReleaseMemObject(clMem.clSumOutputUint);

		clFinish(clState::clq->cmdQueue);
		
		delete[] zeroes;
		
		int xShift;
		int yShift;
		float subXShift;
		float subYShift;
		float maxHeight1;

		// Translate linear array index into row and column.
		PCPCFLib::GetShiftsMIPreConditioned(xShift,yShift,subXShift,subYShift,maxHeight1,mapdata+imagenumber*miSize*miSize,miSize,miSize,options.maxdrift,averagexshift*sgn(currentimage-referenceimage),averageyshift*sgn(currentimage-referenceimage));

		xShiftVals[currentimage] =  preshiftx + xShift;
		yShiftVals[currentimage] = preshifty + yShift;
		subXShifts[currentimage] = preshiftx + xShift + subXShift;
		subYShifts[currentimage] = preshifty + yShift + subYShift;
		defocusshifts[currentimage] = expectedDF;

		MIMap.GetImageDisplay(0).SetDisplayedLayers(imagenumber,imagenumber);
}

void Registration::MakeDriftCorrectedSeries(std::vector<float> &rotscaledseries, size_t* globalWorkSize)
{
	// Now all should be registered - produce the final restored image and any other relevant data like Shift Graphs, and Drift Corrected Stack
	DigitalMicrograph::Image driftcorrected = DigitalMicrograph::RealImage("Drift Corrected Stack",4,width,height,NumberOfImages);
	Gatan::PlugIn::ImageDataLocker driftlocker = Gatan::PlugIn::ImageDataLocker(driftcorrected);
	float* driftstack = (float*) driftlocker.get();

	std::vector<cl_float2> cropimage(width*height);

	for(int i = 0 ; i < NumberOfImages ; i++)
	{
		cl_k_PadCrop.SetArgT(4,padLeft);
		cl_k_PadCrop.SetArgT(5,padRight);
		cl_k_PadCrop.SetArgT(6,padTop);
		cl_k_PadCrop.SetArgT(7,padBottom);
		cl_k_PadCrop.SetArgT(8,subXShifts[i]);
		cl_k_PadCrop.SetArgT(9,subYShifts[i]);
				
		// Copy correct image into full image..
		std::vector<cl_float2> image(fullwidth*fullheight);

		// Now add image to reconstruction - Make WTF - WTFminus			
		float defocus = defocusshifts[i];
			
		for(int j = 0 ; j < fullwidth*fullheight ; j++)
		{
			image[j].s[0] = rotscaledseries[i*fullwidth*fullheight + j];
			image[j].s[1] = 0;
		}

		// Copy full image to GPU and run crop kernel
		clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.fullImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&image[0],0,NULL,NULL);
		cl_k_PadCrop.Enqueue(globalWorkSize);

		// Copy image data back to host
		clEnqueueReadBuffer(clState::clq->cmdQueue,clMem.clImage1,CL_TRUE,0,width*height*sizeof(cl_float2),&cropimage[0],0,NULL,NULL);
		
		// Copy to drift corrected stack
		for(int k = 0 ; k < width*height ; k++)
		{
			driftstack[i*width*height + k] = cropimage[k].s[0];
		}

		image.clear();
	}

	cropimage.clear();
	driftlocker.~ImageDataLocker();

	// Add tags to drift corrected
	DigitalMicrograph::TagGroup drifttags = driftcorrected.GetTagGroup();
	driftcorrected.SetDimensionScale(0,pixelscale);
	driftcorrected.SetDimensionScale(1,pixelscale);
	driftcorrected.SetDimensionUnitString(0,"nm");
	driftcorrected.SetDimensionUnitString(1,"nm");
	drifttags.SetTagAsString("Focal Series:Adjusted focalstep",Lex(options.focalstep));
	drifttags.SetTagAsString("Microscope Info:Voltage",Lex(voltage));
	drifttags.SetTagAsString("Focal Series:Normalized","True");
	drifttags.SetTagAsString("Focal Series:No Align","True");

	// Display drift corrected
	driftcorrected.GetOrCreateImageDocument().Show();
}

void Registration::SetupPCI(size_t* globalWorkSize)
{

	float newsnr = 1/options.snr;
	
	if(gotMTF&&gotNPS)
	{
		// Set Kernel Arguments
		cl_k_Q2.SetArgT(0,clMem.clW);
		cl_k_Q2.SetArgT(1,clMem.clWminus);
		cl_k_Q2.SetArgT(2,clMem.clV);
		cl_k_Q2.SetArgT(3,clMem.clQ);
		cl_k_Q2.SetArgT(4,width);
		cl_k_Q2.SetArgT(5,height);
		cl_k_Q2.SetArgT(6,options.snr);
		cl_k_Q2.SetArgT(7,clMem.clNPS);
		cl_k_Q2.SetArgT(8,scaleMTFx);
		cl_k_Q2.SetArgT(9,scaleMTFy);
		cl_k_Q2.SetArgT(10,mtflength);
	}
	else
	{
		// Set Kernel Arguments
		cl_k_Q.SetArgT(0,clMem.clW);
		cl_k_Q.SetArgT(1,clMem.clWminus);
		cl_k_Q.SetArgT(2,clMem.clV);
		cl_k_Q.SetArgT(3,clMem.clQ);
		cl_k_Q.SetArgT(4,width);
		cl_k_Q.SetArgT(5,height);
		cl_k_Q.SetArgT(6,newsnr);
	}

	cl_k_MinusWavefunction.SetArgT(0,clMem.clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	cl_k_MinusWavefunction.SetArgT(1,clMem.clRestoredMinus);
	cl_k_MinusWavefunction.SetArgT(2,width);
	cl_k_MinusWavefunction.SetArgT(3,height);

	if(gotMTF&&gotNPS)
	{
		cl_k_Q2.Enqueue(globalWorkSize);
	}
	else
	{
		cl_k_Q.Enqueue(globalWorkSize);
	}
	cl_k_MinusWavefunction.Enqueue(globalWorkSize);


}

void Registration::PCITrials(float sumQ, float startdf, float stepdf, float A1r, float A1i, cl_mem &clxFrequencies, cl_mem &clyFrequencies, size_t* globalWorkSize,
							 size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize, float &DfGuess, float &A1rGuess, float &A1iGuess, bool first)
{

	// Note wasnt adding on previous iterations A1r,A1i values to iterate results....

	// Trial lots of C1 Values
	int numberoftrials = 50;

	float* Fpci = new float[numberoftrials];
	float* Fpcic = new float[numberoftrials];
	float* Fpcis = new float[numberoftrials];
	float* Fpcim = new float[numberoftrials];
	float* Fpci2 = new float[numberoftrials];

	float zero2 = 0.0f;

	float kmax = Aberrations.kmax;

	if(first)
		kmax = 3.5f;


	// Set initial arguments for kernel
	cl_k_CalculatePCI.SetArgT(0,clMem.clQ);
	cl_k_CalculatePCI.SetArgT(1,clMem.clFFTImage1); // Was clRestored but this was inverse FFT for display not frequency space -BUGFIX-
	cl_k_CalculatePCI.SetArgT(2,clMem.clRestoredMinus);
	cl_k_CalculatePCI.SetArgT(3,clxFrequencies);
	cl_k_CalculatePCI.SetArgT(4,clyFrequencies);
	cl_k_CalculatePCI.SetArgT(5,clMem.clPCI);
	cl_k_CalculatePCI.SetArgT(6,clMem.clPCIC);
	cl_k_CalculatePCI.SetArgT(7,clMem.clPCIM);
	cl_k_CalculatePCI.SetArgT(8,clMem.clPCIS);
	cl_k_CalculatePCI.SetArgT(9,width);
	cl_k_CalculatePCI.SetArgT(10,height);
	cl_k_CalculatePCI.SetArgT(12,zero2);
	cl_k_CalculatePCI.SetArgT(13,A1r);
	cl_k_CalculatePCI.SetArgT(14,A1i);
	cl_k_CalculatePCI.SetArgT(15,wavelength);
	cl_k_CalculatePCI.SetArgT(16,kmax);

	for(int trial = 0 ; trial < numberoftrials ; trial++)
	{
		float C1 = startdf+trial*stepdf;
		cl_k_CalculatePCI.SetArgT(11,C1);
		cl_k_CalculatePCI.Enqueue(globalWorkSize);

		float fpcitrial = SumReduction(clMem.clPCI,globalSizeSum,localSizeSum,nGroups,totalSize);
		float fpcictrial = SumReduction(clMem.clPCIC,globalSizeSum,localSizeSum,nGroups,totalSize);
		float fpcistrial = SumReduction(clMem.clPCIS,globalSizeSum,localSizeSum,nGroups,totalSize);
		float fpcimtrial = SumReduction(clMem.clPCIM,globalSizeSum,localSizeSum,nGroups,totalSize);

		// Not sure if this is ok, should do sumQ inside loop.
		Fpci[trial] = fpcitrial/sumQ;
		Fpcic[trial] = fpcictrial/sumQ;
		Fpcis[trial] = fpcistrial/sumQ;
		Fpcim[trial] = fpcimtrial/sumQ;
	}

	// Find the best value.
	// Find max

	float maxheight = 0;
	float C1a = 0;
	int besttrial = 0;

	for(int trial = 0 ; trial < numberoftrials ; trial++)
	{
		if(Fpcim[trial] > maxheight)
		{
			maxheight = Fpcim[trial];
			C1a = startdf + trial*stepdf;
			besttrial = trial;
		}
	}

	float astigaxis = atan2(Fpcis[besttrial],Fpcic[besttrial])/2;

	for(int trial = 0 ; trial < numberoftrials ; trial++)
	{
		Fpci2[trial] = Fpci[trial] - Fpcic[trial] * cos(2*astigaxis) - Fpcis[trial]*sin(2*astigaxis);
	}

	float maxheight2 = 0;
	float C1b = 0;

	for(int trial = 0 ; trial < numberoftrials ; trial++)
	{
		if(Fpci2[trial] > maxheight2)
		{
			maxheight2 = Fpci2[trial];
			C1b = startdf + trial*stepdf;
		}
	}

	DfGuess = (C1a + C1b) /2;
	A1rGuess = A1r + 1.5*cos(2*astigaxis)*(C1a-C1b)/2;
	A1iGuess = A1i + 1.5*sin(2*astigaxis)*(C1a-C1b)/2;

	Utility::SetProgressWindow(
		"defocus: " + Lex(DfGuess),
		"A1 real: " + Lex(A1rGuess),
		"A1 imag: " + Lex(A1iGuess)
		);

	DigitalMicrograph::Image pcigraph = DigitalMicrograph::RealImage("PCI",4,numberoftrials);
	Gatan::PlugIn::ImageDataLocker pciLocker(pcigraph);
	float* pciData = (float*) pciLocker.get();

	for(int trial = 0 ; trial < numberoftrials ; trial++)
	{
		pciData[trial] = Fpci[trial];
	}

	pciLocker.~ImageDataLocker();
	pcigraph.GetOrCreateImageDocument().Show();
	pcigraph.SetDimensionCalibration(0,startdf,stepdf,"nm",0);
	DigitalMicrograph::UpdateImage(pcigraph);


	delete[] Fpci;
	delete[] Fpci2;
	delete[] Fpcic;
	delete[] Fpcim;
	delete[] Fpcis;
}

void Registration::SetFixedArgs(ROIPositions ROIpos,cl_mem &clxFrequencies, cl_mem &clyFrequencies, float wavelength)
{
	int init = 1;
	int noinit = 0;

	int width = ROIpos.Width();
	int height = ROIpos.Height();

	// SNR from options is now just Ps....
	// Get SNR for no NPS by using Pn/Ps for Pn = 1;

	float newsnr = 1/options.snr;

	cl_k_PCPCF.SetArgT(0,clMem.clFFTImage1);
	cl_k_PCPCF.SetArgT(1,clMem.clFFTImage2);
	cl_k_PCPCF.SetArgT(2,clMem.clPCPCFResult);
	cl_k_PCPCF.SetArgT(3,clxFrequencies);
	cl_k_PCPCF.SetArgT(4,clyFrequencies);
	cl_k_PCPCF.SetArgT(5,width);
	cl_k_PCPCF.SetArgT(6,height);
	cl_k_PCPCF.SetArgT(8,wavelength);
	cl_k_PCPCF.SetArgT(9,options.pcpcfkmax);

	cl_k_PadCrop.SetArgT(0,clMem.fullImage);
	cl_k_PadCrop.SetArgT(1,clMem.clImage1);
	cl_k_PadCrop.SetArgT(2,fullwidth); // long not int
	cl_k_PadCrop.SetArgT(3,fullheight); // long not int
	cl_k_PadCrop.SetArgT(10,width);
	cl_k_PadCrop.SetArgT(11,height);
	cl_k_PadCrop.SetArgT(12,ROIpos.iTop);
	cl_k_PadCrop.SetArgT(13,ROIpos.iLeft);

	cl_k_WaveTransferFunction.SetArgT(0,clMem.clw);
	cl_k_WaveTransferFunction.SetArgT(1,clxFrequencies);
	cl_k_WaveTransferFunction.SetArgT(2,clyFrequencies);
	cl_k_WaveTransferFunction.SetArgT(3,width);
	cl_k_WaveTransferFunction.SetArgT(4,height);
	cl_k_WaveTransferFunction.SetArgT(5,wavelength);
	cl_k_WaveTransferFunction.SetArgT(6,Aberrations.beta);
	cl_k_WaveTransferFunction.SetArgT(7,Aberrations.delta);
	cl_k_WaveTransferFunction.SetArgT(8,Aberrations.A1r);
	cl_k_WaveTransferFunction.SetArgT(9,Aberrations.A1i);
	cl_k_WaveTransferFunction.SetArgT(11,Aberrations.Cs);
	cl_k_WaveTransferFunction.SetArgT(12,Aberrations.kmax);

	cl_k_WaveTransferFunctionMinus.SetArgT(0,clMem.clwminus);
	cl_k_WaveTransferFunctionMinus.SetArgT(1,clxFrequencies);
	cl_k_WaveTransferFunctionMinus.SetArgT(2,clyFrequencies);
	cl_k_WaveTransferFunctionMinus.SetArgT(3,width);
	cl_k_WaveTransferFunctionMinus.SetArgT(4,height);
	cl_k_WaveTransferFunctionMinus.SetArgT(5,wavelength);
	cl_k_WaveTransferFunctionMinus.SetArgT(6,Aberrations.beta);
	cl_k_WaveTransferFunctionMinus.SetArgT(7,Aberrations.delta);
	cl_k_WaveTransferFunctionMinus.SetArgT(8,Aberrations.A1r);
	cl_k_WaveTransferFunctionMinus.SetArgT(9,Aberrations.A1i);
	cl_k_WaveTransferFunctionMinus.SetArgT(11,Aberrations.Cs);
	cl_k_WaveTransferFunctionMinus.SetArgT(12,Aberrations.kmax);

	cl_k_WienerW.SetArgT(0,clMem.clW);
	cl_k_WienerW.SetArgT(1,clMem.clw);
	cl_k_WienerW.SetArgT(2,width);
	cl_k_WienerW.SetArgT(3,height);
	cl_k_WienerW.SetArgT(4,init);

	cl_k_WienerWMinus.SetArgT(0,clMem.clWminus);
	cl_k_WienerWMinus.SetArgT(1,clMem.clwminus);
	cl_k_WienerWMinus.SetArgT(2,width);
	cl_k_WienerWMinus.SetArgT(3,height);
	cl_k_WienerWMinus.SetArgT(4,init);

	cl_k_WienerV.SetArgT(0,clMem.clV);
	cl_k_WienerV.SetArgT(1,clMem.clw);
	cl_k_WienerV.SetArgT(2,clMem.clwminus);
	cl_k_WienerV.SetArgT(3,width);
	cl_k_WienerV.SetArgT(4,height);
	cl_k_WienerV.SetArgT(5,init);

	cl_k_WienerT.SetArgT(0,clMem.clT);
	cl_k_WienerT.SetArgT(1,clMem.clw);
	cl_k_WienerT.SetArgT(2,clMem.clFFTImage1);
	cl_k_WienerT.SetArgT(3,width);
	cl_k_WienerT.SetArgT(4,height);
	cl_k_WienerT.SetArgT(5,init);

	cl_k_WienerU.SetArgT(0,clMem.clU);
	cl_k_WienerU.SetArgT(1,clMem.clwminus);
	cl_k_WienerU.SetArgT(2,clMem.clFFTImage1);
	cl_k_WienerU.SetArgT(3,width);
	cl_k_WienerU.SetArgT(4,height);
	cl_k_WienerU.SetArgT(5,init);

	cl_k_MakeRestored.SetArgT(0,clMem.clW);
	cl_k_MakeRestored.SetArgT(1,clMem.clWminus);
	cl_k_MakeRestored.SetArgT(2,clMem.clV);
	cl_k_MakeRestored.SetArgT(3,clMem.clT);
	cl_k_MakeRestored.SetArgT(4,clMem.clU);
	cl_k_MakeRestored.SetArgT(5,clMem.clFFTImage1);
	cl_k_MakeRestored.SetArgT(6,width);
	cl_k_MakeRestored.SetArgT(7,height);
	cl_k_MakeRestored.SetArgT(8,newsnr);

	cl_k_RotScale.SetArgT(0,clMem.fullImage);
	cl_k_RotScale.SetArgT(1,clMem.rotScaleImage);
	cl_k_RotScale.SetArgT(2,fullwidth);
	cl_k_RotScale.SetArgT(3,fullheight);
	cl_k_RotScale.SetArgT(4,options.magcal);
	cl_k_RotScale.SetArgT(5,options.rotcal);

	cl_k_Abs.SetArgT(0,clMem.clRestored);
	cl_k_Abs.SetArgT(1,width);
	cl_k_Abs.SetArgT(2,height);

	if(gotMTF&&gotNPS)
	{
		// Set arguments for MTFNPS only kernels
		cl_k_MakeRestoredMTFNPS.SetArgT(0,clMem.clW);
		cl_k_MakeRestoredMTFNPS.SetArgT(1,clMem.clWminus);
		cl_k_MakeRestoredMTFNPS.SetArgT(2,clMem.clV);
		cl_k_MakeRestoredMTFNPS.SetArgT(3,clMem.clT);
		cl_k_MakeRestoredMTFNPS.SetArgT(4,clMem.clU);
		cl_k_MakeRestoredMTFNPS.SetArgT(5,clMem.clFFTImage1);
		cl_k_MakeRestoredMTFNPS.SetArgT(6,width);
		cl_k_MakeRestoredMTFNPS.SetArgT(7,height);
		cl_k_MakeRestoredMTFNPS.SetArgT(8,clMem.clMTF);
		cl_k_MakeRestoredMTFNPS.SetArgT(9,clMem.clNPS);
		cl_k_MakeRestoredMTFNPS.SetArgT(10,mtflength);
		cl_k_MakeRestoredMTFNPS.SetArgT(11,scaleMTFx);
		cl_k_MakeRestoredMTFNPS.SetArgT(12,scaleMTFy);
		cl_k_MakeRestoredMTFNPS.SetArgT(13,options.snr);

		cl_k_WaveTransferFunctionMTF.SetArgT(0,clMem.clw);
		cl_k_WaveTransferFunctionMTF.SetArgT(1,clxFrequencies);
		cl_k_WaveTransferFunctionMTF.SetArgT(2,clyFrequencies);
		cl_k_WaveTransferFunctionMTF.SetArgT(3,width);
		cl_k_WaveTransferFunctionMTF.SetArgT(4,height);
		cl_k_WaveTransferFunctionMTF.SetArgT(5,wavelength);
		cl_k_WaveTransferFunctionMTF.SetArgT(6,Aberrations.beta);
		cl_k_WaveTransferFunctionMTF.SetArgT(7,Aberrations.delta);
		cl_k_WaveTransferFunctionMTF.SetArgT(8,Aberrations.A1r);
		cl_k_WaveTransferFunctionMTF.SetArgT(9,Aberrations.A1i);
		cl_k_WaveTransferFunctionMTF.SetArgT(11,Aberrations.Cs);
		cl_k_WaveTransferFunctionMTF.SetArgT(12,Aberrations.kmax);
		cl_k_WaveTransferFunctionMTF.SetArgT(13,clMem.clMTF);
		cl_k_WaveTransferFunctionMTF.SetArgT(14,mtflength);
		cl_k_WaveTransferFunctionMTF.SetArgT(15,scaleMTFx);
		cl_k_WaveTransferFunctionMTF.SetArgT(16,scaleMTFy);

		cl_k_WaveTransferFunctionMinusMTF.SetArgT(0,clMem.clwminus);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(1,clxFrequencies);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(2,clyFrequencies);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(3,width);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(4,height);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(5,wavelength);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(6,Aberrations.beta);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(7,Aberrations.delta);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(8,Aberrations.A1r);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(9,Aberrations.A1i);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(11,Aberrations.Cs);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(12,Aberrations.kmax);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(13,clMem.clMTF);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(14,mtflength);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(15,scaleMTFx);
		cl_k_WaveTransferFunctionMinusMTF.SetArgT(16,scaleMTFy);
	}

	cl_k_Predicted.SetArgT(0,clMem.clFFTImage1);
	cl_k_Predicted.SetArgT(1,clMem.clFFTImage2);
	cl_k_Predicted.SetArgT(2,width);
	cl_k_Predicted.SetArgT(3,height);
}

void Registration::RotationScale(float* seriesdata, size_t* fullWorkSize, std::vector<float> &rotscaledseries)
{
	std::vector<std::complex<float>> copyImage(fullwidth * fullheight) ;

	// Rotate and Scale this image by expected defocus..
	// To prevent edges going missing should start at end that should be smallest so all images are only scaled upwards....
	if(options.rotscale)
	{
		int startingimage = 0;

		if(options.magcal > 1)
		{
			startingimage = 0;
		}
		else
		{
			startingimage = NumberOfImages -1;
		}


		for(int im = 0; im < NumberOfImages ; im++ )
		{
			float expecteddifference = (im-startingimage)*options.focalstep;
			cl_k_RotScale.SetArgT(6,expecteddifference);

			// Copy image into fullimage
			for(int j = 0; j < fullheight; j++)
				for(int i = 0; i < fullwidth; i++)
				{
					copyImage[i+j*fullwidth] = seriesdata[im*fullwidth*fullheight + i + (j)*fullwidth];
				}

			clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.fullImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			cl_k_RotScale.Enqueue(fullWorkSize);

			clEnqueueReadBuffer(clState::clq->cmdQueue,clMem.rotScaleImage,CL_TRUE,0,fullwidth*fullheight*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			for(int j = 0; j < fullheight;j++)
				for(int i = 0; i < fullwidth;i++)
				{
					rotscaledseries[im*fullwidth*fullheight + i + j*fullwidth] = copyImage[i + (j)*fullwidth].real();
				}
		}		
	} 
	else
	{
		// Fill new series with old series data.
		for(int im = 0; im < NumberOfImages ; im++ )
			for(int j = 0; j < fullheight;j++)
				for(int i = 0; i < fullwidth;i++)
				{
					rotscaledseries[im*fullwidth*fullheight +i + j*fullwidth] = seriesdata[im*fullwidth*fullheight +i + j*fullwidth];
				}
	}

	copyImage.clear();
}

void Registration::IterateImageNumber()
{
	currentimage += (gonext*nextdir);

	while(! (currentimage < NumberOfImages && currentimage >=0)) // loop over adjacent images.. i.e    10,11,9,12,8,13,7,14,6etc..
	{
		gonext++;
		nextdir *=-1;
		currentimage += (gonext*nextdir);
	}

	// for next time increase difference and change direction
	gonext++;
	nextdir *=-1;

}

float Registration::CalculateExpectedDF(int imageone, int imagetwo)
{
	return (imagetwo-imageone) * options.focalstep;
}

float Registration::CalculateTrialDF(int trialnumber, int numberoftrials, float expectedDF)
{
	return expectedDF - (options.searchpercentage/100.0f)*options.focalstep 
					+ ((float)trialnumber/((float)numberoftrials-1.0f))*2.0f*(options.searchpercentage/100.0f)*options.focalstep; 
}

void Registration::CopyImageData(int imagenumber, int iLeft, int iTop,std::vector<std::complex<float>> &data, std::vector<float> &seriesdata, int preshiftx, int preshifty)
{
	for(int j = 0; j < height; j++)
		for(int i = 0; i < width; i++)
		{
			data[i+j*width] = seriesdata[imagenumber*fullwidth*fullheight + i + iLeft + preshiftx + (j+iTop + preshifty)*fullwidth];
		}
}

void Registration::DeterminePadding(ROIPositions ROIpos)
{
	// Find maximum negative and positive shifts
	// All shifts have to be stored relative to the reference image.
	float maxnegx = 0;
	float maxnegy = 0;
	float maxposx = 0;
	float maxposy = 0;

	// Could just loop over this list of currently registered images instead...
	for(int i = 1; i <= ImageList.size(); i++)
	{
		if(subXShifts[ImageList[i-1]] < maxnegx)
			maxnegx = subXShifts[ImageList[i-1]];
		if(subYShifts[ImageList[i-1]] < maxnegy)
			maxnegy = subYShifts[ImageList[i-1]];
		if(subXShifts[ImageList[i-1]] > maxposx)
			maxposx = subXShifts[ImageList[i-1]];
		if(subYShifts[ImageList[i-1]] > maxposy)
			maxposy = subYShifts[ImageList[i-1]];
	}

	// Determine amount to pad on either direction
	padTop = 0;
	padLeft = 0;
	padRight = 0;
	padBottom = 0;

	if (abs(maxnegx)-ROIpos.iLeft > 0)
		padLeft = ceil(abs(maxnegx)-ROIpos.iLeft);

	if (maxposx-ROIpos.iRight > 0)
		padRight = ceil(maxposx-ROIpos.iRight);

	if (abs(maxnegy)-ROIpos.iTop > 0)
		padTop = ceil(abs(maxnegy)-ROIpos.iTop);

	if (maxposy-ROIpos.iBottom > 0)
		padBottom = ceil(maxposy-ROIpos.iBottom);	
}

float Registration::SumReduction(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize)
{
	// Create host array to store reduction results.
	std::vector< std::complex< float > > sums( nGroups );

	
	cl_k_SumReduction.SetArgT(0,Array);

	// Only really need to do these 3 once...
	cl_k_SumReduction.SetArgT(1,clMem.clSumOutput);
	cl_k_SumReduction.SetArgT(2,totalSize);
	cl_k_SumReduction.SetArgLocalMemory(3,256,clFloat2);

	cl_k_SumReduction.Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clState::clq->cmdQueue, clMem.clSumOutput, CL_TRUE, 0, nGroups*sizeof(std::complex<float>) , &sums[0], 0, NULL, NULL );

	// Find out which numbers to read back
	float sum = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sum += sums[i].real();
	}

	return sum;

}

float Registration::SumReductionFloat(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize, int offset)
{
	// Create host array to store reduction results.
	std::vector< float > sums( nGroups );

	
	cl_k_SumReductionFloat.SetArgT(0,Array);

	// Only really need to do these 3 once...
	cl_k_SumReductionFloat.SetArgT(1,clMem.clSumOutputFloat);
	cl_k_SumReductionFloat.SetArgT(2,totalSize);
	cl_k_SumReductionFloat.SetArgLocalMemory(3,256,clFloat);
	cl_k_SumReductionFloat.SetArgT(4,offset);

	cl_k_SumReductionFloat.Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clState::clq->cmdQueue, clMem.clSumOutputFloat, CL_TRUE, 0, nGroups*sizeof(float) , &sums[0], 0, NULL, NULL );

	// Find out which numbers to read back
	float sum = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sum += sums[i];
	}

	return sum;

}


float Registration::SumReductionUint(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize, int offset)
{
	// Create host array to store reduction results.
	std::vector< unsigned int > sums( nGroups );

	
	cl_k_SumReductionUint.SetArgT(0,Array);

	// Only really need to do these 3 once...
	cl_k_SumReductionUint.SetArgT(1,clMem.clSumOutputUint);
	cl_k_SumReductionUint.SetArgT(2,totalSize);
	cl_k_SumReductionUint.SetArgLocalMemory(3,256,clUInt);
	cl_k_SumReductionUint.SetArgT(4,offset);

	cl_k_SumReductionUint.Enqueue3D(globalSizeSum,localSizeSum);

	// Now copy back 
	clEnqueueReadBuffer( clState::clq->cmdQueue, clMem.clSumOutputUint, CL_TRUE, 0, nGroups*sizeof(unsigned int) , &sums[0], 0, NULL, NULL );

	// Find out which numbers to read back
	float sum = 0;

	for(int i = 0 ; i < nGroups; i++)
	{
		sum += sums[i];
	}

	return sum;

}

void Registration::Window(cl_mem &Image, int width, int height)
{
	// Apply a hanning window to image...
	cl_k_HanningWindow.SetArgT(0,Image);
	cl_k_HanningWindow.SetArgT(1,width);
	cl_k_HanningWindow.SetArgT(2,height);

	// Calculate work group sizes
	size_t* globalWorkSize = new size_t[3];
	globalWorkSize[0] = width;
	globalWorkSize[1] = height;
	globalWorkSize[2] = 1;

	cl_k_HanningWindow.Enqueue(globalWorkSize);

}

void Registration::KernelCleanUp()
{
	cl_k_PadCrop.~clKernel();
	cl_k_WienerW.~clKernel();
	cl_k_WienerWMinus.~clKernel();
	cl_k_WienerV.~clKernel();
	cl_k_WienerT.~clKernel();
	cl_k_WienerU.~clKernel();
	cl_k_WaveTransferFunction.~clKernel();
	cl_k_WaveTransferFunctionMinus.~clKernel();
	cl_k_MakeRestored.~clKernel();

	cl_k_Q.~clKernel();
	cl_k_SumReduction.~clKernel();
	cl_k_MinusWavefunction.~clKernel();
	cl_k_CalculatePCI.~clKernel();

	if(gotMTF&&gotNPS)
	{
		cl_k_WaveTransferFunctionMTF.~clKernel();
		cl_k_WaveTransferFunctionMinusMTF.~clKernel();
		cl_k_MakeRestoredMTFNPS.~clKernel();
	}
}

void Registration::LoadMTFNPS(DigitalMicrograph::Image &MTFImage, DigitalMicrograph::Image &NPSImage)
{
	// Set mtf scale in setup function after we know width and height
	
	float* mtfdata;
	float* npsdata;
	
	// class member now.
	//int mtflength;
	int npslength;

	Gatan::PlugIn::ImageDataLocker mtfLocker;
	Gatan::PlugIn::ImageDataLocker npsLocker;

	gotMTF = true;
	gotNPS = true;
	
	mtfLocker = Gatan::PlugIn::ImageDataLocker(MTFImage);
	mtfdata = (float*) mtfLocker.get();
	mtflength = MTFImage.GetDimensionSize(0);
	
	npsLocker = Gatan::PlugIn::ImageDataLocker(NPSImage);
	npsdata = (float*) npsLocker.get();
	npslength = NPSImage.GetDimensionSize(0);
	

	clMem.SetupMTFNPS(mtflength,npslength);

	// Now load data into clMemory.
	clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.clMTF,CL_TRUE,0,mtflength*sizeof(cl_float),&mtfdata[0],0,NULL,NULL);
	clEnqueueWriteBuffer(clState::clq->cmdQueue,clMem.clNPS,CL_TRUE,0,npslength*sizeof(cl_float),&npsdata[0],0,NULL,NULL);

	Debug("MTF and NPS Loaded");
}


// clRegistrationMemories

void clRegistrationMemories::SetupGroupOne(int width, int height, int fullwidth, int fullheight)
{
	clImage1		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float2), 0, &clState::status);
	clImage2		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float2), 0, &clState::status);
	clFFTImage1		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float2), 0, &clState::status);
	clFFTImage2		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float2), 0, &clState::status);
	clPCPCFResult	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float2), 0, &clState::status);

	fullImage		= clCreateBuffer( clState::context, CL_MEM_READ_WRITE, fullwidth * fullheight * sizeof(cl_float2), 0, &clState::status);
	rotScaleImage	= clCreateBuffer( clState::context, CL_MEM_READ_WRITE, fullwidth * fullheight * sizeof(cl_float2), 0, &clState::status);

	clW				= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( float ), 0, &clState::status);
	clWminus		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( float ), 0, &clState::status);
	clw				= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
	clwminus		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
	clT				= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
	clU				= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
	clV				= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
	clRestored		= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width *height * sizeof( std::complex< float > ), 0, &clState::status);
}

void clRegistrationMemories::SetupGroupTwo(int width, int height, int nGroups)
{
	clRestoredMinus	= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clQ				= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clPCI			= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clPCIC			= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clPCIM			= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clPCIS			= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,width*height*sizeof(std::complex<float>),0,&clState::status);
	clSumOutput		= clCreateBuffer(clState::context,CL_MEM_READ_WRITE,nGroups*sizeof(std::complex<float>),0,&clState::status);

}

void clRegistrationMemories::SetupMTFNPS(int mtflength, int npslength)
{
	clMTF = clCreateBuffer(clState::context,CL_MEM_READ_WRITE,mtflength*sizeof(float),0,&clState::status);
	clNPS = clCreateBuffer(clState::context,CL_MEM_READ_WRITE,npslength*sizeof(float),0,&clState::status);
}

void clRegistrationMemories::CleanUp(bool mtfnps)
{
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


	if(mtfnps)
	{
		clReleaseMemObject(clMTF);
		clReleaseMemObject(clNPS);
	}
}