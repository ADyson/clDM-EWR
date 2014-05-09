#include "stdafx.h"
#include "MutualInformation.h"
#include "Utility.h"
#include "clState.h"
#include "clFourier.h"
#include "PCPCFLib.h"

const char* _RotScalesource2 = 
"__kernel void clRotScale(__global const float2* Input, __global float2* Output, int width, int height, float magscale, float rotscale, float df) \n"
"{ \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	int centrex = width/2;	\n"
"	int centrey = height/2; \n"
"	if(xid < width && yid < height) \n"
"	{	\n"
"		int Index = xid + yid*width; \n"
"		float xp = (float)(xid-centrex); \n"
"		float yp = (float)(yid-centrey); \n"
"		float r = hypot( xp , yp); \n"
"		float theta = atan2( yp, xp); \n"
"		float newx = centrex + r*((magscale*df -df) +1.0f) * cos(theta + df*rotscale); \n"
"		float newy = centrey + r*((magscale*df -df) +1.0f) * sin(theta + df*rotscale); \n"
"		int xs2 = floor(newx); \n"
"		int ys2 = floor(newy); \n"
"		if( (xs2  >= width - 1  || xs2 < 0) || (ys2  >= height -1  || ys2 < 0) ) \n"
"		{ \n"
"			Output[Index].x = 0.0f; \n"
"			Output[Index].y = 0.0f; \n"
"		}  else { \n"
"			float subx = newx - xs2; \n"
"			float suby = newy - ys2; \n"
"			float v1 = (1-subx)*(1-suby)*Input[xs2 + ys2*width].x; \n"
"			float v2 = (subx)*(1-suby)*Input[xs2 + 1 + ys2*width].x; \n"
"			float v3 = (1-subx)*(suby)*Input[xs2 + (ys2 + 1)*width].x; \n"
"			float v4 = (subx)*(suby)*Input[xs2 + 1 + (ys2 + 1)*width].x; \n"
"			Output[Index].x = v1 + v2 + v3 + v4; \n"
"			Output[Index].y = 0.0f; \n"
"		} \n"
"	}	\n"
"}	\n"
;


const char* _pcpcfsource2 = 
"__kernel void clPCPCF(__global const float2* fft1, __global const float2* fft2, __global float2* fftresult, __global const float* CLxFrequencies, __global const float* CLyFrequencies, int sizeX, int sizeY, float focalstep, float wavelength, float pcpcfkmax) \n"
"{		\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		float frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]); \n"
"		float compensation = cos(3.1415926534f * focalstep * wavelength * frequency * frequency); \n"
"		float c1r = fft1[Index].x; \n"
"		float c1i = fft1[Index].y; \n"
"		float c2r = fft2[Index].x; \n"
"		float c2i = fft2[Index].y; \n"
"		float denom2 = sqrt( compensation*compensation   * (c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i) \n" 
"		+ 0.001 * 0.001 + 2 * 0.001 * compensation * (c1r*c2r + c1i*c2i)); \n"
"		fftresult[Index].x = (frequency<=pcpcfkmax)*compensation*(c1r*c2r + c1i*c2i)/denom2; \n"
"		fftresult[Index].y = (frequency<=pcpcfkmax)*compensation*(c2i*c1r - c1i*c2r)/denom2; \n"
"	}	\n"
"}		\n"


;
const char* code_clLogPolar = 
"__kernel void clLogPolar(__global const float* Input, __global float* Output, int width, int height, int polarwidth, int polarheight) \n"
"{ \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	int centrex = width/2;	\n"
"	int centrey = height/2; \n"
"	if(xid < polarwidth && yid < polarheight) \n"
"	{	\n"
"		float ph = polarheight;\n"
"		float pw = polarwidth;\n"
"		float lb = log(1.044f); \n"
"		float lres = (log(0.5f * max(width,height))/lb)/pw; \n"
"		float theta = (2 * yid * 3.141592f ) / ph; \n"
"		float xpos = pow(1.044f,lres * xid+1) * cos(theta) + centrex; \n"
"		float ypos = pow(1.044f,lres * xid+1) * sin(theta) + centrey; \n"
"		int Index = xid + yid*polarwidth; \n"
"		int xpos2 = floor(xpos); \n"
"		int ypos2 = floor(ypos); \n"
"		float interpx = xpos - xpos2; \n"
"		float interpy = ypos - ypos2; \n"
"		if(xpos2>width-2||xpos2 < 0||ypos2>height-2||ypos2 < 0) {\n"
"			Output[Index] = 0;\n"
"		} else {\n"
"			Output[Index] = (1-interpx)*(1-interpy) * Input[(ypos2) * width + xpos2] + (interpx)*(1-interpy) * Input[(ypos2) * width + xpos2 + 1] \n"
"			+ (1-interpx)*(interpy)* Input[(ypos2 + 1) * width + xpos2] + (interpx)*(interpy)* Input[(ypos2 + 1) * width + xpos2 + 1]; \n"
"		} \n"
"	} \n"
"} \n"
;
const char* code_clJointHistogram = 
"__kernel void clJointHistogram(__global const float* ImageData1, __global const float* ImageData2, __global uint* JointHistogram, int sizeX, int sizeY, float max, float min, float max2, float min2, int xs, int ys) \n"
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


const char* code_clPartialHistogram = 
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

const char* cl_code_MultiCorrelation = 
"__kernel void clMultiCorrelation(__global const float2* fft1, __global const float2* fft2, __global float2* PCPCF,__global float2* PCF,__global float2* XCF, __global const float* CLxFrequencies, __global const float* CLyFrequencies, int sizeX, int sizeY, float focalstep, float wavelength, float pcpcfkmax) \n"
"{		\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid<sizeX&&yid<sizeY) \n"
"	{	\n"
"		int Index = xid + yid*sizeX; \n"
"		float frequency = sqrt(CLxFrequencies[xid]*CLxFrequencies[xid] + CLyFrequencies[yid]*CLyFrequencies[yid]); \n"
"		float compensation = cos(3.1415926534f * focalstep * wavelength * frequency * frequency); \n"
"		float c1r = fft1[Index].x; \n"
"		float c1i = fft1[Index].y; \n"
"		float c2r = fft2[Index].x; \n"
"		float c2i = fft2[Index].y; \n"
"		float denom2 = sqrt( compensation*compensation   * (c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i) \n" 
"		+ 0.001 * 0.001 + 2 * 0.001 * compensation * (c1r*c2r + c1i*c2i)); \n"
"		PCPCF[Index].x = (frequency<=pcpcfkmax)*compensation*(c1r*c2r + c1i*c2i)/denom2; \n"
"		PCPCF[Index].y = (frequency<=pcpcfkmax)*compensation*(c2i*c1r - c1i*c2r)/denom2; \n"
"		PCF[Index].x = (frequency<=pcpcfkmax)*(c1r*c2r + c1i*c2i)/sqrt((c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i)); \n"
"		PCF[Index].y = (frequency<=pcpcfkmax)*(c2i*c1r - c1i*c2r)/sqrt((c1r*c1r*c2r*c2r + c1i*c1i*c2i*c2i + c1r*c1r*c2i*c2i + c2r*c2r*c1i*c1i)); \n"
"		XCF[Index].x = (frequency<=pcpcfkmax)*compensation*(c1r*c2r + c1i*c2i); \n"
"		XCF[Index].y = (frequency<=pcpcfkmax)*compensation*(c2i*c1r - c1i*c2r); \n"
"	}	\n"
"}		\n"
;

const char* fftShift = 
"__kernel void clfftShift(__global const float2* Input, __global float2* Output, int width, int height) \n"
"{		\n"
"	//Get the work items ID \n"
"	int xid = get_global_id(0);	\n"
"	int yid = get_global_id(1); \n"
"	if(xid < width && yid < height) \n"
"	{	\n"
"		int Index = xid + yid*width; \n"
"		int Yshift = width*height/2; \n"
"		int Xshift = width/2; \n"
"		int Xmid = width/2; \n"
"		int Ymid = height/2; \n"
"		if( xid < Xmid && yid < Ymid ) \n"
"		{ \n"
"			Output[Index+Yshift+Xshift].x = Input[Index].x; \n"
"			Output[Index+Yshift+Xshift].y = Input[Index].y; \n"	
"		} \n"
"		else if( xid >= Xmid && yid < Ymid ) \n"
"		{ \n"
"			Output[Index+Yshift-Xshift].x = Input[Index].x; \n"
"			Output[Index+Yshift-Xshift].y = Input[Index].y; \n"	
"		} \n"
"		else if( xid < Xmid && yid >= Ymid ) \n"
"		{ \n"
"			Output[Index-Yshift+Xshift].x = Input[Index].x; \n"
"			Output[Index-Yshift+Xshift].y = Input[Index].y; \n"	
"		} \n"
"		else if( xid >= Xmid && yid >= Ymid ) \n"
"		{ \n"
"			Output[Index-Yshift-Xshift].x = Input[Index].x; \n"
"			Output[Index-Yshift-Xshift].y = Input[Index].y; \n"	
"		} \n"	
"	}	\n"
"}	\n"
;

					  


float Gauss(float val, float wid)
{
	return exp(-val*val/(2*wid*wid));
}

float MutualInformation( DM_ImageToken image_token, DM_ImageToken image_token2, long a, long b )
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;
	
	float Au(0);
	float Av(0);
	float Auv(0);
	float Bu(0);
	float Bv(0);
	float Buv(0);

	// check that we have been passed a reference to a byte image..
	// To pass a 'CImage' to a '_DM' function, use the '.get()' method.
	if (!DigitalMicrograph::IsByteImage(image.get()))
	{
		DigitalMicrograph::Result("This routine requires a byte image...");
		return 0;
	}
	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );

	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		// prepare to manipulate the data...
		Gatan::uint8 *data = (Gatan::uint8 *) imageL.get();
		Gatan::uint8 *data2 = (Gatan::uint8 *) imageL2.get();
	
		// manipulate the data...
		for (int i2 = 0; i2 < 100; i2++)
			for (int j2 = 0; j2 < 100; j2++)
			{
				i2+=a;
				j2+=b;

				if(i2<0)
					i2+=width;
				if(i2>=width)
					i2-=width;

				if(j2<0)
					j2+=height;
				if(j2>=height)
					j2-=height;

				Au=0;
				Av=0;
				Auv=0;

				for (int i = 0; i < 100; i++)
					for (int j = 0; j < 100; j++)
					{
						float Gu = Gauss(data[i2+width*j2]-data[i+width*j],100);
						float Gv = Gauss(data2[i2+width*j2]-data2[i+width*j],100);
						Au+=Gu;
						Av+=Gv;
						Auv+=(Gu*Gv);
					}
				Bu+=log(Au/(width*height));
				Bv+=log(Av/(width*height));
				Buv+=log(Auv/(width*height));
			}

		float mi = -(Bu + Bv - Buv)/(width*height);
	
		return mi;
		// tell DigitalMicrograph the image is changed...
	}
	
	PLUG_IN_EXIT
}

float MutualInformation2( DM_ImageToken image_token, DM_ImageToken image_token2, long a, long b )
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;

	DigitalMicrograph::Result(Lex(a)+"\n");
	DigitalMicrograph::Result(Lex(b)+"\n");

	// check that we have been passed a reference to a byte image..
	// To pass a 'CImage' to a '_DM' function, use the '.get()' method.
	if (!DigitalMicrograph::IsByteImage(image.get()))
	{
		DigitalMicrograph::Result("This routine requires a byte image...");
		return 0;
	}
	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );

	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		// prepare to manipulate the data...
		Gatan::uint8 *data = (Gatan::uint8 *) imageL.get();
		Gatan::uint8 *data2 = (Gatan::uint8 *) imageL2.get();
	
		// manipulate the data...

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


		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				int i2 = i+a;
				int j2 = j+b;
			
				if(i2<0)
					i2+=width;
				if(i2>=width)
					i2-=width;

				if(j2<0)
					j2+=height;
				if(j2>=height)
					j2-=height;

				hA[data[i2+width*j2]]++;
				hB[data2[i2+width*j2]]++;
				hAB[(data[i+width*j]+256*data2[i2+width*j2])]++;
			}


		Utility::DisplayArray(&hA[0],"A",256,1);
		Utility::DisplayArray(&hB[0],"B",256,1);
		Utility::DisplayArray(&hAB[0],"AB",256,256);


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

		DigitalMicrograph::Result(Lex(eA)+"\n");
		DigitalMicrograph::Result(Lex(eB)+"\n");
		DigitalMicrograph::Result(Lex(eAB)+"\n");
		DigitalMicrograph::Result(Lex(eA + eB - eAB)+"\n");




		
		return eA + eB - eAB;
		// tell DigitalMicrograph the image is changed...
	}
	
	PLUG_IN_EXIT
}

// This one is for floating point images....
float MutualInformationCL( DM_ImageToken image_token, DM_ImageToken image_token2, long a, long b )
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;

	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );

	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		// Get image mins and maxes for histogram range...

		// prepare to manipulate the data...
		float *data = (float *) imageL.get();
		float *data2 = (float *) imageL2.get();
	

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
		cl_k_JointHistogram.SetCodeAndName(code_clJointHistogram,"clJointHistogram");
		cl_k_JointHistogram.BuildKernel();

		//clKernel cl_k_PartialHistogram;
		//cl_k_PartialHistogram.SetCodeAndName(code_clPartialHistogram,"clPartialHistogram");
		//cl_k_PartialHistogram.BuildKernel();

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

		float min, max, min2, max2;

		ImageDisplayGetContrastLimits(image.GetImageDisplay(0),&min,&max);
		ImageDisplayGetContrastLimits(image2.GetImageDisplay(0),&min2,&max2);


		DigitalMicrograph::Result(Lex(min)+"\n");
		DigitalMicrograph::Result(Lex(max)+"\n");

		// Upload the two images

		clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage1,CL_TRUE,0,width*height*sizeof(float),data,0,NULL,NULL);
		clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage2,CL_TRUE,0,width*height*sizeof(float),data2,0,NULL,NULL);


		imageL.~ImageDataLocker();
		imageL2.~ImageDataLocker();
		// Set Kernel Arguments
		int a2 = a;
		int b2 = b;

		cl_k_JointHistogram.SetArgT(0,clImage1);
		cl_k_JointHistogram.SetArgT(1,clImage2);
		cl_k_JointHistogram.SetArgT(2,clJH);
		cl_k_JointHistogram.SetArgT(3,width);
		cl_k_JointHistogram.SetArgT(4,height);
		cl_k_JointHistogram.SetArgT(5,max);
		cl_k_JointHistogram.SetArgT(6,min);
		cl_k_JointHistogram.SetArgT(7,max2);
		cl_k_JointHistogram.SetArgT(8,min2);
		cl_k_JointHistogram.SetArgT(9,a2);
		cl_k_JointHistogram.SetArgT(10,b2);

		cl_k_JointHistogram.Enqueue(globalWorkSize);

		// Retrieve Joint Histogram;

		clEnqueueReadBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,0,256*256*sizeof(unsigned int),&hAB[0],0,NULL,NULL);

		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				int bin1 = floor((data[i+width*j]-min)/(max-min) * 255.0f);
				int bin2 = floor((data2[i+width*j]-min2)/(max2-min2) * 255.0f);

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



		Utility::DisplayArray(&hA[0],"A",256,1);
		Utility::DisplayArray(&hB[0],"B",256,1);
		Utility::DisplayArray(&hAB[0],"AB",256,256);


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


	//	DigitalMicrograph::Result(Lex(eA)+"\n");
	//	DigitalMicrograph::Result(Lex(eB)+"\n");
	//	DigitalMicrograph::Result(Lex(eAB)+"\n");
		DigitalMicrograph::Result(Lex(eA + eB - eAB)+"\n");

		clReleaseMemObject(clImage1);
		clReleaseMemObject(clImage2);
		clReleaseMemObject(clJH);
		clReleaseMemObject(clHA);
		clReleaseMemObject(clHB);
		//cl_k_JointHistogram.~clKernel();
		//cl_k_PartialHistogram.~clKernel();
		

		
		return eA + eB - eAB;
		// tell DigitalMicrograph the image is changed...


	}
	
	PLUG_IN_EXIT
}


void MutualInformationMap( DM_ImageToken image_token, DM_ImageToken image_token2, long a, long b )
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;

	unsigned int* zeroes = new unsigned int[256*256];

	for (int i = 0 ; i < 256*256; i++)
	{
		zeroes[i] = 0;
	}

	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );

	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		// Create map result image

		DigitalMicrograph::Image map = DigitalMicrograph::RealImage("MAP",4,30,30);
		Gatan::PlugIn::ImageDataLocker imageMAP( map );

		// prepare to manipulate the data...
		float *data = (float *) imageL.get();
		float *data2 = (float *) imageL2.get();
		float *mapdata = (float *) imageMAP.get();

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
		cl_k_JointHistogram.SetCodeAndName(code_clJointHistogram,"clJointHistogram");
		cl_k_JointHistogram.BuildKernel();

		//clKernel cl_k_PartialHistogram;
		//cl_k_PartialHistogram.SetCodeAndName(code_clPartialHistogram,"clPartialHistogram");
		//cl_k_PartialHistogram.BuildKernel();

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

		float min, max, min2, max2;

		ImageDisplayGetContrastLimits(image.GetImageDisplay(0),&min,&max);
		ImageDisplayGetContrastLimits(image2.GetImageDisplay(0),&min2,&max2);


		DigitalMicrograph::Result(Lex(min)+"\n");
		DigitalMicrograph::Result(Lex(max)+"\n");

		// Upload the two images

		clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage1,CL_TRUE,0,width*height*sizeof(float),data,0,NULL,NULL);
		clEnqueueWriteBuffer(clState::clq->cmdQueue,clImage2,CL_TRUE,0,width*height*sizeof(float),data2,0,NULL,NULL);
		
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				int bin1 = floor((data[i+width*j]-min)/(max-min) * 255.0f);
				int bin2 = floor((data2[i+width*j]-min2)/(max2-min2) * 255.0f);

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

		imageL.~ImageDataLocker();
		imageL2.~ImageDataLocker();
		// Set Kernel Arguments

		for(int k = 0; k < 30; k++)
			for(int l = 0; l < 30; l++)
			{


			int a2 = a-15 + k;
			int b2 = b - 15 + l;

			clEnqueueWriteBuffer(clState::clq->cmdQueue,clJH,CL_TRUE,0,256*256*sizeof(unsigned int),zeroes,0,NULL,NULL);

			cl_k_JointHistogram.SetArgT(0,clImage1);
			cl_k_JointHistogram.SetArgT(1,clImage2);
			cl_k_JointHistogram.SetArgT(2,clJH);
			cl_k_JointHistogram.SetArgT(3,width);
			cl_k_JointHistogram.SetArgT(4,height);
			cl_k_JointHistogram.SetArgT(5,max);
			cl_k_JointHistogram.SetArgT(6,min);
			cl_k_JointHistogram.SetArgT(7,max2);
			cl_k_JointHistogram.SetArgT(8,min2);
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



			//DigitalMicrograph::Result(Lex(eA)+"\n");
			//DigitalMicrograph::Result(Lex(eB)+"\n");
		//	DigitalMicrograph::Result(Lex(eAB)+"\n");
			DigitalMicrograph::Result(Lex(eA + eB - eAB)+"\n");

			float mi = eA + eB - eAB;

			mapdata[k+30*l] = mi;

		}

		imageMAP.~ImageDataLocker();
		map.GetOrCreateImageDocument().Show();

		clReleaseMemObject(clImage1);
		clReleaseMemObject(clImage2);
		clReleaseMemObject(clJH);
		clReleaseMemObject(clHA);
		clReleaseMemObject(clHB);
		//cl_k_JointHistogram.~clKernel();
		//cl_k_PartialHistogram.~clKernel();
		
		// tell DigitalMicrograph the image is changed...


	}
	
	PLUG_IN_EXIT
}


void XCFPCFPCPCF( DM_ImageToken image_token, DM_ImageToken image_token2, long df )
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;
	
	// check that we have been passed a reference to a byte image..
	// To pass a 'CImage' to a '_DM' function, use the '.get()' method.
	if (!DigitalMicrograph::IsFloatImage(image.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}

	if (!DigitalMicrograph::IsFloatImage(image2.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}
	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );
	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		float pixelscale = image.GetDimensionScale(0);

		int w2 = 2048;
		int h2 = 2048;

		//while(w2*2 <= width && h2*2 <= height)
		//{
		//	w2*=2;
		//	h2*=2;
		//}

		
		// prepare to manipulate the data...
		float *data = (float*) imageL.get();
		float *data2 = (float*) imageL2.get();
	
		// copy to a couple of buffers

		// Build Histogram Kernels..
		clKernel cl_k_MultiCorrelation;
		cl_k_MultiCorrelation.SetCodeAndName(cl_code_MultiCorrelation,"clMultiCorrelation");
		cl_k_MultiCorrelation.BuildKernel();


		// Calculate work group sizes
		size_t* globalWorkSize = new size_t[3];
		globalWorkSize[0] = w2;
		globalWorkSize[1] = h2;
		globalWorkSize[2] = 1;

		std::vector<std::complex<float>> h_Image1(w2*h2);
		std::vector<std::complex<float>> h_Image2(w2*h2);

		for(int i = 0; i < w2; i++)
			for(int j = 0; j < h2; j++)
			{
				h_Image1[i + w2*j] = data[i + width*j];
				h_Image2[i + w2*j] = data2[i + width*j];
			}
	
		// Create memory buffers
		cl_mem clImage1	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * h2 * sizeof(cl_float2), 0, &clState::status);
		cl_mem clImage2	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * h2 * sizeof(cl_float2), 0, &clState::status);

		cl_mem clXCF	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * h2 * sizeof(cl_float2), 0, &clState::status);
		cl_mem clPCF	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * h2 * sizeof(cl_float2), 0, &clState::status);
		cl_mem clPCPCF	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * h2 * sizeof(cl_float2), 0, &clState::status);

		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage1, CL_TRUE, 0, w2*h2*sizeof(cl_float2), &h_Image1[0], 0, 0, 0);
		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage2, CL_TRUE, 0, w2*h2*sizeof(cl_float2), &h_Image2[0], 0, 0, 0);


			// Setup two frequency arrays for this size.
		float* xFrequencies = new float[w2];
		float* yFrequencies = new float[h2];

		int midX = ceil(float(w2)/2);
		int midY = ceil(float(h2)/2);
	
		for(int i = 1 ; i <= w2 ; i++)
		{
			if(i <= midX)
				xFrequencies[i-1] = (i - 1)/(pixelscale * w2);
			else xFrequencies[i-1] = (i - 1 - w2)/(pixelscale * w2);
		}

		for(int i = 1 ; i <= h2 ; i++)
		{
			if(i <= midY)
				yFrequencies[i-1] = (i - 1)/(pixelscale * h2);
			else yFrequencies[i-1] = (i - 1 - h2)/(pixelscale * h2);
		}
		// Upload frequencies to GPU!
		cl_mem clxFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, w2 * sizeof(cl_float), 0, &clState::status);
		cl_mem clyFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, h2 * sizeof(cl_float), 0, &clState::status);

		clEnqueueWriteBuffer( clState::clq->cmdQueue, clxFrequencies, CL_FALSE, 0, w2*sizeof(float), &xFrequencies[0], 
				0, NULL, NULL );
		clEnqueueWriteBuffer( clState::clq->cmdQueue, clyFrequencies, CL_TRUE, 0, h2*sizeof(float), &yFrequencies[0], 
				0, NULL, NULL );

		// Fourier them...

		DigitalMicrograph::TagGroup imagetags = image.GetTagGroup();

		float voltage;
		imagetags.GetTagAsFloat("Microscope Info:Voltage",&voltage);

		float e = 1.6e-19;
		float wavelength = 6.63e-34*3e+8*1e+9 / sqrt((e*voltage*(2*9.11e-31*9e+16 + e*voltage)));

		float pcpcf =4.0f;

		float df2 = df;
		clFourier* FFT = new clFourier(clState::context, clState::clq);
		FFT->Setup(w2,h2);
	
		FFT->Enqueue(clImage1,clImage1,CLFFT_FORWARD);
		FFT->Enqueue(clImage2,clImage2,CLFFT_FORWARD);

		cl_k_MultiCorrelation.SetArgT(0,clImage1);
		cl_k_MultiCorrelation.SetArgT(1,clImage2);
		cl_k_MultiCorrelation.SetArgT(2,clPCPCF);
		cl_k_MultiCorrelation.SetArgT(3,clPCF);
		cl_k_MultiCorrelation.SetArgT(4,clXCF);
		cl_k_MultiCorrelation.SetArgT(5,clxFrequencies);
		cl_k_MultiCorrelation.SetArgT(6,clyFrequencies);
		cl_k_MultiCorrelation.SetArgT(7,w2);
		cl_k_MultiCorrelation.SetArgT(8,h2);
		cl_k_MultiCorrelation.SetArgT(9,df2);
		cl_k_MultiCorrelation.SetArgT(10,wavelength);
		cl_k_MultiCorrelation.SetArgT(11,pcpcf);

		cl_k_MultiCorrelation.Enqueue(globalWorkSize);

		FFT->Enqueue(clPCPCF,clPCPCF,CLFFT_BACKWARD);
		FFT->Enqueue(clPCF,clPCF,CLFFT_BACKWARD);
		FFT->Enqueue(clXCF,clXCF,CLFFT_BACKWARD);

		clKernel cl_k_fftShift;
		cl_k_fftShift.SetCodeAndName(fftShift,"clfftShift");
		cl_k_fftShift.BuildKernel();

		cl_k_fftShift.SetArgT(0,clPCPCF);
		cl_k_fftShift.SetArgT(1,clImage1);
		cl_k_fftShift.SetArgT(2,w2);
		cl_k_fftShift.SetArgT(3,h2);
		cl_k_fftShift.Enqueue(globalWorkSize);


		Utility::PrintCLMemToImage(clImage1,"PCPCF",w2,h2,clFloat2,clState::clq);

		cl_k_fftShift.SetArgT(0,clPCF);
		cl_k_fftShift.SetArgT(1,clImage2);
		cl_k_fftShift.Enqueue(globalWorkSize);


		Utility::PrintCLMemToImage(clImage2,"PCF",w2,h2,clFloat2,clState::clq);

		cl_k_fftShift.SetArgT(0,clXCF);
		cl_k_fftShift.SetArgT(1,clPCPCF);
		cl_k_fftShift.Enqueue(globalWorkSize);


		Utility::PrintCLMemToImage(clPCPCF,"XCF",w2,h2,clFloat2,clState::clq);
	}
	
	PLUG_IN_EXIT
}


void Magnification( DM_ImageToken image_token, DM_ImageToken image_token2)
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;
	
	// check that we have been passed a reference to a byte image..
	// To pass a 'CImage' to a '_DM' function, use the '.get()' method.
	if (!DigitalMicrograph::IsFloatImage(image.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}

	if (!DigitalMicrograph::IsFloatImage(image2.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}
	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );
	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		float pixelscale = image.GetDimensionScale(0);
	
		// prepare to manipulate the data...
		float *data = (float*) imageL.get();
		float *data2 = (float*) imageL2.get();
	
		int ph = 1024;
		int pw = 1024;

		// Build Histogram Kernels..
		clKernel cl_k_LogPolar;
		cl_k_LogPolar.SetCodeAndName(code_clLogPolar,"clLogPolar");
		cl_k_LogPolar.BuildKernel();


		// Calculate work group sizes
		size_t* globalWorkSize = new size_t[3];
		globalWorkSize[0] = 1024;
		globalWorkSize[1] = 1024;
		globalWorkSize[2] = 1;

		std::vector<float> h_Image1(width*height);
		std::vector<float> h_Image2(width*height);

		for(int i = 0; i < width; i++)
			for(int j = 0; j < height; j++)
			{
				h_Image1[i + width*j] = data[i + width*j];
				h_Image2[i + width*j] = data2[i + width*j];
			}
	
		// Create memory buffers
		cl_mem clImage1	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);
		cl_mem clImage2	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width * height * sizeof(cl_float), 0, &clState::status);

		cl_mem clLP1	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 1024 * 1024 * sizeof(cl_float), 0, &clState::status);
		cl_mem clLP2	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, 1024 * 1024 * sizeof(cl_float), 0, &clState::status);

		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage1, CL_TRUE, 0, width * height *sizeof(cl_float), &h_Image1[0], 0, 0, 0);
		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage2, CL_TRUE, 0, width * height *sizeof(cl_float), &h_Image2[0], 0, 0, 0);

		cl_k_LogPolar.SetArgT(0,clImage1);
		cl_k_LogPolar.SetArgT(1,clLP1);
		cl_k_LogPolar.SetArgT(2,width);
		cl_k_LogPolar.SetArgT(3,height);
		cl_k_LogPolar.SetArgT(4,pw);
		cl_k_LogPolar.SetArgT(5,ph);

		cl_k_LogPolar.Enqueue(globalWorkSize);

		cl_k_LogPolar.SetArgT(0,clImage2);
		cl_k_LogPolar.SetArgT(1,clLP2);

		cl_k_LogPolar.Enqueue(globalWorkSize);

				// Now should have two log polar transformed images which we can cross correlate as usual and then work back to find scale and rotation.
		Utility::PrintCLMemToImage(clImage1,"1",width,height,clFloat,clState::clq);
		Utility::PrintCLMemToImage(clLP1,"LP1",pw,ph,clFloat,clState::clq);
		Utility::PrintCLMemToImage(clLP2,"LP2",pw,ph,clFloat,clState::clq);

	}
	
	PLUG_IN_EXIT
}


void MagnificationTest( DM_ImageToken image_token, DM_ImageToken image_token2, float expectedDF, long numberoftrials, float minscale, float maxscale)
{
	// PLUG_IN_ENTRY and PLUG_IN_EXIT are required in any script function to handle
	// C++ exceptions properly.
	PLUG_IN_ENTRY
	
	DigitalMicrograph::Image image( image_token );
	DigitalMicrograph::Image image2( image_token2 );

	long width, height, imageSize;
	
	// check that we have been passed a reference to a byte image..
	// To pass a 'CImage' to a '_DM' function, use the '.get()' method.
	if (!DigitalMicrograph::IsFloatImage(image.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}

	if (!DigitalMicrograph::IsFloatImage(image2.get()))
	{
		DigitalMicrograph::Result("This routine requires a float image...");
		return;
	}
	
	// get the size of the image and a handle..
	DigitalMicrograph::GetSize( image.get(), &width, &height );
	{
		Gatan::PlugIn::ImageDataLocker imageL( image );
		Gatan::PlugIn::ImageDataLocker imageL2( image2 );

		float pixelscale = image.GetDimensionScale(0);
	
		// prepare to manipulate the data...
		float *data1 = (float*) imageL.get();
		float *data2 = (float*) imageL2.get();
	
		// Get the ROI from the first image

		int numROI;
		DigitalMicrograph::ImageDisplay imDisp;
		imDisp = image.GetImageDisplay(0);
	
		numROI = imDisp.CountROIs();

		DigitalMicrograph::ROI regionOfInterest;

		if(numROI==0)
		{
			// Create image display
			regionOfInterest = DigitalMicrograph::NewROI();
			regionOfInterest.SetRectangle(0,0,height,width);
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
			DigitalMicrograph::OpenAndSetProgressWindow("ROI is not a rectangle","","");
			return;
		}
	

		regionOfInterest.GetRectangle(&top,&left,&bottom,&right);

		// Do Some Checks
		int iTop = boost::numeric_cast<int>(top);
		int iLeft = boost::numeric_cast<int>(left);
		int iBottom = boost::numeric_cast<int>(bottom);
		int iRight = boost::numeric_cast<int>(right);

		int sizeX = iRight - iLeft;
		int sizeY = iBottom - iTop;

		std::vector<std::complex<float>> h_Image1(sizeX*sizeY);
		std::vector<std::complex<float>> h_Image2(sizeX*sizeY);
		std::vector<std::complex<float>> copyImage(width * height) ;

		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
			{
				h_Image1[i + sizeX*j] = data1[iLeft + i + width*(j+iTop)];
				//h_Image2[i + sizeX*j] = data2[iLeft + i + width*(j+iTop)];
			}

		for(int i = 0; i < width; i++)
			for(int j = 0; j < height; j++)
			{
				copyImage[i + width*j] = data2[i + width*j];
			}
	
		// Create memory buffers
		cl_mem clImage1	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX*sizeY* sizeof(cl_float2), 0, &clState::status);
		cl_mem clImage2	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX* sizeY * sizeof(cl_float2), 0, &clState::status);

		cl_mem clFFTImage1	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX*sizeY* sizeof(cl_float2), 0, &clState::status);
		cl_mem clFFTImage2	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX* sizeY * sizeof(cl_float2), 0, &clState::status);


		cl_mem clfullImage	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width*height* sizeof(cl_float2), 0, &clState::status);
		cl_mem clrotScaleImage	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, width*height * sizeof(cl_float2), 0, &clState::status);

		cl_mem clPCPCFResult	= clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX*sizeY* sizeof(cl_float2), 0, &clState::status);


		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX * sizeY *sizeof(cl_float2), &h_Image1[0], 0, 0, 0);
		clEnqueueWriteBuffer(clState::clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX * sizeY *sizeof(cl_float2), &h_Image2[0], 0, 0, 0);


		clFourier* OpenCLFFT = new clFourier(clState::context, clState::clq);
		OpenCLFFT->Setup(sizeX,sizeY);

		DigitalMicrograph::TagGroup imagetags = image.GetTagGroup();

		float voltage;
		imagetags.GetTagAsFloat("Microscope Info:Voltage",&voltage);

			// Setup two frequency arrays for this size.
		float* xFrequencies = new float[sizeX];
		float* yFrequencies = new float[sizeY];

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


				// Upload frequencies to GPU!
		cl_mem clxFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX * sizeof(cl_float), 0, &clState::status);
		cl_mem clyFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeY * sizeof(cl_float), 0, &clState::status);

		clEnqueueWriteBuffer( clState::clq->cmdQueue, clxFrequencies, CL_FALSE, 0, sizeX*sizeof(float), &xFrequencies[0], 
					0, NULL, NULL );
		clEnqueueWriteBuffer( clState::clq->cmdQueue, clyFrequencies, CL_TRUE, 0, sizeY*sizeof(float), &yFrequencies[0], 
					0, NULL, NULL );
		
		// What needs fixing
		// Need to make pcpcf and rotscale kernels.
		// Need to make a float2 copy of ROI pos from both images..
		// Need frequencies for the pcpcpcpcpcf kernel.

		size_t* fullWorkSize = new size_t[3];
		fullWorkSize[0] = width;
		fullWorkSize[1] = height;
		fullWorkSize[2] = 1;

		size_t* globalWorkSize = new size_t[3];
		globalWorkSize[0] = sizeX;
		globalWorkSize[1] = sizeY;
		globalWorkSize[2] = 1;

		// Put into full, image, call rotscale, get from rotscale image.. arg 4 = magcal...
		std::vector<std::complex<float>> returnImage(width * height) ;

		std::vector<std::complex<float>> data(sizeX * sizeY) ;

		float bestheight = 0.0f;
		float bestscale = minscale;

		float four = 4.0f;

		clKernel cl_k_PCPCF;
		cl_k_PCPCF.SetCodeAndName(_pcpcfsource2,"clPCPCF");
		cl_k_PCPCF.BuildKernel();

		cl_k_PCPCF.SetArgT(0,clFFTImage1);
		cl_k_PCPCF.SetArgT(1,clFFTImage2);
		cl_k_PCPCF.SetArgT(2,clPCPCFResult);
		cl_k_PCPCF.SetArgT(3,clxFrequencies);
		cl_k_PCPCF.SetArgT(4,clyFrequencies);
		cl_k_PCPCF.SetArgT(5,sizeX);
		cl_k_PCPCF.SetArgT(6,sizeY);
		cl_k_PCPCF.SetArgT(8,wavelength);
		cl_k_PCPCF.SetArgT(9,four);

		clKernel cl_k_RotScale;
		cl_k_RotScale.SetCodeAndName(_RotScalesource2,"clRotScale");
		cl_k_RotScale.BuildKernel();

		float zero = 0.0f;
		cl_k_RotScale.SetArgT(0,clfullImage);
		cl_k_RotScale.SetArgT(1,clrotScaleImage);
		cl_k_RotScale.SetArgT(2,width);
		cl_k_RotScale.SetArgT(3,height);
		cl_k_RotScale.SetArgT(5,zero);

		// Make graph to store results in...

		DigitalMicrograph::Image ResultGraph = DigitalMicrograph::RealImage("Results",4,numberoftrials);
		Gatan::PlugIn::ImageDataLocker ResultLocker = Gatan::PlugIn::ImageDataLocker(ResultGraph);
		float* results = (float*) ResultLocker.get();


		for(int trial = 0; trial < numberoftrials; trial++)
		{
			// Rescale one of the images between the min and max scales to test.
			float scale = minscale + trial * (maxscale-minscale)/numberoftrials;	
			// Retrieve and put back in copy image...

			clEnqueueWriteBuffer(clState::clq->cmdQueue,clfullImage,CL_TRUE,0,width*height*sizeof(cl_float2),&copyImage[0],0,NULL,NULL);

			cl_k_RotScale.SetArgT(6,expectedDF);
			cl_k_RotScale.SetArgT(4,scale);
			cl_k_RotScale.Enqueue(fullWorkSize);

			clEnqueueReadBuffer(clState::clq->cmdQueue,clrotScaleImage,CL_TRUE,0,width*height*sizeof(cl_float2),&returnImage[0],0,NULL,NULL);

			// Now get the roi region and determine shifts with it...

			for(int j = 0; j < sizeY; j++)
				for(int i = 0; i < sizeX; i++)
				{
					data[i+j*sizeX] = returnImage[i + iLeft + (j+iTop)*width];
				}

			clEnqueueWriteBuffer( clState::clq->cmdQueue , clImage1, CL_FALSE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &h_Image1[ 0 ], 
						0, NULL, NULL );
			clEnqueueWriteBuffer( clState::clq->cmdQueue, clImage2, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &data[ 0 ], 
						0, NULL, NULL );

			
			OpenCLFFT->Enqueue(clImage1,clFFTImage1,CLFFT_FORWARD);
			OpenCLFFT->Enqueue(clImage2,clFFTImage2,CLFFT_FORWARD);

		
			cl_k_PCPCF.SetArgT(7,expectedDF);
			cl_k_PCPCF.Enqueue(globalWorkSize);

			OpenCLFFT->Enqueue(clPCPCFResult,clImage1,CLFFT_BACKWARD);

			clEnqueueReadBuffer( clState::clq->cmdQueue, clImage1, CL_TRUE, 0, sizeX*sizeY*sizeof(std::complex<float>) , &data[ 0 ], 
					0, NULL, NULL );

			float maxheight = 0.0f;
			int maxPosition1 = 0;

			for(int j = 0; j< sizeX*sizeY;j++)
			{	
				// Use real part of ifft to get peak heights, not absolute....
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
	
			float xShift;
			float yShift;
			float subXShift;
			float subYShift;

			// Shift is positive or negative depending on image quadrant it appears in.
			if(maxindexr > floor(sizeY/2 +.5))
				yShift = maxindexr - (sizeY) -1;
			else
				yShift = maxindexr - 1;

			if(maxindexc > floor(sizeX/2 +.5))
				xShift = maxindexc - (sizeX) -1;
			else
				xShift = maxindexc - 1;


			PCPCFLib::FindVertexParabola(maxPosition1,sizeX,sizeY, data, subXShift, subYShift, maxheight);

			if(maxheight > bestheight) {
					bestheight = maxheight;
					bestscale = scale;
			}

			Debug("Height = "+Lex(maxheight)+" for "+Lex(scale));
			results[trial] = maxheight;
		}

		Debug("BEST SCALE IS "+Lex(bestscale));
		ResultLocker.~ImageDataLocker();

		ResultGraph.SetDimensionScale(0,(maxscale-minscale)/numberoftrials);

		// Origin is the point at which the zero poisition would be.

		ResultGraph.SetDimensionOrigin(0,-1.0f/((maxscale-minscale)/numberoftrials));
		ResultGraph.GetOrCreateImageDocument().Show();
	}
	
	PLUG_IN_EXIT
}



