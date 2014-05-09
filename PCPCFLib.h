#pragma once
#include <complex>
#include "Utility.h"
//#include "PCPCFFunction.h"
#include "clState.h"
#include "clFourier.h"


class PCPCFLib
{
public:
	PCPCFLib(void);
	~PCPCFLib(void);

	size_t* globalWorkSize;
	clKernel* PCPCF;
	clFourier* FFT;
public:
	
	static void FindVertexParabola(int maxposition1, int sizeX, int sizeY, std::vector<std::complex<float>> &data,
				float &xoffset, float &yoffset, float &peakheight);

	static void FindVertexParabolaMI(int maxposition1, int sizeX, int sizeY, float *data,
				float &xoffset, float &yoffset, float &peakheight);

	static void GetPeak(float* peak,int maxindexc,int maxindexr,int sizeX, int sizeY, std::vector<std::complex<float>> &data,
				float &xpos, float &ypos);

	static void GetShifts(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight,
				std::vector<std::complex<float>> &data, int sizeX, int sizeY);

	static void GetShifts(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight,
				std::vector<std::complex<float>> &data, int sizeX, int sizeY, float maxshift);

	static void GetShiftsMI(int &xShift, int &yShift, float &subXShift, float &subYShift,float &maxheight,
				float* data, int sizeX, int sizeY, float maxshift);

	static void PrepareImages(float* Data, int width, int height, int numberofimages);

	//void PhaseCompensatedPCF(int numberoftrials, float expectedDF, PCPCFOptions options, cl_mem &clPCPCFResult, cl_mem &clImage1, int sizeX, int sizeY, int imagenumber, std::vector<std::complex<float>> &dataOne, int* xShiftVals, int* yShiftVals, float* subXShifts, float* subYShifts, float* defocusshifts, int preshiftx, int preshifty);
};

