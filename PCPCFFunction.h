#include "clKernel.h"
#include "clFourier.h"

struct PCPCFOptions
{
	bool determinefocus;
	bool pcfrecon;
	bool rotscale;
	float focalstep;
	int steps;
	int reference;
	float pcpcfkmax;
	float searchpercentage;
	float magcal;
	float rotcal;
	float maxdrift;
};

struct Abberrations
{
	float beta;
	float delta;
	float Cs;
	float kmax;
	float A1r;
	float A1i;
};

void PCPCFWrapper2(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice *cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb);

void PCPCFWrapper(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice *cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb);

// Sometimes getting funny results after a certain point in code, not sure which step is failing but it starts giving same defocus for every image which should be impossible...
void PCPCFWrapper3(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb);

// Includes moving registration window
void PCPCFWrapper4(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb);

// Includes moving registration window
//void PCPCFWrapperReWrite(float* seriesdata, int numberOfImages, int xDim, int yDim, int iTop, int iLeft, int iBottom, int iRight, int* xShiftVals, int* yShiftVals, 
//				  float* subXShifts, float* subYShifts, float* defocusshifts, cl_context &context, clDevice* cldev, clFourier* FFT, cl_mem &clxFrequencies, cl_mem &clyFrequencies, clQueue* clq, PCPCFOptions &options, float wavelength, Abberrations &Abb);