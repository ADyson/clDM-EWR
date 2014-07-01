#include "clFourier.h"
#include "clKernel.h"
#include "ROIPositions.h"

namespace PCFRegistration
{
	struct PCFOptions
	{
		bool determinefocus;
		bool pcfrecon;
		bool rotscale;
		bool MI;
		float focalstep;
		int steps;
		int reference;
		float pcpcfkmax;
		float searchpercentage;
		float magcal;
		float rotcal;
		float maxdrift;
		float startdf;
		float enddf;
		float snr;
		bool knownfocus;
		std::vector<float> focuslist;
	};

	struct Abberrations2
	{
		float beta;
		float delta;
		float Cs;
		float kmax;
		float A1r;
		float A1i;
	};
};

class clRegistrationMemories
{
public:
	// First Group
	cl_mem clImage1;
	cl_mem clImage2;
	cl_mem clFFTImage1;
	cl_mem clFFTImage2;
	cl_mem clPCPCFResult;
	cl_mem fullImage;
	cl_mem rotScaleImage;
	cl_mem clW;
	cl_mem clWminus;
	cl_mem clw;
	cl_mem clwminus;
	cl_mem clT;
	cl_mem clU;
	cl_mem clV;
	cl_mem clRestored;
	cl_mem clMTF;
	cl_mem clNPS;

	void SetupGroupOne(int width, int height, int fullwidth, int fullheight);
	void SetupMTFNPS(int mtflength, int npslength);

	// Second Group

	cl_mem clRestoredMinus;
	cl_mem clQ;
	cl_mem clPCI;
	cl_mem clPCIC;
	cl_mem clPCIM;
	cl_mem clPCIS;
	cl_mem clSumOutput;
	cl_mem clSumOutputFloat;
	cl_mem clSumOutputUint;

	void SetupGroupTwo(int width, int height, int nGroups);
	void CleanUp(bool mtfnps);
};

#pragma once
class Registration
{
private:
	int gonext;
	int nextdir;

	int padTop;
	int padLeft;
	int padRight;
	int padBottom;
	int currentimage;
	bool gotMTF;
	bool gotNPS;
	int mtflength;
	float scaleMTFx;
	float scaleMTFy;
	float voltage;

	

public:
	Registration(void);
	~Registration(void);

	cl_mem ReductionResult;
	cl_mem ReductionPosition;


	bool noalign;
	// Vector to hold IDs of registered images as they are completed.
	std::vector<int> ImageList;

	// Open CL Memories
	clRegistrationMemories clMem;

	// Other options/settings
	PCFRegistration::PCFOptions options;
	PCFRegistration::Abberrations2 Aberrations;

	// Open CL Kernels
	clFourier* OpenCLFFT;
	clKernel cl_k_PCPCF;
	clKernel cl_k_PadCrop;
	clKernel cl_k_WaveTransferFunction;
	clKernel cl_k_WaveTransferFunctionMinus;
	clKernel cl_k_WienerW;
	clKernel cl_k_WienerWMinus;
	clKernel cl_k_WienerT;
	clKernel cl_k_WienerU;
	clKernel cl_k_WienerV;
	clKernel cl_k_MakeRestored;
	clKernel cl_k_RotScale;
	clKernel cl_k_Q;
	clKernel cl_k_Q2;
	clKernel cl_k_SumReduction;
	clKernel cl_k_SumReductionFloat;
	clKernel cl_k_SumReductionUint;
	clKernel cl_k_MinusWavefunction;
	clKernel cl_k_CalculatePCI;
	clKernel cl_k_Abs;
	clKernel cl_k_WaveTransferFunctionMTF;
	clKernel cl_k_WaveTransferFunctionMinusMTF;
	clKernel cl_k_MakeRestoredMTFNPS;
	clKernel cl_k_HanningWindow;
	clKernel cl_k_Predicted;
	clKernel cl_k_MaxReduction;

	void KernelCleanUp();

	// Things to set before using class
	void Setup(clFourier* FFT,int ref, int num, int wid, int hei, int fullwid, int fullhei, int* xs, int* ys, float* sxs, float* sys, float* dfs, float wavel, float pixelscale, float volts, float min, float max);
	int referenceimage;
	
	int NumberOfImages;
	int width;
	int height;
	int fullwidth;
	int fullheight;
	int* xShiftVals;
	int* yShiftVals;
	float* subXShifts;
	float* subYShifts;
	float* defocusshifts;
	float wavelength;
	float pixelscale;
	float min;
	float max;

	void IterateImageNumber();

	void BuildKernels();
	void LoadMTFNPS(DigitalMicrograph::Image &MTFImage, DigitalMicrograph::Image &NPSImage);
	void SetFixedArgs(ROIPositions ROIpos, cl_mem &clxFrequencies, cl_mem &clyFrequencies, float wavelength);
	void RotationScale(float* seriesdata,size_t* fullWorkSize,std::vector<float> &rotscaledseries);
	void RegisterSeries(float* seriesdata, ROIPositions ROIpos, cl_mem &clxFrequencies, cl_mem &clyFrequencies);
	void MIRegisterSeries(float* seriesdata, ROIPositions ROIpos, cl_mem &clxFrequencies, cl_mem &clyFrequencies);
	void MagnificationPCF(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, int preshiftx, int preshifty, size_t* globalWorkSize, std::vector<float> &seriesdata, int iLeft, int iTop, int im );
	void PhaseCompensatedPCF(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, int preshiftx, int preshifty, size_t* globalWorkSize);
	void MutualInformation(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, std::vector<std::complex<float>> &dataTwo, int preshiftx, int preshifty, size_t* globalWorkSize, int miSize, DigitalMicrograph::Image &MIMap, int imagenumber);
	void MutualInformationFast(int numberoftrials, float expectedDF, std::vector<std::complex<float>> &dataOne, std::vector<std::complex<float>> &dataTwo, int preshiftx, int preshifty, size_t* globalWorkSize, int miSize, DigitalMicrograph::Image &MIMap, int imagenumber);
	void AddToReconstruction(std::vector<float> &rotscaledseries, size_t* globalWorkSize, float DfGuess, float A1rGuess, float A1iGuess);
	void MakeDriftCorrectedSeries(std::vector<float> &rotscaledseries, size_t* globalWorkSize);
	void SetupPCI(size_t* globalWorkSize);
	void Window(cl_mem &Image, int width, int height);
	void PCITrials(float sumQ, float startdf, float stepdf,float A1r, float A1i, cl_mem &clxFrequencies, cl_mem &clyFrequencies, size_t* globalWorkSize, size_t* globalSizeSum,
			size_t* localSizeSum, int nGroups, int totalSize, float &DfGuess, float &A1rGuess, float &A1iGuess,bool first);

	// Little utility functions or wrappers
	float CalculateExpectedDF(int imageone, int imagetwo);
	float CalculateTrialDF(int trialnumber, int numberoftrials, float expectedDF);
	void CopyImageData(int imagenumber, int iLeft, int iTop,std::vector<std::complex<float>> &data, std::vector<float> &seriesdata, int preshiftx, int preshifty);
	void DeterminePadding(ROIPositions ROIpos);
	float SumReduction(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize);
	float SumReductionFloat(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize, int offset);
	float SumReductionUint(cl_mem &Array, size_t* globalSizeSum, size_t* localSizeSum, int nGroups, int totalSize, int offset);
	
};


