// AdjustRecon.cpp : implementation file
//

#include "stdafx.h"
#include "AdjustRecon.h"
#include "FTSR.h"
#include "boost\lexical_cast.hpp"
#include "Utility.h"
#include "ClKernelsAdjust.h"
#include "clState.h"


// AdjustRecon dialog

IMPLEMENT_DYNAMIC(AdjustRecon, CDialog)

AdjustRecon::AdjustRecon(CWnd* pParent /*=NULL*/)
	: CDialog(AdjustRecon::IDD, pParent)
{
	mParent = pParent;
	
	C1 = "0";
	C3 = "0";
	A1r = "0";
	A1i = "0";
	A2r = "0";
	B2r = "0";
	A3r = "0";
	S3r = "0";
	A2i = "0";
	B2i = "0";
	A3i = "0";
	S3i = "0";

}

AdjustRecon::~AdjustRecon()
{
}

void AdjustRecon::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_C1, m_C1);
	DDX_Text(pDX,IDC_C1,C1);
	DDX_Control(pDX, IDC_C3, m_C3);
	DDX_Text(pDX,IDC_C3,C3);
	DDX_Control(pDX, IDC_A1R, m_A1r);
	DDX_Text(pDX,IDC_A1R,A1r);
	DDX_Control(pDX, IDC_A1I, m_A1i);
	DDX_Text(pDX,IDC_A1I,A1i);
	DDX_Control(pDX, IDC_A2R, m_A2r);
	DDX_Text(pDX,IDC_A2R,A2r);
	DDX_Control(pDX, IDC_B2R, m_B2r);
	DDX_Text(pDX,IDC_B2R,B2r);
	DDX_Control(pDX, IDC_A3R, m_A3r);
	DDX_Text(pDX,IDC_A3R,A3r);
	DDX_Control(pDX, IDC_S3R, m_S3r);
	DDX_Text(pDX,IDC_S3R,S3r);
	DDX_Control(pDX, IDC_A2I, m_A2i);
	DDX_Text(pDX,IDC_A2I,A2i);
	DDX_Control(pDX, IDC_B2I, m_B2i);
	DDX_Text(pDX,IDC_B2I,B2i);
	DDX_Control(pDX, IDC_A3I, m_A3i);
	DDX_Text(pDX,IDC_A3I,A3i);
	DDX_Control(pDX, IDC_S3I, m_S3i);
	DDX_Text(pDX,IDC_S3I,S3i);
}


BEGIN_MESSAGE_MAP(AdjustRecon, CDialog)
	ON_BN_CLICKED(IDC_MAKECORRECTED, &AdjustRecon::OnBnClickedAdjust)
END_MESSAGE_MAP()


// AdjustRecon message handlers

// This is apparently needed for Modeless dialogs to not have a memory leak when they are closed.
void AdjustRecon::PostNcDestroy()
{
	CDialog::PostNcDestroy();
	delete this;
}

void AdjustRecon::OnBnClickedAdjust()
{
	// Get Reconstructed Image - > Make a new one with aberrations corrected
	UpdateData(true);

	EWR* ParentClass = (EWR*)mParent;

	clQueue* clq = clState::clq;
	cl_context context = clState::context;
	clDevice* cldev = clState::cldev;

	DigitalMicrograph::Image Recon = DigitalMicrograph::GetFrontImage();
	DigitalMicrograph::TagGroup ReconTags = Recon.GetTagGroup();

	float wavelength;
	float kmax;

	try
	{
		ReconTags.GetTagAsFloat("Reconstruction:Wavelength",&wavelength);
		ReconTags.GetTagAsFloat("Reconstruction:Kmax",&kmax);
	}
	catch(...)
	{
		DigitalMicrograph::Result("Are you sure this is a reconstructed image, some tags are missing \n");
	}

	// Get image size and calibrations.

	Gatan::uint32 xDim; // Dont use these in a sum that can go negative, will break things i.e.  (1-xdim)/2 wont work...
	Gatan::uint32 yDim; // Dont use these in a sum that can go negative, will break things i.e.  (1-xdim)/2 wont work...

	Recon.GetDimensionSizes(xDim,yDim);
	float reconpixelscale = Recon.GetDimensionScale(0);

	int xSize = xDim;
	int ySize = yDim;


	// Get abberation values from dialog.

	float fC1 = boost::lexical_cast<float>(C1);
	float fC3 = boost::lexical_cast<float>(C3);
	float fA1r = boost::lexical_cast<float>(A1r);
	float fA1i = boost::lexical_cast<float>(A1i);
	float fA2r = boost::lexical_cast<float>(A2r);
	float fB2r = boost::lexical_cast<float>(B2r);
	float fA3r = boost::lexical_cast<float>(A3r);
	float fS3r = boost::lexical_cast<float>(S3r);
	float fA2i = boost::lexical_cast<float>(A2i);
	float fB2i = boost::lexical_cast<float>(B2i);
	float fA3i = boost::lexical_cast<float>(A3i);
	float fS3i = boost::lexical_cast<float>(S3i);

	// Make frequency arrays

	float* xFrequencies = new float[xDim];
	float* yFrequencies = new float[yDim];

	int midX = ceil(float(xDim)/2);
	int midY = ceil(float(yDim)/2);
	

	for(int i = 1 ; i <= xDim ; i++)
	{
		if(i <= midX)
			xFrequencies[i-1] = (i - 1)/(reconpixelscale * xSize);
		else xFrequencies[i-1] = (i - 1 - xSize)/(reconpixelscale * xSize);
	}



	for(int i = 1 ; i <= yDim ; i++)
	{
		if(i <= midY)
			yFrequencies[i-1] = (i - 1)/(reconpixelscale * ySize);
		else 
		{
			float a = (i - 1 - ySize);
			yFrequencies[i-1] = a/(reconpixelscale * ySize);
		}
	}

	// Upload to GPU

	cl_mem clxFrequencies = clCreateBuffer ( context, CL_MEM_READ_WRITE, xDim * sizeof(cl_float), 0, &status);
	cl_mem clyFrequencies = clCreateBuffer ( context, CL_MEM_READ_WRITE, yDim * sizeof(cl_float), 0, &status);

	clEnqueueWriteBuffer( clq->cmdQueue, clxFrequencies, CL_FALSE, 0, xDim*sizeof(float), &xFrequencies[0], 
				0, NULL, NULL );
	clEnqueueWriteBuffer( clq->cmdQueue, clyFrequencies, CL_TRUE, 0, yDim*sizeof(float), &yFrequencies[0], 
				0, NULL, NULL );

	cl_mem clReconstruction = clCreateBuffer ( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem clReconstruction2 = clCreateBuffer ( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	cl_mem clCorrected = clCreateBuffer ( context, CL_MEM_READ_WRITE, xDim * yDim * sizeof(cl_float2), 0, &status);
	;



	// Upload the image to the GPU

	// FFT it

	


	Gatan::PlugIn::ImageDataLocker ReconLocker(Recon);
	Gatan::complex64* data = (Gatan::complex64*)ReconLocker.get();

	std::vector<cl_float2> datavec(xSize*ySize);

	
	for(int i = 1 ; i < xDim * yDim ; i++)
	{
		datavec[i].s[0] = data[i].real() - 1;
		datavec[i].s[1] = data[i].imag();
	}
	
	clFourier* FFT = new clFourier(context, clq);
	FFT->Setup(xSize,ySize);

	
	clEnqueueWriteBuffer(clq->cmdQueue, clReconstruction , CL_TRUE, 0 , xDim * yDim * sizeof(cl_float2) , &datavec[0] , 0 , NULL , NULL);

	ReconLocker.~ImageDataLocker();

	FFT->Enqueue(clReconstruction,clReconstruction2,CLFFT_FORWARD);

	

	// Create kernel to do reconstruction correction
	clKernel* clThirdOrderCorrection = new clKernel(ThirdOrderCorrectionsource,context,cldev,"clThirdOrderCorrection",clq);
	clThirdOrderCorrection->BuildKernel();

	clThirdOrderCorrection->SetArgT(0,clReconstruction2);
	clThirdOrderCorrection->SetArgT(1,clCorrected);
	clThirdOrderCorrection->SetArgT(2,clxFrequencies);
	clThirdOrderCorrection->SetArgT(3,clyFrequencies);
	clThirdOrderCorrection->SetArgT(4,xSize);
	clThirdOrderCorrection->SetArgT(5,ySize);
	clThirdOrderCorrection->SetArgT(6,wavelength);
	clThirdOrderCorrection->SetArgT(7,fA1r);
	clThirdOrderCorrection->SetArgT(8,fA1i);
	clThirdOrderCorrection->SetArgT(9,fA2r);
	clThirdOrderCorrection->SetArgT(10,fA2i);
	clThirdOrderCorrection->SetArgT(11,fA3r);
	clThirdOrderCorrection->SetArgT(12,fA3i);
	clThirdOrderCorrection->SetArgT(13,fB2r);
	clThirdOrderCorrection->SetArgT(14,fB2i);
	clThirdOrderCorrection->SetArgT(15,fS3r);
	clThirdOrderCorrection->SetArgT(16,fS3i);
	clThirdOrderCorrection->SetArgT(17,fC1);
	clThirdOrderCorrection->SetArgT(18,fC3);
	clThirdOrderCorrection->SetArgT(19,kmax);

	size_t* globalWorkSize = new size_t[3];

	globalWorkSize[0] = xSize;
	globalWorkSize[1] = ySize;
	globalWorkSize[2] = 1;

	clThirdOrderCorrection->Enqueue(globalWorkSize);


	// FFT back
	FFT->Enqueue(clCorrected,clReconstruction,CLFFT_BACKWARD);

	DigitalMicrograph::Image CorrectedImage = Utility::PrintCLMemToImagePlusOne(clReconstruction,"Adjusted Reconstruction",xSize,ySize,clFloat2,clq);

	CorrectedImage.GetImageDisplay(0).SetComplexMode(5);

	clThirdOrderCorrection->~clKernel();
	clReleaseMemObject(clCorrected);
	clReleaseMemObject(clReconstruction);
	clReleaseMemObject(clReconstruction2);
}
