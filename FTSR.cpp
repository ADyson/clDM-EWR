// EWR.cpp : implementation file
//
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include "stdafx.h"
#include "FTSR.h"
#include "boost/lexical_cast.hpp"
#include <sstream>
#include "DMPluginCamera.h"
#include "standardfunctions.h"
#include "CLKernels.h"
#include "clFourier.h"
#include "OptionsDialog.h"
#include "utility.h"
#include "clState.h"
#include "AdjustRecon.h"
#include "PCPCFLib.h"




// GPU usually platform 0 dev 0, cpu is 0,1 with AMD Gfx or 1,0 with nVidia,
// cant see Gfx over Remote Desktop though so CPU 0,0

// EWR dialog

IMPLEMENT_DYNAMIC(EWR, CDialog)

EWR::EWR(CWnd * pParent /*=NULL*/)
	: CDialog(EWR::IDD, pParent)
	, m_Beta_Text(_T("0.5"))
	, m_Cs_Text(_T("1"))
	, m_Delta_Text(_T("3"))
	, m_Inf_Text(_T("5"))
	, m_SearchEnd_Text(_T("50"))
	, m_SearchStart_Text(_T("-50"))
	, m_SearchSteps_Text(_T("50")) {

	// Set default starting parameters
	PCFWrapper.options.determinefocus = true;
	PCFWrapper.options.MI = false;
	PCFWrapper.options.steps = 5;
	PCFWrapper.options.reference = 1; // Starts at zero
	PCFWrapper.options.pcpcfkmax = 4;
	PCFWrapper.options.searchpercentage = 30;
	PCFWrapper.options.pcfrecon = true;
	PCFWrapper.options.rotscale = false;
	PCFWrapper.options.magcal = 1.00;
	PCFWrapper.options.rotcal = 0;
	PCFWrapper.options.startdf = -50;
	PCFWrapper.options.enddf = 50;
	PCFWrapper.options.snr = 0.2;
	PCFWrapper.options.maxdrift = 0;

	// Setting default values to pass into options.
	maxdrift = "";
	GotMTF = false;
	GotNPS = false;

}

EWR::~EWR() {
}

void EWR::DoDataExchange(CDataExchange * pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT1, m_Inf);
	DDX_Text(pDX, IDC_EDIT1, m_Inf_Text);
	DDX_Control(pDX, IDC_EDIT2, m_SearchStart);
	DDX_Text(pDX, IDC_EDIT2, m_SearchStart_Text);
	DDX_Control(pDX, IDC_EDIT3, m_SearchEnd);
	DDX_Text(pDX, IDC_EDIT3, m_SearchEnd_Text);
	// This one is an integer only, no control required
	DDX_Text(pDX, IDC_EDIT4, m_SearchSteps_Text);
	DDX_Control(pDX, IDC_EDIT5, m_Cs);
	DDX_Text(pDX, IDC_EDIT5, m_Cs_Text);
	DDX_Control(pDX, IDC_EDIT6, m_Beta);
	DDX_Text(pDX, IDC_EDIT6, m_Beta_Text);
	DDX_Control(pDX, IDC_EDIT7, m_Delta);
	DDX_Text(pDX, IDC_EDIT7, m_Delta_Text);
}


BEGIN_MESSAGE_MAP(EWR, CDialog)
	ON_BN_CLICKED(IDC_BUTTON1, &EWR::ResetOpenCLDevice)
	ON_BN_CLICKED(IDC_BUTTON2, &EWR::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_OPTIONSBUTTON, &EWR::OnBnClickedOptions)
	//ON_BN_CLICKED(IDC_BUTTON3, &EWR::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_OPTIONSBUTTON2, &EWR::OnBnClickedOptionsbutton2)
	ON_BN_CLICKED(IDC_Adjust, &EWR::OnBnClickedAdjust)
	ON_BN_CLICKED(IDC_BUTTONNPS, &EWR::OnBnClickedButtonnps)
	ON_BN_CLICKED(IDC_BUTTONMTF, &EWR::OnBnClickedButtonmtf)
END_MESSAGE_MAP()

BOOL EWR::PreTranslateMessage(MSG * pMsg) {
	if (pMsg->message == WM_KEYDOWN) {
		if ((pMsg->wParam == VK_RETURN) || (pMsg->wParam == VK_ESCAPE))
			pMsg->wParam = VK_TAB;
	}

	return CDialog::PreTranslateMessage(pMsg);
}


// From Worker Class
void EWR::DoWork() {
	if (!clState::OpenCLAvailable) {
		Utility::SetResultWindow("OpenCL is not available on this system\n");
	}

	// At this point we want the front image to be a stack of any size;
	DigitalMicrograph::Image frontImage;

	try {
		frontImage = DigitalMicrograph::GetFrontImage();
	} catch (...) {
		DigitalMicrograph::Result("No front image \n");
		return;
	}

	// Get dimensions so we know how big to make the slices;
	unsigned long xDim;
	unsigned long yDim;
	unsigned long numberOfImages;

	frontImage.GetDimensionSizes(xDim, yDim, numberOfImages);

	// Check for ROI and if not create one.

	int numROI;
	DigitalMicrograph::ImageDisplay imDisp;
	imDisp = frontImage.GetImageDisplay(0);

	numROI = imDisp.CountROIs();

	DigitalMicrograph::ROI regionOfInterest;

	if (numROI == 0) {
		// Create image display
		regionOfInterest = DigitalMicrograph::NewROI();
		regionOfInterest.SetRectangle(0, 0, yDim, xDim);
		imDisp.AddROI(regionOfInterest);
	}

	else {
		regionOfInterest = imDisp.GetROI(0);
	}

	// Get Region Positions
	float top;
	float left;
	float bottom;
	float right;

	// Check its a rectangle ROI
	if (!regionOfInterest.IsRectangle()) {
		DigitalMicrograph::OpenAndSetProgressWindow("ROI is not a rectangle", "", "");
		return;
	}


	regionOfInterest.GetRectangle(&top, &left, &bottom, &right);

	// Do Some Checks
	int iTop = boost::numeric_cast<int>(top);
	int iLeft = boost::numeric_cast<int>(left);
	int iBottom = boost::numeric_cast<int>(bottom);
	int iRight = boost::numeric_cast<int>(right);

	int sizeX = iRight - iLeft;
	int sizeY = iBottom - iTop;

	ROIPositions ROIposition (iTop, iLeft, iBottom, iRight);

	// Check its a power of two size aswell

	//TODO: Add padding option to do other sizes
	if (!(sizeX == 64 || sizeX == 128 || sizeX == 256 || sizeX == 512 || sizeX == 1024 || sizeX == 2048 || sizeX == 4096)) {
		DigitalMicrograph::Result("ROI must be power of 2 size");
		return;
	}

	if (!(sizeY == 64 || sizeY == 128 || sizeY == 256 || sizeY == 512 || sizeY == 1024 || sizeY == 2048 || sizeY == 4096)) {
		DigitalMicrograph::Result("ROI must be power of 2 size");
		return;
	}

	// Need pixelscale, and voltage, and focalstep to get right parameters.
	DigitalMicrograph::TagGroup imagetags = frontImage.GetTagGroup();

	float fStep;
	imagetags.GetTagAsFloat("Focal Series:Adjusted focalstep", &fStep);

	float voltage;
	imagetags.GetTagAsFloat("Microscope Info:Voltage", &voltage);

	float pixelscale = frontImage.GetDimensionScale(0);

	PCFWrapper.options.focalstep = fStep;

	// Setup two frequency arrays for this size.
	float * xFrequencies = new float[sizeX];
	float * yFrequencies = new float[sizeY];

	int midX = ceil(float(sizeX) / 2);
	int midY = ceil(float(sizeY) / 2);

	for (int i = 1 ; i <= sizeX ; i++) {
		if (i <= midX)
			xFrequencies[i - 1] = (i - 1) / (pixelscale * sizeX);
		else xFrequencies[i - 1] = (i - 1 - sizeX) / (pixelscale * sizeX);
	}

	for (int i = 1 ; i <= sizeY ; i++) {
		if (i <= midY)
			yFrequencies[i - 1] = (i - 1) / (pixelscale * sizeY);
		else yFrequencies[i - 1] = (i - 1 - sizeY) / (pixelscale * sizeY);
	}

	// Calculate Electron Wavelength in nm

	float e = 1.6e-19;
	float wavelength = 6.63e-34 * 3e+8 * 1e+9 / sqrt((e * voltage * (2 * 9.11e-31 * 9e+16 + e * voltage)));

	// Setup arrays to store shift values.
	int * xShiftVals = new int[numberOfImages];
	int * yShiftVals = new int[numberOfImages];
	float * subXShifts = new float[numberOfImages];
	float * subYShifts = new float[numberOfImages];
	float * defocusshifts = new float[numberOfImages];

	xShiftVals[0] = 0;
	yShiftVals[0] = 0;
	subXShifts[0] = 0;
	subYShifts[0] = 0;
	defocusshifts[0] = 0; // Doesnt really have one.

	// Setup a fourier transform

	FFT = new clFourier(clState::context, clState::clq);
	FFT->Setup(sizeX, sizeY);

	// Upload frequencies to GPU!
	cl_mem clxFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeX * sizeof(cl_float), 0, &clState::status);
	cl_mem clyFrequencies = clCreateBuffer ( clState::context, CL_MEM_READ_WRITE, sizeY * sizeof(cl_float), 0, &clState::status);

	clEnqueueWriteBuffer( clState::clq->cmdQueue, clxFrequencies, CL_FALSE, 0, sizeX * sizeof(float), &xFrequencies[0],
	                      0, NULL, NULL );
	clEnqueueWriteBuffer( clState::clq->cmdQueue, clyFrequencies, CL_TRUE, 0, sizeY * sizeof(float), &yFrequencies[0],
	                      0, NULL, NULL );


	// Get pointer to image data
	Gatan::PlugIn::ImageDataLocker seriesLocker(frontImage);
	float * seriesdata = (float *) seriesLocker.get();

	// TODO: Get these options from dialog etc... actually use them in code (kmax)
	PCFWrapper.Aberrations.A1i = 0;
	PCFWrapper.Aberrations.A1r = 0;
	PCFWrapper.Aberrations.beta = boost::lexical_cast<float>(m_Beta_Text) / 1000;
	PCFWrapper.Aberrations.delta = boost::lexical_cast<float>(m_Delta_Text);
	PCFWrapper.Aberrations.Cs = boost::lexical_cast<float>(m_Cs_Text);
	PCFWrapper.Aberrations.kmax = boost::lexical_cast<float>(m_Inf_Text);
	PCFWrapper.options.startdf = boost::lexical_cast<float>(m_SearchStart_Text);
	PCFWrapper.options.enddf = boost::lexical_cast<float>(m_SearchEnd_Text);

	// Reference Image has 0 shifts and focus etc...
	// If not doing reconstruction PCPCF then the reference should always be 1 (Need to make sure defocus difference is correctly reported aswell - distance from reference or distance from neighbour
	xShiftVals[PCFWrapper.options.reference] = 0;
	yShiftVals[PCFWrapper.options.reference] = 0;
	subXShifts[PCFWrapper.options.reference] = 0;
	subYShifts[PCFWrapper.options.reference] = 0;
	defocusshifts[PCFWrapper.options.reference] = 0; // Doesnt really have one.


	if (maxdrift == "") {
		PCFWrapper.options.maxdrift = 0;
	} else {
		PCFWrapper.options.maxdrift = boost::lexical_cast<float>(maxdrift);
	}

	std::string shouldalign;
	imagetags.GetTagAsString("Focal Series:No Align", shouldalign);

	if (shouldalign == "True" || shouldalign == "true")
		PCFWrapper.noalign = true;
	else
		PCFWrapper.noalign = false;

	std::string isnormalized;
	imagetags.GetTagAsString("Focal Series:Normalized", isnormalized);

	std::string focuslist;
	imagetags.GetTagAsString("Focal Series:Known Focus", focuslist);

	if (focuslist == "true" || focuslist == "True") {
		PCFWrapper.options.knownfocus = true;

		for (int i = 1; i < numberOfImages; i++) {
			std::string focusvalue;
			std::string path = "Focal Series:Focus Values:" + Lex(i);
			imagetags.GetTagAsString(path.c_str(), focusvalue);
			PCFWrapper.options.focuslist.push_back( boost::lexical_cast<float>(focusvalue));
		}
	} else {
		PCFWrapper.options.knownfocus = false;
	}


	if (!(isnormalized == "true" || isnormalized == "True")) {
		// Normalise images and perform hot pixel removal.
		PCPCFLib::PrepareImages(seriesdata, xDim, yDim, numberOfImages);
		imagetags.SetTagAsString("Focal Series:Normalized", "true");
	}

	frontImage.DataChanged();

	// Get contrast limits incase we are goign to do MI...
	float min, max;
	ImageDisplayGetContrastLimits(frontImage.GetImageDisplay(0), &min, &max);

	// Expand min and max slightly to make sure it can cover every possible value.
	//min*=1.05;
	//max*=1.05;

	// If we have an mtf and nps can pass these to PCFWrapper.
	// Must be before setup because setup sets clKernel arguments.
	if (GotMTF && GotNPS) {
		PCFWrapper.LoadMTFNPS(MTFImage, NPSImage);
	}

	PCFWrapper.Setup(FFT, PCFWrapper.options.reference, numberOfImages, sizeX, sizeY, xDim, yDim, xShiftVals, yShiftVals, subXShifts,
	                 subYShifts, defocusshifts, wavelength, pixelscale, voltage, min, max);

	// Now to reconstruction
	if (PCFWrapper.options.MI) {
		PCFWrapper.MIRegisterSeries(seriesdata, ROIposition, clxFrequencies, clyFrequencies);
	} else {
		PCFWrapper.RegisterSeries(seriesdata, ROIposition, clxFrequencies, clyFrequencies);
	}

	// Now we have relative shifts of all the images.
	// Print them to screen

	for (int i = 0 ; i < numberOfImages ; i++) {
		Utility::SetResultWindow("X Shift for image " + Lex(i + 1) + " : " + Lex(subXShifts[i]) + "\n");
		Utility::SetResultWindow("Y Shift for image " + Lex(i + 1) + " : " + Lex(subYShifts[i]) + "\n");
		Utility::SetResultWindow("Defocus for image " + Lex(i + 1) + " : " + Lex(defocusshifts[i]) + "\n");
	}


	//TODO: convert this to line plot graph with 2 lines...
	DigitalMicrograph::Image Xgraph = DigitalMicrograph::RealImage("X Shift", 4, numberOfImages);
	Gatan::PlugIn::ImageDataLocker XLocker(Xgraph);
	float * XData = (float *) XLocker.get();

	for (int i = 0 ; i < numberOfImages ; i++) {
		XData[i] = subXShifts[i];
	}

	DigitalMicrograph::Image Ygraph = DigitalMicrograph::RealImage("Y Shift", 4, numberOfImages);
	Gatan::PlugIn::ImageDataLocker YLocker(Ygraph);
	float * YData = (float *) YLocker.get();

	for (int i = 0 ; i < numberOfImages ; i++) {
		YData[i] = subYShifts[i];
	}

	XLocker.~ImageDataLocker();
	Xgraph.GetOrCreateImageDocument().Show();
	YLocker.~ImageDataLocker();
	Ygraph.GetOrCreateImageDocument().Show();

	// Clear non needed data.
	seriesLocker.~ImageDataLocker();
	clReleaseMemObject(clxFrequencies);
	clReleaseMemObject(clyFrequencies);
	FFT->~clFourier();

	delete[] xShiftVals;
	delete[] yShiftVals;
	delete[] subXShifts;
	delete[] subYShifts;
	delete[] defocusshifts;
}



// Remember to check for OpenCL availability at the start of function calls
void EWR::OnBnClickedButton2() {
	// Ensure data from controls is up to date
	UpdateData(TRUE);

	// Don't run two at once - that would be bad
	if (GetStatus() != Working) {
		Start();
	}
}

//void EWR::OnBnClickedCancel()
//{
//	if(GetStatus() == Working)
//	{
//		// Have to add points at which thread checks if it has been asked to stop.
//		Stop();
//	}
//}

void EWR::OnBnClickedOptions() {
	int a = 0;
	int b = 0;
	int c = 0;
	int d = 0;

	if (PCFWrapper.options.determinefocus) {
		a = 1;
	}

	if (PCFWrapper.options.pcfrecon) {
		b = 1;
	}

	if (PCFWrapper.options.rotscale) {
		c = 1;
	}

	if (PCFWrapper.options.MI) {
		d = 1;
	}

	COptionsDialog dlg(Lex(PCFWrapper.options.pcpcfkmax), Lex(PCFWrapper.options.steps), Lex(PCFWrapper.options.reference),
	                   Lex(PCFWrapper.options.searchpercentage), a, b, c, Lex(PCFWrapper.options.magcal), Lex(PCFWrapper.options.rotcal), maxdrift,
	                   Lex(PCFWrapper.options.snr), d, this);
	dlg.DoModal();
}

void EWR::UpdateOptions(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage,
                        int determinedf, int reconpcf) {
	if (determinedf == 1) {
		PCFWrapper.options.determinefocus = true;
	} else
		PCFWrapper.options.determinefocus = false;

	if (reconpcf == 1) {
		PCFWrapper.options.pcfrecon = true;
	} else
		PCFWrapper.options.pcfrecon = false;


	PCFWrapper.options.steps = boost::lexical_cast<int>(pcftrials);
	PCFWrapper.options.reference = boost::lexical_cast<int>(reference);
	PCFWrapper.options.pcpcfkmax = boost::lexical_cast<float>(pcfkmax);
	PCFWrapper.options.searchpercentage = boost::lexical_cast<float>(searchpercentage);
}

void EWR::UpdateOptions2(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage,
                         int determinedf, int reconpcf, int magrotfix, std::string magscale, std::string rotscale, std::string maxdrift, std::string snr,
                         int mi) {
	if (determinedf == 1) {
		PCFWrapper.options.determinefocus = true;
	} else
		PCFWrapper.options.determinefocus = false;

	if (reconpcf == 1) {
		PCFWrapper.options.pcfrecon = true;
	} else
		PCFWrapper.options.pcfrecon = false;

	if (magrotfix == 1) {
		PCFWrapper.options.rotscale = true;
	} else
		PCFWrapper.options.rotscale = false;

	if (mi == 1) {
		PCFWrapper.options.MI = true;
	} else
		PCFWrapper.options.MI = false;


	PCFWrapper.options.steps = boost::lexical_cast<int>(pcftrials);
	PCFWrapper.options.reference = boost::lexical_cast<int>(reference);
	PCFWrapper.options.pcpcfkmax = boost::lexical_cast<float>(pcfkmax);
	PCFWrapper.options.searchpercentage = boost::lexical_cast<float>(searchpercentage);
	PCFWrapper.options.magcal = boost::lexical_cast<float>(magscale);
	PCFWrapper.options.rotcal = boost::lexical_cast<float>(rotscale);
	PCFWrapper.options.snr = boost::lexical_cast<float>(snr);
	this->maxdrift = maxdrift;
}

// Open Image Series as Stack Button
void EWR::OnBnClickedOptionsbutton2() {
	std::vector<std::string> filepaths;
	std::vector<DigitalMicrograph::Image> images;

	CFileDialog FileDialog(TRUE, "*.*", NULL, OFN_HIDEREADONLY | OFN_ALLOWMULTISELECT,
	                       "All files (*.*)|*.*|DM3 Files (*.dm3)|*.dm3|DM4 Files (*.dm4)|*.dm4||");

	CString filenames;


	FileDialog.m_ofn.lpstrFile = filenames.GetBuffer(2048);
	FileDialog.m_ofn.nMaxFile  = 2048;

	if (FileDialog.DoModal() == IDOK) {
		POSITION pos = FileDialog.GetStartPosition();

		if (pos) {
			CString PathName = "";

			do {
				PathName = FileDialog.GetNextPathName(pos);
				std::string pname  = PathName;
				filepaths.push_back(pname);
			} while (pos);
		}
	}

	// Now have a vector of all filepaths..

	for (int i = 0 ; i < filepaths.size() ; i++) {
		Utility::SetResultWindow(filepaths[i]);
		images.push_back( DigitalMicrograph::OpenImage(filepaths[i]));
	}

	// All already in correct order :)

	// Open one to collect information...

	float xcalibration;
	float xorigin;
	DigitalMicrograph::String xunits;
	float ycalibration;
	float yorigin;
	DigitalMicrograph::String yunits;

	DigitalMicrograph::TagGroup ImageTags = images[0].GetTagGroup();

	images[0].GetDimensionCalibration(0, &xorigin, &xcalibration, xunits, 0);
	images[0].GetDimensionCalibration(1, &yorigin, &ycalibration, yunits, 0);

	Gatan::uint32 xsize;
	Gatan::uint32 ysize;

	images[0].GetDimensionSizes(xsize, ysize);

	// Check for consistency of data.
	for (int i = 0 ; i < images.size() ; i++) {
		Gatan::uint32 checkxsize;
		Gatan::uint32 checkysize;

		images[i].GetDimensionSizes(checkxsize, checkysize);

		if (checkxsize != xsize || checkysize != ysize) {
			Utility::SetResultWindow("Not all images have the same dimensions");

			// Get rid of all original images now...
			for (int i = 0 ; i < images.size() ; i++) {
				// Still need to close all the images even if exiting early.
				DigitalMicrograph::CloseImage(images[i]);
			}

			return;
		}
	}

	// Then make an image of the correct size and loop through adding all data...

	DigitalMicrograph::Image StackImage = DigitalMicrograph::RealImage("Focal Series Stack", 4, xsize, ysize, images.size());
	DigitalMicrograph::TagGroup StackTags = StackImage.GetTagGroup();

	Gatan::PlugIn::ImageDataLocker StackLocker(StackImage);
	float * stackdata = (float *) StackLocker.get();

	// Note - images could be not floating point type images....
	for (int i = 0 ; i < images.size() ; i++) {
		Gatan::PlugIn::ImageDataLocker ImageLocker(images[i]);
		float * imagedata = (float *) ImageLocker.get();

		for (int xpixel = 0 ; xpixel < xsize ; xpixel++)
			for (int ypixel = 0 ; ypixel < ysize ; ypixel++) {
				stackdata[i * xsize * ysize + ypixel * xsize + xpixel] = imagedata[ypixel * xsize + xpixel];
			}

		ImageLocker.~ImageDataLocker();
	}

	StackLocker.~ImageDataLocker();

	StackImage.GetOrCreateImageDocument().Show();

	// Copy tags and calibrations...
	DigitalMicrograph::TagGroupReplaceTagsWithCopy(StackTags, ImageTags);
	StackImage.SetDimensionCalibration(0, xorigin, xcalibration, xunits, 0);
	StackImage.SetDimensionCalibration(1, yorigin, ycalibration, yunits, 0);

	// PROFIT

	// Get rid of all original images now...
	for (int i = 0 ; i < images.size() ; i++) {
		DigitalMicrograph::CloseImage(images[i]);
	}
}


void EWR::OnBnClickedAdjust() {
	// TODO: Add your control notification handler code here
	// Launch the adjust recon dialog? Anyway for it to be launched without slowing everything down
	AdjustRecon * dlg = new AdjustRecon(this);
	dlg->Create(AdjustRecon::IDD);
	dlg->ShowWindow(SW_SHOW);
}

void EWR::ResetOpenCLDevice() {
	clAmdFftTeardown();
	clReleaseContext(clState::context);
	clReleaseCommandQueue(clState::clq->cmdQueue);

	DigitalMicrograph::TagGroup PersistentTags = DigitalMicrograph::GetPersistentTagGroup();
	DigitalMicrograph::String clPlatformTag;
	DigitalMicrograph::String clDeviceTag;


	DigitalMicrograph::TagGroupGetTagAsString(PersistentTags, "OpenCL:Platform", clPlatformTag);
	DigitalMicrograph::TagGroupGetTagAsString(PersistentTags, "OpenCL:Device", clDeviceTag);

	int platformnumber;
	int devicenumber;

	try {
		platformnumber = boost::lexical_cast<int>(clPlatformTag);
		devicenumber = boost::lexical_cast<int>(clDeviceTag);
	} catch (boost::bad_lexical_cast e) {
		DigitalMicrograph::Result("Open CL Device and Platform incorrectly set in global tags \n");
		DigitalMicrograph::Result("Set Platform and Device 0,1,2 etc in tags OpenCl:Platform and OpenCL:Device \n");
		DigitalMicrograph::Result("Run clinfo at command prompt to find available devices and platforms");
		return;
	}

	cl_uint numDevices;
	cl_device_id * devices;

	//Setup OpenCL
	clState::OpenCLAvailable = false;
	clState::context = NULL;
	numDevices = 0;
	devices = NULL;

	// Maybe Can Do OpenCL setup and device registering here - Print to Ouput with device data?
	// Discover and initialize available platforms
	cl_uint numPlatforms = 0;
	cl_platform_id * platforms = NULL;

	// Use clGetPlatformIds() to retrieve the number of platforms
	clState::status = clGetPlatformIDs(0, NULL, &numPlatforms);

	// Allocate enough space for each platform
	platforms = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));

	// Fill in platforms with clGetPlatformIDs()
	clState::status = clGetPlatformIDs(numPlatforms, platforms, NULL);

	// Discover and initialize available devices
	// use clGetDeviceIDs() to retrieve number of devices present
	clState::status = clGetDeviceIDs(platforms[platformnumber], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);

	// Allocate enough space for each device
	devices = (cl_device_id *)malloc(numDevices * sizeof(cl_device_id));

	// Fill in devices with clGetDeviceIDs()
	clState::status = clGetDeviceIDs(platforms[platformnumber], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);

	// Most of initialisation is done, would be nice to print device information...
	//Getting the device name
	size_t deviceNameLength = 4096;
	size_t actualSize;
	char * tempDeviceName = (char *)malloc(4096);
	char * deviceName;
	clState::status |= clGetDeviceInfo(devices[devicenumber], CL_DEVICE_NAME, deviceNameLength, tempDeviceName, &actualSize);

	if (clState::status == CL_SUCCESS) {
		deviceName = (char *)malloc(actualSize);
		memcpy(deviceName, tempDeviceName, actualSize);
		free(tempDeviceName);

		std::string devName(deviceName);
		DigitalMicrograph::Result("Using OpenCL on device " + devName + " - OCL - EWR\n");
		DigitalMicrograph::Result("To change edit the Global Tags OpenCL:Platform and OpenCL:Device then restart DM\n");
		clState::OpenCLAvailable = true;

		clState::context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &clState::status);
		clState::clq = new clQueue();
		clState::clq->SetupQueue(clState::context, devices[devicenumber]);
		clState::cldev = new clDevice(numDevices, devices[devicenumber]);
	}

	if (clState::status != CL_SUCCESS) {
		DigitalMicrograph::Result("Could not setup OpenCL on this computer, run clinfo to check availability\n");
	}

}

void EWR::OnBnClickedButtonnps() {
	// TODO: Add your control notification handler code here
	DigitalMicrograph::String npslocation;
	DigitalMicrograph::OpenDialog(npslocation);

	//Check file exists
	if (!DigitalMicrograph::DoesFileExist(npslocation)) {
		DigitalMicrograph::OpenAndSetProgressWindow("error", "loading", "file");
		return;
	}

	NPSImage = DigitalMicrograph::NewImageFromFile(npslocation);
	GotNPS = true;


	// These don't seem to be working to display the text...
	//	std::string mtfpath = mtflocation;
	//m_mtfpath = mtfpath.c_str();

	//m_mtflocation.SetWindowText(mtfpath.c_str());
	//UpdateData(TRUE);

	// Check it is at least loading the file
	NPSImage.GetOrCreateImageDocument().Show();
}


void EWR::OnBnClickedButtonmtf() {
	// TODO: Add your control notification handler code here
	DigitalMicrograph::String mtflocation;
	DigitalMicrograph::OpenDialog(mtflocation);

	//Check file exists
	if (!DigitalMicrograph::DoesFileExist(mtflocation)) {
		DigitalMicrograph::OpenAndSetProgressWindow("error", "loading", "file");
		return;
	}

	MTFImage = DigitalMicrograph::NewImageFromFile(mtflocation);
	GotMTF = true;


	// These don't seem to be working to display the text...
	//	std::string mtfpath = mtflocation;
	//m_mtfpath = mtfpath.c_str();

	//m_mtflocation.SetWindowText(mtfpath.c_str());
	//UpdateData(TRUE);

	// Check it is at least loading the file
	MTFImage.GetOrCreateImageDocument().Show();
}
