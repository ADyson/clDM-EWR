// ManualAlignment.cpp : implementation file
//

#include "stdafx.h"
#include "ManualAlignment.h"
#include "afxdialogex.h"
#include "Utility.h"


// ManualAlignment dialog

IMPLEMENT_DYNAMIC(ManualAlignment, CDialog)

ManualAlignment::ManualAlignment(CWnd* pParent /*=NULL*/)
	: CDialog(ManualAlignment::IDD, pParent)
{
	coarse = false;
}

ManualAlignment::~ManualAlignment()
{
}

void ManualAlignment::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(ManualAlignment, CDialog)
	ON_BN_CLICKED(IDC_BUTTONUP, &ManualAlignment::OnBnClickedButtonup)
	ON_BN_CLICKED(IDC_COARSE, &ManualAlignment::OnBnClickedCoarse)
	ON_BN_CLICKED(IDC_BUTTONRIGHT, &ManualAlignment::OnBnClickedButtonright)
	ON_BN_CLICKED(IDC_BUTTONDOWN, &ManualAlignment::OnBnClickedButtondown)
	ON_BN_CLICKED(IDC_BUTTONLEFT, &ManualAlignment::OnBnClickedButtonleft)
	ON_BN_CLICKED(IDC_BUTTONSTART, &ManualAlignment::OnBnClickedButtonstart)
	ON_BN_CLICKED(IDC_BUTTONSAVE, &ManualAlignment::OnBnClickedButtonsave)
	ON_BN_CLICKED(IDC_BUTTONFINISH, &ManualAlignment::OnBnClickedButtonfinish)
	ON_BN_CLICKED(IDC_BUTTONCLEAR, &ManualAlignment::OnBnClickedButtonclear)
	ON_BN_CLICKED(IDC_BUTTONSHOW, &ManualAlignment::OnBnClickedButtonshow)
END_MESSAGE_MAP()


// ManualAlignment message handlers


void ManualAlignment::OnBnClickedButtonup()
{
	// TODO: Add your control notification handler code here

	img = DigitalMicrograph::GetFrontImage();
	imgdisp = img.GetImageDisplay(0);

	// Find current slice
	long startslice;
	long endslice;

	imgdisp.GetDisplayedLayers(&startslice,&endslice);

	Gatan::uint32 sizeX;
	Gatan::uint32 sizeY;
	Gatan::uint32 sizeZ;

	img.GetDimensionSizes(sizeX,sizeY,sizeZ);


	Gatan::PlugIn::ImageDataLocker imgLocker(img);
	float* imgdata = (float*) imgLocker.get();

	// Copy
	std::vector<float> slicedata(sizeX*sizeY);

	// Rewrite Adjusted Version
	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			int i2 = i;
			int j2 = j + 1 +19*(coarse); // up 1 or 20 depending on coarse...

			if(j2 < 0)
				j2 += sizeY;
			else if(j2 >= sizeY)
				j2-=sizeY;


			slicedata[i + j * sizeX] = imgdata[startslice*sizeX*sizeY + i2 + j2*sizeX];

		}
	}

	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			imgdata[startslice*sizeX*sizeY + i + j*sizeX] = slicedata[i + j * sizeX];
		}
	}

	imgLocker.~ImageDataLocker();

	img.DataChanged();
}


void ManualAlignment::OnBnClickedCoarse()
{
	// TODO: Add your control notification handler code here
	if(coarse)
		coarse = false;
	else
		coarse = true;
}


void ManualAlignment::OnBnClickedButtonright()
{
	img = DigitalMicrograph::GetFrontImage();
	imgdisp = img.GetImageDisplay(0);

	// TODO: Add your control notification handler code here
	// Find current slice
	long startslice;
	long endslice;

	imgdisp.GetDisplayedLayers(&startslice,&endslice);

	Gatan::uint32 sizeX;
	Gatan::uint32 sizeY;
	Gatan::uint32 sizeZ;

	img.GetDimensionSizes(sizeX,sizeY,sizeZ);


	Gatan::PlugIn::ImageDataLocker imgLocker(img);
	float* imgdata = (float*) imgLocker.get();

	// Copy
	std::vector<float> slicedata(sizeX*sizeY);

	// Rewrite Adjusted Version
	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			int i2 = i - 1 -19*(coarse);
			int j2 = j; // up 1 or 20 depending on coarse...

			if(i2 < 0)
				i2 += sizeX;
			else if(i2 >= sizeX)
				i2-=sizeX;

			slicedata[i + j * sizeX] = imgdata[startslice*sizeX*sizeY + i2 + j2*sizeX];

		}
	}

	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			imgdata[startslice*sizeX*sizeY + i + j*sizeX] = slicedata[i + j * sizeX];
		}
	}

	imgLocker.~ImageDataLocker();
	
	img.DataChanged();
}


void ManualAlignment::OnBnClickedButtondown()
{
	img = DigitalMicrograph::GetFrontImage();
	imgdisp = img.GetImageDisplay(0);

	// TODO: Add your control notification handler code here
		// Find current slice
	long startslice;
	long endslice;

	imgdisp.GetDisplayedLayers(&startslice,&endslice);

	Gatan::uint32 sizeX;
	Gatan::uint32 sizeY;
	Gatan::uint32 sizeZ;

	img.GetDimensionSizes(sizeX,sizeY,sizeZ);


	Gatan::PlugIn::ImageDataLocker imgLocker(img);
	float* imgdata = (float*) imgLocker.get();

	// Copy
	std::vector<float> slicedata(sizeX*sizeY);

	// Rewrite Adjusted Version
	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			int i2 = i;
			int j2 = j - 1 -19*(coarse); // up 1 or 20 depending on coarse...

			if(j2 < 0)
				j2 += sizeY;
			else if(j2 >= sizeY)
				j2-=sizeY;


			slicedata[i + j * sizeX] = imgdata[startslice*sizeX*sizeY + i2 + j2*sizeX];

		}
	}

	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			imgdata[startslice*sizeX*sizeY + i + j*sizeX] = slicedata[i + j * sizeX];
		}
	}

	imgLocker.~ImageDataLocker();
	
	img.DataChanged();

}


void ManualAlignment::OnBnClickedButtonleft()
{
	img = DigitalMicrograph::GetFrontImage();
	imgdisp = img.GetImageDisplay(0);

	// TODO: Add your control notification handler code here
		// Find current slice
	long startslice;
	long endslice;

	imgdisp.GetDisplayedLayers(&startslice,&endslice);

	Gatan::uint32 sizeX;
	Gatan::uint32 sizeY;
	Gatan::uint32 sizeZ;

	img.GetDimensionSizes(sizeX,sizeY,sizeZ);


	Gatan::PlugIn::ImageDataLocker imgLocker(img);
	float* imgdata = (float*) imgLocker.get();

	// Copy
	std::vector<float> slicedata(sizeX*sizeY);

	// Rewrite Adjusted Version
	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			int i2 = i + 1 +19*(coarse);
			int j2 = j ; // up 1 or 20 depending on coarse...

			if(i2 < 0)
				i2 += sizeX;
			else if(i2 >= sizeX)
				i2-=sizeX;

			slicedata[i + j * sizeX] = imgdata[startslice*sizeX*sizeY + i2 + j2*sizeX];

		}
	}

	for(int i = 0 ; i < sizeX ; i++)
	{
		for(int j = 0 ; j < sizeY ; j++)
		{
			imgdata[startslice*sizeX*sizeY + i + j*sizeX] = slicedata[i + j * sizeX];
		}
	}

	imgLocker.~ImageDataLocker();
	
	img.DataChanged();
}


void ManualAlignment::OnBnClickedButtonstart()
{
	// TODO: Add your control notification handler code here

	// Initialise a vector with the length of the current image stack to hold the roi positions and whether they have		   been recorded

	DigitalMicrograph::Image StackImage = DigitalMicrograph::GetFrontImage();
	Gatan::uint32 xDim;
	Gatan::uint32 yDim;
	Gatan::uint32 zDim;

	// Get stack size
	StackImage.GetDimensionSizes(xDim,yDim,zDim);

	// Resize all vectors
	ROITop.resize(zDim);
	ROIBottom.resize(zDim);
	ROILeft.resize(zDim);
	ROIRight.resize(zDim);
	ROIRec.resize(zDim);

	// Ensure all bools are set to false
	for(int i = 0 ; i < zDim ; i++)
	{
		ROIRec[i] = false;
	}

	// Finished allocation

}


void ManualAlignment::OnBnClickedButtonsave()
{
	DigitalMicrograph::Image StackImage = DigitalMicrograph::GetFrontImage();
	Gatan::uint32 xDim;
	Gatan::uint32 yDim;
	Gatan::uint32 zDim;

	// Get stack size
	StackImage.GetDimensionSizes(xDim,yDim,zDim);

	// Get position of current image within the stack
	long start;
	long end;
	StackImage.GetImageDisplay(0).GetDisplayedLayers(&start,&end);

	if(end - start != 0)
	{
		DigitalMicrograph::Result("Ensure image is displaying one slice at a time \n");
		return;
	}

	DigitalMicrograph::ROI RoI = StackImage.GetImageDisplay(0).GetROI(0);

	// Get Region Positions
	float top;
	float left;
	float bottom;
	float right;

	// Check its a rectangle ROI
	if(!RoI.IsRectangle())
	{
		DigitalMicrograph::OpenAndSetProgressWindow("ROI is not a rectangle","","");
		return;
	}
	

	RoI.GetRectangle(&top,&left,&bottom,&right);

	// Do Some Checks
	int iTop = boost::numeric_cast<int>(top);
	int iLeft = boost::numeric_cast<int>(left);
	int iBottom = boost::numeric_cast<int>(bottom);
	int iRight = boost::numeric_cast<int>(right);
	

	// Save to vectors
	ROILeft[start] = iLeft;
	ROIRight[start] = iRight;
	ROITop[start] = iTop;
	ROIBottom[start] = iBottom;
	ROIRec[start] = true;

}


void ManualAlignment::OnBnClickedButtonfinish()
{
	// Run this code in background thread
	Start();
}


void ManualAlignment::OnBnClickedButtonclear()
{
	// Here we need to reset the vectors so we can try again.

	ROILeft.clear();
	ROIRight.clear();
	ROITop.clear();
	ROIBottom.clear();
	ROIRec.clear();
}

void ManualAlignment::DoWork()
{
	// Background Threaded Function

	// Create a new stack image with the correct size.
	DigitalMicrograph::Image StackImage = DigitalMicrograph::GetFrontImage();
	Gatan::uint32 xDim;
	Gatan::uint32 yDim;
	Gatan::uint32 zDim;

	// Get stack size
	StackImage.GetDimensionSizes(xDim,yDim,zDim);

	// Check that all ROI's are accounted for and of the same dimensions

	int initialsizex = 0;
	int initialsizey = 0;

	for(int i = 0 ; i < zDim ; i++)
	{
		if(!ROIRec[i])
		{
			// TODO - be more specific - put in a class
			DigitalMicrograph::Result("ROI Missing for one or more images \n");
			return;
		}

		// Set ROI size for first image.
		if(i==0)
		{
			initialsizex = ROIRight[i] - ROILeft[i];
			initialsizey = ROIBottom[i] - ROITop[i];
		}

		// Check all ROIs are matching in size with the initial one.
		if(ROIRight[i] - ROILeft[i] != initialsizex && ROIBottom[i] - ROITop[i] != initialsizey)
		{
			DigitalMicrograph::Result("ROI's not of equal size \n");
			return;
		}
	}

	// Now prepare an image to hold all all ROIs

	DigitalMicrograph::Image AlignedStack = DigitalMicrograph::RealImage("Aligned Stack",4,initialsizex,initialsizey,zDim);

	Gatan::PlugIn::ImageDataLocker AlignedLocker(AlignedStack);
	float* AlignedData = (float*) AlignedLocker.get();

	Gatan::PlugIn::ImageDataLocker StackLocker(StackImage);
	float* StackData = (float*) StackLocker.get();

	for(int image = 0 ; image < zDim; image++)
	{
		for(int j = 0 ; j < initialsizey; j++)
			for(int i = 0; i < initialsizex; i++)
			{
				AlignedData[Utility::LinearAccess3D(i,j,image,initialsizex,initialsizey)] = StackData[Utility::LinearAccess3D(i,ROILeft[image],j,ROITop[image],image,xDim,yDim)];
			}
	}

	// Now show the new stack image and copy relevant calibration data...

	AlignedStack.CopyCalibrationFrom(StackImage);
	DigitalMicrograph::TagGroup stacktags = StackImage.GetTagGroup();
	DigitalMicrograph::TagGroup alignedtags = AlignedStack.GetTagGroup();

	DigitalMicrograph::TagGroupCopyTagsFrom(alignedtags,stacktags);

	alignedtags.SetTagAsString("Focal Series:No Align","True");
	AlignedStack.GetOrCreateImageDocument().Show();
	
}

void ManualAlignment::OnBnClickedButtonshow()
{
	// Move ROI to its recorded position

	DigitalMicrograph::Image StackImage = DigitalMicrograph::GetFrontImage();

	// Get position of current image within the stack
	long start;
	long end;
	StackImage.GetImageDisplay(0).GetDisplayedLayers(&start,&end);

	if(end - start != 0)
	{
		DigitalMicrograph::Result("Ensure image is displaying one slice at a time \n");
		return;
	}

	DigitalMicrograph::ROI RoI = StackImage.GetImageDisplay(0).GetROI(0);

	// Check its a rectangle ROI
	if(!RoI.IsRectangle())
	{
		DigitalMicrograph::OpenAndSetProgressWindow("ROI is not a rectangle","","");
		return;
	}

	// Check if position has been recorded.
	if(ROIRec[start])
	{
			RoI.SetRectangle(ROITop[start],ROILeft[start],ROIBottom[start],ROIRight[start]);
	}
	else
	{
		DigitalMicrograph::Result("No ROI position recorded for this slice");
		return;
	}

}
