// DMDialog.cpp : implementation file
//

#include "stdafx.h"
#include "TemplateMFCDialogPlugIn.h"
#include "DMPluginStubs.h"
#include "SeriesFunctions.h"

// To call cl Initialisation Routine
#include "clState.h"




IMPLEMENT_DYNCREATE(SeriesFunctionsDialog, CDialog)
SeriesFunctionsDialog::SeriesFunctionsDialog(CWnd* pParent /*=NULL*/)
	: CDialog(SeriesFunctionsDialog::IDD, pParent)
{
	EditCAL = "0";
	EditV = "0";
	EditDF = "0";
	NormalCheck = false;
	AlignCheck = false;
	EditXDrift = "0";
	EditYDrift = "0";
}

SeriesFunctionsDialog::~SeriesFunctionsDialog()
{
}

void SeriesFunctionsDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST2, ImageList);
	DDX_Control(pDX, IDC_EDITKV, NumEditV);
	DDX_Text(pDX, IDC_EDITKV, EditV);
	DDX_Control(pDX, IDC_EDITDF, NumEditDF);
	DDX_Text(pDX, IDC_EDITDF, EditDF);
	DDX_Control(pDX, IDC_EDITCAL, NumEditCAL);
	DDX_Text(pDX, IDC_EDITCAL, EditCAL);
	DDX_Control(pDX, IDC_EDITXDRIFT, NumEditXDrift);
	DDX_Text(pDX, IDC_EDITXDRIFT, EditXDrift);
	DDX_Control(pDX, IDC_EDITYDRIFT, NumEditYDrift);
	DDX_Text(pDX, IDC_EDITYDRIFT, EditYDrift);
	DDX_Check(pDX, IDC_ALIGN, AlignCheck);
	DDX_Check(pDX, IDC_NORMAL, NormalCheck);
}


BEGIN_MESSAGE_MAP(SeriesFunctionsDialog, CDialog)
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_BN_CLICKED(IDC_GET, &SeriesFunctionsDialog::OnBnClickedGet)
	ON_BN_CLICKED(IDC_DISABLE, &SeriesFunctionsDialog::OnBnClickedDisable)
	ON_BN_CLICKED(IDC_ENABLE, &SeriesFunctionsDialog::OnBnClickedEnable)
	ON_BN_CLICKED(IDC_NORMAL, &SeriesFunctionsDialog::OnBnClickedNormal)
	ON_BN_CLICKED(IDC_ALIGN, &SeriesFunctionsDialog::OnBnClickedAlign)
	ON_BN_CLICKED(IDC_SETKV, &SeriesFunctionsDialog::OnBnClickedSetkv)
	ON_BN_CLICKED(IDC_SETFSTEP, &SeriesFunctionsDialog::OnBnClickedSetfstep)
	ON_BN_CLICKED(IDC_SETCAL, &SeriesFunctionsDialog::OnBnClickedSetcal)
	ON_BN_CLICKED(IDC_SETX, &SeriesFunctionsDialog::OnBnClickedSetx)
END_MESSAGE_MAP()


// SeriesFunctionsDialog message handlers

BOOL SeriesFunctionsDialog::Create(UINT templateID, CWnd* pParentWnd)
{
	// TODO: Add your specialized code here and/or call the base clase
	return CDialog::Create(IDD, pParentWnd);
}



BOOL SeriesFunctionsDialog::OnInitDialog()
{
	CDialog::OnInitDialog();
	// TODO:  Add extra initialization here
	ModifyStyle(0, WS_GROUP | WS_TABSTOP);

	// Add columns to the ListCtrl
	CRect rect;
	ImageList.GetClientRect(&rect);
	ImageList.InsertColumn(0, _T("Image Number"), LVCFMT_LEFT, rect.Width()*0.5);
	ImageList.InsertColumn(1, _T("Include"), LVCFMT_LEFT, rect.Width()-(rect.Width()*0.6));
	// Leave a gap because of scrollbar taking some space, dont want a horizontal scroll to appear.

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}


void SeriesFunctionsDialog::OnPaint()
{
	CDialog::OnPaint();
}

BOOL SeriesFunctionsDialog::OnEraseBkgnd(CDC* pDC)
{
	// TODO: Add your message handler code here and/or call default

	return CDialog::OnEraseBkgnd(pDC);
}

void SeriesFunctionsDialog::PopulateList()
{
	// Clear current list first...
	ImageList.DeleteAllItems();

	//One item per row. Other columns are subitems...

	// Use the LV_ITEM structure to insert the items
	LVITEM lvi;
	CString strItem;

	for (int i = 0; i < SelectedSeries.size(); i++)
	{
		// Insert the first item
		lvi.mask = LVIF_TEXT;
		strItem.Format(_T("Image %i"), i);
		lvi.iItem = i;
		lvi.iSubItem = 0;
		lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
		//lvi.iImage = i % 8;		// There are 8 images in the image list
		ImageList.InsertItem(&lvi);
		// Set subitem 1
		if (SelectedSeries[i].second)
			strItem = "true";
		else strItem = "false";

		lvi.iSubItem = 1;
		lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
		ImageList.SetItem(&lvi);
	}
}


void SeriesFunctionsDialog::OnBnClickedGet()
{
	SelectedSeries.clear();

	try 
	{
		SeriesImage = DigitalMicrograph::GetFrontImage();
	}
	CatchAll("No Front Image");

	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	// Get number of slices in the image
	unsigned long xPix, yPix, zPix;
	SeriesImage.GetDimensionSizes(xPix, yPix, zPix);

	// Ideally would read in whether they already have selections stored...
	for (int i = 0; i < zPix; i++)
	{
		SelectedSeries.push_back(std::make_pair(i, true));
	}

	ReadSelectedTags();
	PopulateList();

	// Get normalization state and adjust dialog
	std::string normalized;
	if (imagetags.DoesTagExist("Focal Series::Normalized"))
	{
		imagetags.GetTagAsString("Focal Series::Normalized", normalized);
		if (normalized == "True" || normalized == "true")
			NormalCheck = TRUE;
		else
			NormalCheck = FALSE;
	}
	else
	{
		NormalCheck = false;
	}

	// Get force aligned state and adjust dialog
	std::string noalign;
	if (imagetags.DoesTagExist("Focal Series::No Align"))
	{
		imagetags.GetTagAsString("Focal Series::No Align", noalign);
		if (noalign == "True" || noalign == "true")
			AlignCheck = TRUE;
		else
			AlignCheck = FALSE;
	}
	else
	{
		AlignCheck = false;
	}

	// Get voltage and adjust dialog
	float voltage;
	if (imagetags.DoesTagExist("Microscope Info:Voltage"))
	{
		imagetags.GetTagAsFloat("Microscope Info:Voltage", &voltage);
		EditV = Lex(voltage).c_str();
	}
	else
	{
		EditV = "0";
	}

	// Get focalstep and adjust dialog
	float fStep;
	if (imagetags.DoesTagExist("Focal Series::Adjusted focalstep"))
	{
		imagetags.GetTagAsFloat("Focal Series:Adjusted focalstep", &fStep);
		EditDF = Lex(fStep).c_str();
	}
	else
	{
		EditDF = "0";
	}

	// Get manual alignment and adjust dialog
	float xdrift;
	if (imagetags.DoesTagExist("Focal Series:Manual X Drift"))
	{
		imagetags.GetTagAsFloat("Focal Series:Manual X Drift", &xdrift);
		EditXDrift = Lex(xdrift).c_str();
	}
	else
	{
		EditXDrift = "0";
	}

	float ydrift;
	if (imagetags.DoesTagExist("Focal Series:Manual Y Drift"))
	{
		imagetags.GetTagAsFloat("Focal Series:Manual Y Drift", &ydrift);
		EditYDrift = Lex(ydrift).c_str();
	}
	else
	{
		EditYDrift = "0";
	}

	// Get calibration and adjust dialog
	float pixelscale = SeriesImage.GetDimensionScale(0);
	EditCAL = Lex(pixelscale).c_str();

	// For sending changes to dialog use FALSE
	UpdateData(FALSE);

	
}


void SeriesFunctionsDialog::OnBnClickedDisable()
{
	POSITION p = ImageList.GetFirstSelectedItemPosition();
	while (p)
	{
		int nSelected = ImageList.GetNextSelectedItem(p);
		SelectedSeries[nSelected].second = false;
	}

	
	PopulateList();
	WriteSelectedTags();
}


void SeriesFunctionsDialog::OnBnClickedEnable()
{
	POSITION p = ImageList.GetFirstSelectedItemPosition();
	while (p)
	{
		int nSelected = ImageList.GetNextSelectedItem(p);
		SelectedSeries[nSelected].second = true;
	}

	PopulateList();
	WriteSelectedTags();
}

void SeriesFunctionsDialog::ReadSelectedTags()
{
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	for (int i = 1; i <= SelectedSeries.size(); i++)
	{
		bool t = true;
		if (imagetags.DoesTagExist("Focal Series::Selected::" + Lex(i - 1)))
		{
			imagetags.GetTagAsBoolean("Focal Series::Selected::" + Lex(i - 1), &t);
			SelectedSeries[i - 1].second = t;
		}
	}
}

void SeriesFunctionsDialog::WriteSelectedTags()
{
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	for (int i = 1; i <= SelectedSeries.size(); i++)
	{
		std::string tagpath = "Focal Series::Selected::" + Lex(i - 1);
		imagetags.SetTagAsBoolean(tagpath.c_str(), SelectedSeries[i - 1].second);
	}
}

void SeriesFunctionsDialog::OnBnClickedNormal()
{
	UpdateData(TRUE);
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	bool t = NormalCheck;
	imagetags.SetTagAsBoolean("Focal Series::Normalized", t);
}

void SeriesFunctionsDialog::OnBnClickedAlign()
{
	UpdateData(TRUE);
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	bool t = AlignCheck;
	imagetags.SetTagAsBoolean("Focal Series::No Align", t);
}


void SeriesFunctionsDialog::OnBnClickedSetkv()
{
	UpdateData(TRUE);
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	imagetags.SetTagAsString("Microscope Info::Voltage", EditV);
}


void SeriesFunctionsDialog::OnBnClickedSetfstep()
{
	UpdateData(TRUE);
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	imagetags.SetTagAsString("Focal Series::Adjusted focalstep", EditDF);
}


void SeriesFunctionsDialog::OnBnClickedSetcal()
{
	UpdateData(TRUE);
	SeriesImage.SetDimensionScale(0, boost::lexical_cast<float>(EditCAL));
	SeriesImage.SetDimensionScale(1, boost::lexical_cast<float>(EditCAL));
}


void SeriesFunctionsDialog::OnBnClickedSetx()
{
	UpdateData(TRUE);
	DigitalMicrograph::TagGroup imagetags = SeriesImage.GetTagGroup();

	imagetags.SetTagAsString("Focal Series::Manual X Drift", EditXDrift);
	imagetags.SetTagAsString("Focal Series::Manual Y Drift", EditYDrift);
}
