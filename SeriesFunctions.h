#pragma once
#include "amsEdit.h"
#include "resource.h"
#include "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\atlmfc\include\afxwin.h"
#include "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\atlmfc\include\afxcmn.h"


// CDMDialog dialog
class SeriesFunctionsDialog : public CDialog
{
	DECLARE_DYNCREATE(SeriesFunctionsDialog)

public:
	SeriesFunctionsDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~SeriesFunctionsDialog();

// Dialog Data
	enum { IDD = IDD_SERIESDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	void PopulateList();
	void ReadSelectedTags();
	void WriteSelectedTags();

	DigitalMicrograph::Image SeriesImage;
	std::vector<std::pair<int, bool>> SelectedSeries;

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL Create(UINT templateID, CWnd* pParentWnd = NULL);
	virtual BOOL OnInitDialog();
	
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	CListCtrl ImageList;
	afx_msg void OnBnClickedGet();
	afx_msg void OnBnClickedDisable();
	afx_msg void OnBnClickedEnable();
	CString EditV;
	CString EditDF;
	CString EditCAL;
	CAMSNumericEdit NumEditV;
	CAMSNumericEdit NumEditDF;
	CAMSNumericEdit NumEditCAL;
	BOOL NormalCheck;
	BOOL AlignCheck;
	CAMSNumericEdit NumEditXDrift;
	CAMSNumericEdit NumEditYDrift;
	CString EditXDrift;
	CString EditYDrift;
	afx_msg void OnBnClickedNormal();
	afx_msg void OnBnClickedAlign();
	afx_msg void OnBnClickedSetkv();
	afx_msg void OnBnClickedSetfstep();
	afx_msg void OnBnClickedSetcal();
	afx_msg void OnBnClickedSetx();
};
