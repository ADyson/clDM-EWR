#pragma once
#include "Resource.h"
#include "WorkerClass.h"


// ManualAlignment dialog

class ManualAlignment : public CDialog, public WorkerClass
{
	DECLARE_DYNAMIC(ManualAlignment)

public:
	ManualAlignment(CWnd* pParent = NULL);   // standard constructor
	virtual ~ManualAlignment();

// Dialog Data
	enum { IDD = IDD_MANUAL };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonup();

	DigitalMicrograph::Image img;
	DigitalMicrograph::ImageDisplay imgdisp;
	bool coarse;
	afx_msg void OnBnClickedCoarse();
	afx_msg void OnBnClickedButtonright();
	afx_msg void OnBnClickedButtondown();
	afx_msg void OnBnClickedButtonleft();
	afx_msg void OnBnClickedButtonstart();
	afx_msg void OnBnClickedButtonsave();
	afx_msg void OnBnClickedButtonfinish();
	afx_msg void OnBnClickedButtonclear();

	std::vector<int> ROITop;
	std::vector<int> ROIBottom;
	std::vector<int> ROILeft;
	std::vector<int> ROIRight;
	std::vector<bool> ROIRec;



	void DoWork();
	afx_msg void OnBnClickedButtonshow();
};
