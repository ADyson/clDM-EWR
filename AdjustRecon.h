#pragma once
#include "Resource.h"
#include "c:\program files (x86)\microsoft visual studio 9.0\vc\atlmfc\include\afxwin.h"
#include "clKernel.h"
#include "clFourier.h"
#include "amsEdit.h"


// AdjustRecon dialog

class AdjustRecon : public CDialog
{
	DECLARE_DYNAMIC(AdjustRecon)

public:
	AdjustRecon(CWnd* pParent = NULL);   // standard constructor
	virtual ~AdjustRecon();

// Dialog Data
	enum { IDD = IDD_ADJUSTRECON };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

	void PostNcDestroy();
public:
	afx_msg void OnBnClickedAdjust();

	CAMSNumericEdit m_C1;
	CAMSNumericEdit m_C3;
	CAMSNumericEdit m_A1r;
	CAMSNumericEdit m_A1i;
	CAMSNumericEdit m_A2r;
	CAMSNumericEdit m_A2i;
	CAMSNumericEdit m_B2r;
	CAMSNumericEdit m_B2i;
	CAMSNumericEdit m_A3r;
	CAMSNumericEdit m_A3i;
	CAMSNumericEdit m_S3r;
	CAMSNumericEdit m_S3i;

	CString C1;
	CString C3;
	CString A1r;
	CString A1i;
	CString A2r;
	CString A2i;
	CString B2r;
	CString B2i;
	CString A3r;
	CString A3i;
	CString S3r;
	CString S3i;

	void* mParent;
	cl_int status;
};
