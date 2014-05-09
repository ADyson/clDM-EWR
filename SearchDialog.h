#pragma once
#include "numedit.h"
#include "amsEdit.h"
#include "resource.h"
#include "c:\program files (x86)\microsoft visual studio 9.0\vc\atlmfc\include\afxwin.h"
#include "clState.h"

class Acquisition;

// SearchDialog dialog

class SearchDialog : public CDialog
{
	DECLARE_DYNAMIC(SearchDialog)

public:
	SearchDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~SearchDialog();

// Dialog Data
	enum { IDD = IDD_DIALOG3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CAMSNumericEdit m_premovex;
	CAMSNumericEdit m_premovey;
	CString s_premovex;
	CString s_premovey;
	CAMSNumericEdit m_foc;
	CString s_foc;
	CButton b_SetFocus;
	CAMSNumericEdit m_premovedelay;
	CString s_premovedelay;
	afx_msg void OnBnClickedSearchacquire();

	Acquisition* AcqParent;
	afx_msg void OnBnClickedSearchacquireseries();
};
