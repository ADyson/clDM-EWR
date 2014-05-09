#pragma once
#include "mytabctrl.h"


// CDMDialog dialog
class CDMDialog : public CDialog
{
	DECLARE_DYNCREATE(CDMDialog)

public:
	CDMDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDMDialog();

// Dialog Data
	enum { IDD = IDD_DMDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL Create(UINT templateID, CWnd* pParentWnd = NULL);
	virtual BOOL OnInitDialog();
	
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	MyTabCtrl m_tbCtrl;

	
};
