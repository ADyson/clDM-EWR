#pragma once


// MyTabCtrl
class MyTabCtrl : public CTabCtrl
{
	DECLARE_DYNAMIC(MyTabCtrl)

public:
	MyTabCtrl();

public:
	//Array to hold the list of dialog boxes/tab pages for CTabCtrl
	int m_DialogID[3];

	//CDialog Array Variable to hold the dialogs 
	CDialog *m_Dialog[3];

	int m_nPageCount;

	//Function to Create the dialog boxes during startup
	void InitDialogs();

	//Function to activate the tab dialog boxes 
	void ActivateTabDialogs();

	virtual ~MyTabCtrl();

protected:

	afx_msg void OnSelchange(NMHDR* pNMHDR, LRESULT* pResult);

	DECLARE_MESSAGE_MAP()
};


