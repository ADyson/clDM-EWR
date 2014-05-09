#pragma once
#include "resource.h"
#include "numedit.h"
#include "amsEdit.h"
#include "c:\program files (x86)\microsoft visual studio 9.0\vc\atlmfc\include\afxwin.h"
#include "CL/OpenCl.h"
#include "clAmdFft.h"
#include <complex>
#include "clKernel.h"
#include "clFourier.h"
//#include "PCPCFFunction.h"
#include "WorkerClass.h"
#include "Registration.h"


// EWR dialog
class EWR : public CDialog, public WorkerClass
{
	DECLARE_DYNAMIC(EWR)

public:
	EWR(CWnd* pParent = NULL);   // standard constructor
	virtual ~EWR();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
	virtual BOOL PreTranslateMessage(MSG* pMsg);

public:
	
	void UpdateOptions(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf);
	void UpdateOptions2(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf, int magrotfix, std::string magscale, std::string rotscale, std::string maxdrift, std::string snr, int mi);
	
	//PCPCFOptions options;
	//Abberrations Abb;
	
	// Threaded WorkerClassFunction
	void DoWork();
	
	// CL Stuff
	//bool OpenCLAvailable;
	//cl_context context;
	// Use this to check status after every API call
	//cl_int status;
	//cl_uint numDevices;
	//cl_device_id* devices;
	//clDevice* cldev;
	//clQueue* clq;

	clFourier* FFT;
	Registration PCFWrapper;

	// Recheck tags to see if we are still using same device if its changed clear and re initialise
	afx_msg void ResetOpenCLDevice();
	
	CAMSNumericEdit m_Inf;
	CString m_Inf_Text;

	CAMSNumericEdit m_SearchStart;
	CString m_SearchStart_Text;

	CAMSNumericEdit m_SearchEnd;
	CString m_SearchEnd_Text;

	CString m_SearchSteps_Text;

	CAMSNumericEdit m_Cs;
	CString m_Cs_Text;

	CAMSNumericEdit m_Beta;
	CString m_Beta_Text;

	CAMSNumericEdit m_Delta;
	CString m_Delta_Text;

	std::string maxdrift;

	DigitalMicrograph::Image MTFImage;
	DigitalMicrograph::Image NPSImage;
	bool GotMTF;
	bool GotNPS;

	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedOptions();
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedOptionsbutton2();
	afx_msg void OnBnClickedAdjust();
	afx_msg void OnBnClickedButtonnps();
	afx_msg void OnBnClickedButtonmtf();
};
