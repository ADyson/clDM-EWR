#pragma once


#include "TemplateMFCDialogPlugIn.h"
#include "afxwin.h"
#include "standardfunctions.h"
#include "numedit.h"
#include "CL/OpenCl.h"
#include "clAmdFft.h"
#include <complex>
#include "DMPlugInCamera.h"
#include "WorkerClass.h"



// Acquisition dialog
class Acquisition : public CDialog, WorkerClass
{
	DECLARE_DYNAMIC(Acquisition)

public:
	Acquisition(CWnd* pParent = NULL);   // standard constructor
	virtual ~Acquisition();

// Dialog Data
	enum { IDD = IDD_DIALOG2 };



protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
	virtual BOOL PreTranslateMessage(MSG* pMsg);

public:
	// focalstep control
	CNumEdit m_focalStep;
	CString focalstep;
	float floatFocalStep;

	// exposure time control
	CNumEdit m_expTime;
	CString expTime;
	float floatExpTime;

	// number of images control
	CEdit m_numImages;
	CString numImages;
	int intNumImages;

	// number of underfocus images control
	CEdit m_numUnderfocus;
	CString numUnderfocus;
	int intNumUnderfocus;

	// settling time
	CEdit m_settlingTime;
	CString settlingTime;
	float floatSettlingTime;
	
	// binning control
	CComboBox m_binningCombo;
	CString binningText;

	// Voltage Change Mode Switch 1 is VMode, 0 is FMode
	int VoltageModeEnable;

	// progress bar
	CProgressCtrl m_progressBar;
	
	// start acquisition button
	CButton m_acquireButton;
	afx_msg void StartAcquire();
	
	// pause/stop acquisition buton
	CButton m_StopButton;
	afx_msg void StopAcquire();

	// Low Dose Controls

	// Delay before acquisition starts
	CNumEdit m_PreAcqDelay;
	CString PreAcqDelay;
	float floatPreAcqDelay;

	// How many DAC units to move beam... (needs calibration)
	CEdit m_DACShift;
	CString DACShift;
	int intDACShift;

	// Low Dose Mode Switch
	int LowDoseEnable;

	

	// Was going to have a seperate function to call to do Low Dose Acquisition
	//void LowDoseAcquire(float preacq, float postexp);

	// Was going to have a seperate function to call to do Standard Acquisition
	//void StandardAcquire(float preacq, float postexp);

	// From WorkerClass - this function will run in an alternate thread thus not halting UI.
	virtual void DoWork();
	void SingleLowDoseAcquire();
	void SeriesLowDoseAcquire();
	void PreShiftBeam();

	afx_msg void OnBnClickedCheckfocvol();
	afx_msg void OnBnClickedChecklow();
	afx_msg void OnBnClickedSearch();
};
