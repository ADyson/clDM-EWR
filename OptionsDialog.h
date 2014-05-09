#pragma once
#include "numedit.h"
#include "amsEdit.h"
#include "resource.h"

// COptionsDialog dialog

class COptionsDialog : public CDialog
{
	DECLARE_DYNAMIC(COptionsDialog)

public:
	COptionsDialog(CWnd* pParent = NULL);   // standard constructor
	COptionsDialog(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf, CWnd* pParent = NULL);   // standard constructor
	COptionsDialog(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf, int magrotfix, std::string magscale, std::string rotscale, std::string maxdrift, std::string snr, int mi, CWnd* pParent = NULL);   // 3 new options
	virtual void OnOK();
	virtual ~COptionsDialog();

// Dialog Data
	enum { IDD = IDD_OPTIONSDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

	CWnd* mParent;

	CAMSNumericEdit m_Kmax;
	CString m_Kmax_Text;

	CAMSNumericEdit m_MaxDrift;
	CString m_MaxDrift_Text;

	CString m_SearchSteps_Text;

	CAMSNumericEdit m_SearchPercentage;
	CString m_SearchPercentage_Text;

	CString m_Reference_Text;

	CAMSNumericEdit m_MagScale;
	CString m_MagScale_Text;

	CAMSNumericEdit m_RotScale;
	CString m_RotScale_Text;

	CAMSNumericEdit m_SNR;
	CString m_SNR_Text;

	int DetermineDF;
	int ReconPCF;
	int MagRotFix;
	int MI;
};
