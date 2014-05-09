// OptionsDialog.cpp : implementation file
//

#include "stdafx.h"
#include "OptionsDialog.h"
#include "afxdialogex.h"
#include "FTSR.h"


// COptionsDialog dialog

IMPLEMENT_DYNAMIC(COptionsDialog, CDialog)

COptionsDialog::COptionsDialog(CWnd* pParent /*=NULL*/)
	: CDialog(COptionsDialog::IDD, pParent)
	, m_SearchPercentage_Text(_T("30"))
	, m_Kmax_Text(_T("4"))
	, m_Reference_Text(_T("1"))
	, m_SearchSteps_Text(_T("5"))
	, m_MaxDrift_Text(_T(""))
	, m_SNR_Text(_T("0.2"))
{

}

COptionsDialog::COptionsDialog(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf, CWnd* pParent)
	: CDialog(COptionsDialog::IDD, pParent)
{
	// Set dialog parameters during construction
	mParent = pParent;
	m_SearchPercentage_Text = searchpercentage.c_str();
	m_Kmax_Text = pcfkmax.c_str();
	m_Reference_Text = reference.c_str();
	m_SearchSteps_Text = pcftrials.c_str();
	DetermineDF = determinedf;
	ReconPCF = reconpcf;

}


COptionsDialog::COptionsDialog(std::string pcfkmax, std::string pcftrials, std::string reference, std::string searchpercentage, int determinedf, int reconpcf, int magrotfix, std::string magscale, std::string rotscale, std::string maxdrift,std::string snr, int mi, CWnd* pParent)
	: CDialog(COptionsDialog::IDD, pParent)
{
	// Set dialog parameters during construction
	mParent = pParent;
	m_SearchPercentage_Text = searchpercentage.c_str();
	m_Kmax_Text = pcfkmax.c_str();
	m_Reference_Text = reference.c_str();
	m_SearchSteps_Text = pcftrials.c_str();
	DetermineDF = determinedf;
	ReconPCF = reconpcf;
	MagRotFix = magrotfix;
	MI = mi;
	m_MagScale_Text = magscale.c_str();
	m_RotScale_Text = rotscale.c_str();
	m_MaxDrift_Text = maxdrift.c_str();
	m_SNR_Text = snr.c_str();

}

COptionsDialog::~COptionsDialog()
{
}

void COptionsDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX,IDC_CHECK1,DetermineDF);
	DDX_Check(pDX,IDC_CHECK2,ReconPCF);
	DDX_Check(pDX,IDC_CHECKMI,MI);
	DDX_Text(pDX,IDC_EDIT1,m_SearchSteps_Text);
	DDX_Text(pDX,IDC_EDIT8,m_Reference_Text);
	DDX_Control(pDX,IDC_EDIT9,m_Kmax);
	DDX_Text(pDX,IDC_EDIT9,m_Kmax_Text);
	DDX_Control(pDX,IDC_EDIT2,m_SearchPercentage);
	DDX_Text(pDX,IDC_EDIT2,m_SearchPercentage_Text);
	DDX_Control(pDX,IDC_EDIT11,m_MagScale);
	DDX_Text(pDX,IDC_EDIT11,m_MagScale_Text);
	DDX_Control(pDX,IDC_EDIT10,m_RotScale);
	DDX_Text(pDX,IDC_EDIT10,m_RotScale_Text);
	DDX_Check(pDX,IDC_CHECK3,MagRotFix);
	DDX_Control(pDX,IDC_EDITMAXDRIFT,m_MaxDrift);
	DDX_Text(pDX,IDC_EDITMAXDRIFT,m_MaxDrift_Text);
	DDX_Control(pDX,IDC_EDITSNR,m_SNR);
	DDX_Text(pDX,IDC_EDITSNR,m_SNR_Text);
}

void COptionsDialog::OnOK() 
{
   // TODO: Add extra validation here
   
   // Ensure that your UI got the necessary input 
   // from the user before closing the dialog. The 
   // default OnOK will close this.
   //if ( m_nMyValue == 0 ) // Is a particular field still empty?
   //{
    //  AfxMessageBox("Please enter a value for MyValue");
     // return; // Inform the user that he can't close the dialog without
              // entering the necessary values and don't close the 
              // dialog.
   //}

	UpdateData(TRUE);


	//((EWR*)mParent)->UpdateOptions((LPCSTR)m_Kmax_Text, (LPCSTR)m_SearchSteps_Text, (LPCSTR)m_Reference_Text, (LPCSTR) m_SearchPercentage_Text, DetermineDF, ReconPCF);
	((EWR*)mParent)->UpdateOptions2((LPCSTR)m_Kmax_Text, (LPCSTR)m_SearchSteps_Text, (LPCSTR)m_Reference_Text, (LPCSTR) m_SearchPercentage_Text, DetermineDF, ReconPCF,MagRotFix,(LPCSTR)m_MagScale_Text,(LPCSTR)m_RotScale_Text,(LPCSTR)m_MaxDrift_Text,(LPCSTR)m_SNR_Text,MI);

   CDialog::OnOK(); // This will close the dialog and DoModal will return.
}


BEGIN_MESSAGE_MAP(COptionsDialog, CDialog)
END_MESSAGE_MAP()


// COptionsDialog message handlers
