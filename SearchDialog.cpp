// SearchDialog.cpp : implementation file
//

#include "stdafx.h"
#include "SearchDialog.h"
//#include "afxdialogex.h"
#include "DMScriptFunctions.h"
#include "Utility.h"
#include "boost\lexical_cast.hpp"
#include "Acquisition.h"


// SearchDialog dialog

IMPLEMENT_DYNAMIC(SearchDialog, CDialog)

// Default cosntructor set all values to zero.
SearchDialog::SearchDialog(CWnd* pParent /*=NULL*/)
	: CDialog(SearchDialog::IDD, pParent)
	, s_premovex(_T("0"))
	, s_premovey(_T("0"))
	, s_foc(_T("0"))
	, s_premovedelay(_T("0"))
{
	AcqParent = (Acquisition*) pParent;
}

SearchDialog::~SearchDialog()
{
}

void SearchDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PREMOVEX, m_premovex);
	DDX_Control(pDX, IDC_PREMOVEY, m_premovey);
	DDX_Text(pDX, IDC_PREMOVEX, s_premovex);
	DDX_Text(pDX, IDC_PREMOVEY, s_premovey);
	DDX_Control(pDX, IDC_FOCUS, m_foc);
	DDX_Text(pDX, IDC_FOCUS, s_foc);
	DDX_Control(pDX, IDC_SETFOCUS, b_SetFocus);
	DDX_Control(pDX, IDC_PREMOVEDELAY, m_premovedelay);
	DDX_Text(pDX, IDC_PREMOVEDELAY, s_premovedelay);
}


BEGIN_MESSAGE_MAP(SearchDialog, CDialog)
	ON_BN_CLICKED(IDC_SEARCHACQUIRE, &SearchDialog::OnBnClickedSearchacquire)
	ON_BN_CLICKED(IDC_SEARCHACQUIRESERIES, &SearchDialog::OnBnClickedSearchacquireseries)
END_MESSAGE_MAP()


// SearchDialog message handlers


void SearchDialog::OnBnClickedSearchacquire()
{
	UpdateData(true);
	// Want to move stage to correct position, take an image low dose style, then move stage back into position.
	// If selected, want to determine and correct focus in new position first. ( 2 Methods....)

	// Get current stage positions;

	double StagePositionX = 0;
	double StagePositionY = 0;

	StagePositionX = DMScript::EMGetStageX();
	StagePositionY = DMScript::EMGetStageY();

	Utility::SetResultWindow("Current Stage at X: " + boost::lexical_cast<std::string>(StagePositionX) +" Y: " + boost::lexical_cast<std::string>(StagePositionY));

	// Get dialog values for distance to move;
	double StageMoveX = 0;
	double StageMoveY = 0;

	try
	{
		StageMoveX = boost::lexical_cast<double>(s_premovex);
		StageMoveY = boost::lexical_cast<double>(s_premovey);
	}
	catch(boost::bad_lexical_cast e)
	{
		DigitalMicrograph::Result("Cannot parse stage movement amount \n");
		return;
	}

	// Move Stage to New Position

	double NewStagePositionX = StageMoveX + StagePositionX;
	double NewStagePositionY = StageMoveY + StagePositionY;

	AcqParent->PreShiftBeam();

	// Output movement details
	Utility::SetResultWindow("Moving Stage to X: " + boost::lexical_cast<std::string>(NewStagePositionX)+ " Y: " + boost::lexical_cast<std::string>(NewStagePositionY));


	//DMScript::EMSetStageXY(NewStagePositionX,NewStagePositionY);
	DMScript::EMSetStageX(NewStagePositionX);
	DMScript::EMSetStageY(NewStagePositionY);

	// Subroutine to determine focus from quick images?
	// float focus = DetermineFocus();

	// Wait abit to allow it to settle down
	float premovedelay = boost::lexical_cast<double>(s_premovedelay);

	DigitalMicrograph::Sleep(premovedelay);


	// Now call an acquisition lowdosestyle (Use prompts to get settings for single acquisition?)
	AcqParent->SingleLowDoseAcquire();

	// Return to previous position.

	// Output movement details
	Utility::SetResultWindow("Moving Stage to X: " + boost::lexical_cast<std::string>(StagePositionX)+ " Y: " + boost::lexical_cast<std::string>(StagePositionY));

	//DMScript::EMSetStageXY(StagePositionX,StagePositionY);
	DMScript::EMSetStageX(StagePositionX);
	DMScript::EMSetStageY(StagePositionY);

}

void SearchDialog::OnBnClickedSearchacquireseries()
{
	UpdateData(true);

	// Want to move stage to correct position, take an image low dose style, then move stage back into position.
	// If selected, want to determine and correct focus in new position first. ( 2 Methods....)

	// Get current stage positions;

	double StagePositionX = 0;
	double StagePositionY = 0;

	StagePositionX = DMScript::EMGetStageX();
	StagePositionY = DMScript::EMGetStageY();

	Utility::SetResultWindow("Current Stage at X: " + boost::lexical_cast<std::string>(StagePositionX) +" Y: " + boost::lexical_cast<std::string>(StagePositionY));


	// Get dialog values for distance to move;
	double StageMoveX = 0;
	double StageMoveY = 0;

	try
	{
		StageMoveX = boost::lexical_cast<double>(s_premovex);
		StageMoveY = boost::lexical_cast<double>(s_premovey);
	}
	catch(boost::bad_lexical_cast e)
	{
		DigitalMicrograph::Result("Cannot parse stage movement amount \n");
		return;
	}

	// Move Stage to New Position

	double NewStagePositionX = StageMoveX + StagePositionX;
	double NewStagePositionY = StageMoveY + StagePositionY;

	// Move beam out if low dose

	AcqParent->PreShiftBeam();

	// Output movement details
	Utility::SetProgressWindow("Moving Stage to ","X: " + boost::lexical_cast<std::string>(NewStagePositionX),"Y: " + boost::lexical_cast<std::string>(NewStagePositionY));

	bool ResultOK = false;
	std::string prompt = "Move stage to "+boost::lexical_cast<std::string>(NewStagePositionX) + " , " +boost::lexical_cast<std::string>(NewStagePositionY);
	DigitalMicrograph::GetBoolean(prompt.c_str(),false,&ResultOK);

	if(ResultOK) {
		DMScript::EMSetStageX(NewStagePositionX);
		DMScript::EMSetStageY(NewStagePositionY);
	}


	// Subroutine to determine focus from quick images?
	// float focus = DetermineFocus();

	// Wait abit to allow it to settle down
	float premovedelay = boost::lexical_cast<double>(s_premovedelay);

	DigitalMicrograph::Sleep(premovedelay);

	// Now call an acquisition lowdosestyle (Use prompts to get settings for single acquisition?)
	AcqParent->SeriesLowDoseAcquire();

	// Return to previous position.

	// Output movement details
	Utility::SetProgressWindow("Moving Stage to ","X: " + boost::lexical_cast<std::string>(StagePositionX),"Y: " + boost::lexical_cast<std::string>(StagePositionY));
	
	bool ResultOK2 = false;
	std::string prompt2 = "Move stage to "+boost::lexical_cast<std::string>(StagePositionX) + " , " +boost::lexical_cast<std::string>(StagePositionY);
	DigitalMicrograph::GetBoolean(prompt2.c_str(),false,&ResultOK2);

	if(ResultOK2) {
		DMScript::EMSetStageX(StagePositionX);
		DMScript::EMSetStageY(StagePositionY);
	}


}
