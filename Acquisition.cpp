// Acquisition.cpp : implementation file

// Enable double support


#include "stdafx.h"
#include "Acquisition.h"
#include "NumEdit.h"
#include "boost/lexical_cast.hpp"
#include "boost/math/special_functions/round.hpp"
#include "SearchDialog.h"
#include "DMScriptFunctions.h"

// Extra functions required for operation
std::string FormattedMag(float mag) {
	if ( mag < 1000000 ) {
		std::stringstream out;
		out << mag / 1000;
		out << "k" ;
		return out.str();
	}

	else {
		std::stringstream out;
		out << mag / 1000000;
		out << "M" ;
		return out.str();
	};
}


std::string FormattedMagED(float mag) {
	std::stringstream out;
	out << "ED";
	out << mag / 10;
	out << "cm" ;
	return out.str();
}

void AddFocalTags(DigitalMicrograph::Image AcquiredImage, int i, float df, float focalstep, int intNumImages) {
	DigitalMicrograph::TagGroup imagetags = AcquiredImage.GetTagGroup();

	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Focal Series:Image Number", boost::lexical_cast<std::string>(i));
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Focal Series:Series Length",
	        boost::lexical_cast<std::string>(intNumImages));
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Focal Series:Defocus", boost::lexical_cast<std::string>(df));
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Focal Series:Adjusted focalstep",
	        boost::lexical_cast<std::string>(focalstep));

}

// Note this function also properly overwrites the microscope info operator/specimen tags.
void AddAberrationTags(DigitalMicrograph::Image imagein) {
	DigitalMicrograph::TagGroup pers;
	DigitalMicrograph::String a2mag = "";
	DigitalMicrograph::String a2ang = "";
	DigitalMicrograph::String b2mag = "";
	DigitalMicrograph::String b2ang = "";
	DigitalMicrograph::String a3mag = "";
	DigitalMicrograph::String a3ang = "";
	DigitalMicrograph::String c3mag = "";
	DigitalMicrograph::String s3mag = "";
	DigitalMicrograph::String s3ang = "";
	DigitalMicrograph::String b3mag = "";
	DigitalMicrograph::String b3ang = "";
	DigitalMicrograph::String a4mag = "";
	DigitalMicrograph::String a4ang = "";
	DigitalMicrograph::String s4mag = "";
	DigitalMicrograph::String s4ang = "";
	DigitalMicrograph::String a5mag = "";
	DigitalMicrograph::String a5ang = "";
	DigitalMicrograph::String s5mag = "";
	DigitalMicrograph::String s5ang = "";

	pers = DigitalMicrograph::GetPersistentTagGroup();

	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A2(Mag)", a2mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A2(Angle)", a2ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:B2(Mag)", b2mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:B2(Angle)", b2ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A3(Mag)", a3mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A3(Angle)", a3ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:C3(Mag)", c3mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S3(Mag)", s3mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S3(Angle)", s3ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:B3(Mag)", b3mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:B3(Angle)", b3ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A4(Mag)", a4mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A4(Angle)", a4ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S4(Mag)", s4mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S4(Angle)", s4ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A5(Mag)", a5mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:A5(Angle)", a5ang);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S5(Mag)", s5mag);
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Aberrations:S5(Angle)", s5ang);

	DigitalMicrograph::TagGroup imagetags;
	DigitalMicrograph::String SpecimenName("");
	DigitalMicrograph::String OperatorName("");

	// Check for Specimen Name

	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Specimen", SpecimenName);

	// Check for Operator Name

	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Operator", OperatorName);
	// Doesn't really matter if blank.

	imagetags = imagein.GetTagGroup();

	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A2(Mag)", a2mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A2(Angle)", a2ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:B2(Mag)", b2mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:B2(Angle)", b2ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A3(Mag)", a3mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A3(Angle)", a3ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:C3(Mag)", c3mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S3(Mag)", s3mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S3(Angle)", s3ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:B3(Mag)", b3mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:B3(Angle)", b3ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A4(Mag)", a4mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A4(Angle)", a4ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S4(Mag)", s4mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S4(Angle)", s4ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A5(Mag)", a5mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:A5(Angle)", a5ang);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S5(Mag)", s5mag);
	DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Aberrations:S5(Angle)", s5ang);


	try {
		DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Microscope Info:Operator", OperatorName);
		DigitalMicrograph::TagGroupSetTagAsString(imagetags, "Microscope Info:Specimen", SpecimenName);

		// These tags are only going to exist already if the other image information is added.
		//DigitalMicrograph::TagGroupSetTagAsString(imagetags,"Microscope Info:Items:[0]:Value",SpecimenName);
		//DigitalMicrograph::TagGroupSetTagAsString(imagetags,"Microscope Info:Items:[1]:Value",OperatorName);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::Result("Can't Add Abberation Tags\n");
	}
}



// Acquisition dialog
IMPLEMENT_DYNAMIC(Acquisition, CDialog)

Acquisition::Acquisition(CWnd * pParent /*=NULL*/)
	: CDialog(Acquisition::IDD, pParent) {
	// Set default starting parameters
	focalstep = "10";
	expTime = "1.0";
	numImages = "20";
	numUnderfocus = "19";
	m_progressBar.SetRange(0, 20);
	binningText = "1";
	DACShift = "3600";
	settlingTime = "0.5";

	// Check for focalstep tag and create if notlo
	DigitalMicrograph::TagGroup persistent = DigitalMicrograph::GetPersistentTagGroup();

	if (!persistent.DoesTagExist("Focal Series:focalstep (nm)")) {
		DigitalMicrograph::TagGroupSetTagAsString(persistent, "Focal Series:focalstep (nm)", "");
	}

	// Should default to zero anyway but must be zero as buttons start unchecked and only logic is alternating on press.
	VoltageModeEnable = 0;
	LowDoseEnable = 0;
}

Acquisition::~Acquisition() {
}

void Acquisition::DoDataExchange(CDataExchange * pDX) {
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT31, m_focalStep);
	DDX_Control(pDX, IDC_EDIT39, m_expTime);
	DDX_Control(pDX, IDC_EDIT32, m_numImages);
	DDX_Control(pDX, IDC_EDIT38, m_numUnderfocus);
	DDX_Control(pDX, IDC_EDITSETTLING, m_settlingTime);
	DDX_Text(pDX, IDC_EDIT31, focalstep);
	DDX_Text(pDX, IDC_EDIT39, expTime);
	DDX_Text(pDX, IDC_EDIT32, numImages);
	DDX_Text(pDX, IDC_EDIT38, numUnderfocus);
	DDX_Text(pDX, IDC_EDITSETTLING, settlingTime);
	DDX_Text(pDX, IDC_COMBO1, binningText);
	DDX_Control(pDX, IDC_COMBO1, m_binningCombo);
	DDX_Control(pDX, IDC_PROGRESS1, m_progressBar);
	DDX_Control(pDX, IDC_BUTTON4, m_acquireButton);
	DDX_Control(pDX, IDC_BUTTON5, m_StopButton);
	DDX_Control(pDX, IDC_EDIT42, m_DACShift);
	DDX_Text(pDX, IDC_EDIT42, DACShift);
}


BEGIN_MESSAGE_MAP(Acquisition, CDialog)
	ON_BN_CLICKED(IDC_BUTTON4, &Acquisition::StartAcquire)
	ON_BN_CLICKED(IDC_CHECKFOCVOL, &Acquisition::OnBnClickedCheckfocvol)
	ON_BN_CLICKED(IDC_CHECKLOW, &Acquisition::OnBnClickedChecklow)
	ON_BN_CLICKED(IDC_SEARCH, &Acquisition::OnBnClickedSearch)
END_MESSAGE_MAP()

BOOL Acquisition::PreTranslateMessage(MSG * pMsg) {
	if (pMsg->message == WM_KEYDOWN) {
		if ((pMsg->wParam == VK_RETURN) || (pMsg->wParam == VK_ESCAPE))
			pMsg->wParam = VK_TAB;
	}

	return CDialog::PreTranslateMessage(pMsg);
}


void Acquisition::StartAcquire() {
	// Don't run two at once - that would be bad
	if (GetStatus() != Working) {
		Start();
	}

	//Enable stop button
	// Disable this button
	GetDlgItem(IDC_BUTTON5)->EnableWindow(TRUE);
}

void Acquisition::DoWork() {
	// Retrieve Values from Dialog ready to start imaging.
	UpdateData(TRUE);

	// Cancel any currently updating live views... didn't seem to work.
	//DigitalMicrograph::StopAllAcquisitionDaemons();
	Gatan::Camera::StopCurrentCameraViewer(false); // Stop it but dont close it...

	// Get size of shift to be put into low dose mode.
	intDACShift = boost::lexical_cast<int>(DACShift);

	// Get settling time
	floatSettlingTime = boost::lexical_cast<float>(settlingTime);
	// Call this before any beam shifting is to be done.
	DigitalMicrograph::EMPrepareShift();

	// If low dose start by moving beam out of position

	// Technically this bit should be in a seperate function, or stage movement should come afterwards, because we want to take a series somewhere we havent been before.
	if (LowDoseEnable == 1) {
		try {
			DigitalMicrograph::EMBeamShift(intDACShift, intDACShift);
		} catch (...) {
			DigitalMicrograph::Result("Unable to shift beam position \n");
			return;
		}
	}

	// First check we have a save location specified.
	DigitalMicrograph::String filepath("");
	DigitalMicrograph::TagGroup pers = DigitalMicrograph::GetPersistentTagGroup();
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:filepath", filepath);

	// Check path exists
	if (!DigitalMicrograph::DoesDirectoryExist( filepath )) {
		DigitalMicrograph::OpenAndSetProgressWindow("Root Directory Not Available", "Or Not Specified", "             :(");
		return;
	}

	// Check for Specimen Name
	DigitalMicrograph::String SpecimenName("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Specimen", SpecimenName);
	// Doesn't really matter if blank.

	// Get Area Number
	DigitalMicrograph::String AreaNumber("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:Area", AreaNumber);

	// Find out todays date
	CTime theTime;
	theTime = CTime::GetCurrentTime();
	CString date = theTime.Format( "%d %B %Y" );
	std::string datestring = ((LPCTSTR)date);



	// Set progress bar to correct size
	intNumImages = boost::lexical_cast<int>(numImages);
	m_progressBar.SetRange(0, intNumImages);

	intNumUnderfocus = boost::lexical_cast<int>(numUnderfocus);
	floatFocalStep = boost::lexical_cast<float>(focalstep);
	floatExpTime = boost::lexical_cast<float>(expTime);

	int binning(0);

	try {
		binning = boost::lexical_cast<int>(binningText);
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	if (!(binning == 1 || binning == 2 || binning == 4)) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	// Check values make sense i.e underfocus < total images
	if ((!(intNumUnderfocus <= intNumImages)) || intNumImages < 1) {
		DigitalMicrograph::OpenAndSetProgressWindow("Number of Images", "Incorrectly Set", "");
		return;
	}

	DigitalMicrograph::String indicatedMag;
	DigitalMicrograph::String operatingVoltage;

	try {
		// Get current microscope settings
		DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Voltage", operatingVoltage);
		DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Indicated Magnification", indicatedMag);
	} catch (...) {
		DigitalMicrograph::Result("Problem getting current microscope settigns \n");
		return;
	}

	float mag = 0.0f;

	try {
		mag = boost::lexical_cast<float>(indicatedMag.get_string());
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Indicated Magnification", "taggroup", "incorrectly set");
		return;
	}

	DigitalMicrograph::String operationMode;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Operation Mode", operationMode);

	// Check for focalstep
	DigitalMicrograph::String microscopeFStep("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:focalstep (nm)", microscopeFStep);
	std::string microscopeFStep2 = microscopeFStep;

	if (microscopeFStep2.empty()) {
		//Need to know this value
		DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:focalstep");
		return;
	}

	// Check for focalstep voltagemode
	DigitalMicrograph::String microscopeFStepVol("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:voltage focalstep (nm)", microscopeFStepVol);
	std::string microscopeFStep3 = microscopeFStepVol;

	if (VoltageModeEnable) {
		if (microscopeFStep3.empty()) {
			//Need to know this value
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:voltage focalstep");
			return;
		}
	}

	float nmperstep = 1.0f;

	if (VoltageModeEnable == 0) {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	} else {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep3);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	}

	// Setup Variables for Microscope (ARM200)
	// Should add nmperdac entry to be applicable to any microscope
	// Maybe read from tag's and create tag if none already...
	Gatan::uint32 xpixels = 0; //Get from CCD anyway
	Gatan::uint32 ypixels = 0; //Get from CCD anyway

	// Get current DM camera
	Gatan::Camera::Camera camera;

	try {
		camera = Gatan::Camera::GetCurrentCamera();
		Gatan::Camera::CCD_GetSize(camera, &xpixels, &ypixels);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("No Camera Detected", "", "");
		return;
	}

	bool inserted = false;

	try {
		inserted = Gatan::Camera::GetCameraInserted(camera);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Couldn't check camera", "status", "");
		return;
	}

	if (inserted != true) {
		DigitalMicrograph::OpenAndSetProgressWindow("Camera not inserted", "", "");
		return;
	}

	// Want gain normalized imaging unless doing post processing yourself
	Gatan::Camera::AcquisitionProcessing processing =
	    Gatan::Camera::kMaxProcessing; //  kGainNormalized/kMaxProcessing (not sure what difference is).
	Gatan::Camera::AcquisitionParameters acqparams;

	try {
		acqparams = Gatan::Camera::CreateAcquisitionParameters(camera, processing, floatExpTime, binning, binning, 0, 0, ypixels,
		            xpixels);
		Gatan::CM::SetDoContinuousReadout(acqparams, true);
		Gatan::CM::SetQualityLevel(acqparams, 1); // Can't remember if fast or slow :D
		Gatan::Camera::Validate_AcquisitionParameters(camera, acqparams);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Problem with acquisition", "parameters", "");
		return;
	}


	// NEW BIT FOR ALTERNATE ACQUISITION
	Gatan::CM::AcquisitionPtr acq = CreateAcquisition( camera, acqparams );
	// Turn into script object for my dms function
	DigitalMicrograph::ScriptObject acqtok =
	    DigitalMicrograph::ScriptObjectProxy<Gatan::Camera::AcquisitionImp, DigitalMicrograph::DMObject>::to_object_token(acq.get());
	Gatan::CM::FrameSetInfoPtr fsi = DMScript::GetFrameSetInfoPtr(acqtok);
	Gatan::Camera::AcquisitionImageSourcePtr acqsource = Gatan::Camera::AcquisitionImageSource::New(acq, fsi, 0);

	// To add calibration information etc..
	Gatan::Camera::AcquisitionImp * Acquisition = acq.get();


	// Check / Create rest of directory structure

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring);

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string());

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string());

	// Make new folder for this focal series
	// First build directory name
	std::string magnification = FormattedMag(mag);
	std::string stdfocalstep = LPCTSTR(focalstep);

	std::string dirname = "focalseries " + stdfocalstep + "nm " + magnification;

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string() + "\\" + dirname))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string() + "\\" + dirname);

	// Count number of files already in directory
	DigitalMicrograph::TagGroup FilesInDir = DigitalMicrograph::GetFilesInDirectory(filepath.get_string() + "\\" + datestring + "\\" +
	        SpecimenName.get_string() + "\\Area " + AreaNumber.get_string() + "\\" + dirname, 2 );
	int Files = DigitalMicrograph::TagGroupCountTags(FilesInDir);
	Files++;

	std::string filedir = filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
	                      AreaNumber.get_string() + "\\" + dirname + "\\" + t_to_string<int>(Files);

	if (!DigitalMicrograph::DoesDirectoryExist(filedir)) {
		DigitalMicrograph::CreateDirectory(filedir);
	}

	// Find image size and setup new images.
	int imwidth = xpixels / binning;
	int imheight = ypixels / binning;

	// Test if its faster without displaying anything - not really
	// Images wont always be realimage - GIF is signed int..
	DigitalMicrograph::Image view = DigitalMicrograph::RealImage("Through Focus View", 4, imwidth, imheight);
	DigitalMicrograph::Image stack = DigitalMicrograph::RealImage("Defocus Stack", 4, imwidth, imheight, intNumImages);
	DigitalMicrograph::ImageDocument viewdoc = view.GetOrCreateImageDocument();
	DigitalMicrograph::ImageDocument stackdoc = stack.GetOrCreateImageDocument();
	Gatan::PlugIn::ImageDataLocker stackLocker(stack);
	float * stackdata = (float *) stackLocker.get();

	// Maybe this is wrong size, factor binning in here.
	viewdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));

	// Hold image information in vectors so can do all saving at the end (faster).
	std::vector<DigitalMicrograph::Image> Images;
	std::vector<std::string> savelocations;

	// Get HighTension info if using voltage mode.
	double startHT;

	if (VoltageModeEnable == 1 ) {
		startHT = DMScript::EMGetHighTensionOffset();
	}


	if (VoltageModeEnable == 0 ) {
		// Go to most underfocus
		if (focalstep != "0") {
			try {
				DigitalMicrograph::EMChangeFocus(-intNumUnderfocus * boost::math::round<int>(floatFocalStep / nmperstep));
			} catch (...) {
				short error;
				long context;
				DigitalMicrograph::GetException(&error, &context);
				DigitalMicrograph::ErrorDialog(error);
				DigitalMicrograph::OpenAndSetProgressWindow("Problem with changing", "focus", "");
				return;
			}
		}
	} else {
		// Go to most underfocus
		if (focalstep != "0") {
			try {
				DMScript::EMSetHighTensionOffset(startHT - intNumUnderfocus * boost::math::round<int>(floatFocalStep / nmperstep));
			} catch (...) {
				short error;
				long context;
				DigitalMicrograph::GetException(&error, &context);
				DigitalMicrograph::ErrorDialog(error);
				DigitalMicrograph::OpenAndSetProgressWindow("Problem with changing", "focus", "");
				return;
			}
		}
	}


	// Start Acquisition
	acqsource->BeginAcquisition();

	for (int i = 1; i <= intNumImages; i++) {
		// Check if we are requested to stop
		if (GetStatus() == Stopping) {
			if (VoltageModeEnable == 0 ) {
				// Fix by how much we have changed focus so far.
				DigitalMicrograph::EMChangeFocus((intNumUnderfocus - (i - 1))*boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT);
			}

			// Beam position will be restored at end anyway

			// Print useful message
			DigitalMicrograph::Result("Acquisition has been aborted \n");

			// Exit loop
			break;
		}


		// Establish current focus including rounding

		int pos = (-intNumUnderfocus + (i - 1)) * boost::math::round<int>(floatFocalStep / nmperstep);
		float df = pos * nmperstep;

		DigitalMicrograph::Image AcquiredImage;

		try {
			AcquiredImage = Gatan::Camera::CreateImageForAcquire(acq, "Acquired Image");
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::Result("Couldn't Create Image for Acquire\n");
			return;
		}

		// Try to add acquisition information
		try {
			Acquisition->Acquire_WriteInfo(AcquiredImage, 1);
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			std::string desc = DigitalMicrograph::GetExceptionDescription();
			DigitalMicrograph::Result(desc + "\n");
			DigitalMicrograph::Result("Couldn't Add Acquisition Tags\n");
			return;
		}

		// Try to acquire an image


		// Bring beam back into position if using low dose
		if (LowDoseEnable == 1) {
			DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
		}

		bool acqprmchanged = false;

		try {
			// Not sure of the effect of max waiting time

			// Returns a bool if a read finishes within call?
			// true means restart if it was already finishing a read, this is incase parameters hadnt changed from last time yet?
			if (!acqsource->AcquireTo(AcquiredImage, true, 0.5f, acqprmchanged)) {
				// Now wait for it to finish again but dont restart if it finishes durign call....
				while (!acqsource->AcquireTo(AcquiredImage, false, 0.5f, acqprmchanged)) {
					// Waiting for read to finish
				}

				// Has read finished now??
				if (LowDoseEnable == 1) {
					// MOVE BEAM GET OUT THE WAY
					DigitalMicrograph::EMBeamShift(intDACShift, intDACShift);
					// GET OUT THE WAY BEAM
					// GET OUT THE WAY
				}
			}
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::OpenAndSetProgressWindow("Couldn't Acquire Image", "", "");
			return;
		}

		Gatan::PlugIn::ImageDataLocker acquireLocker(AcquiredImage);
		Gatan::PlugIn::ImageDataLocker throughfocusLocker(view);

		// Check if camera records float or integer type images.
		if (DigitalMicrograph::IsFloatImage(AcquiredImage.get())) {
			float * acquireddata = (float *) acquireLocker.get();
			float * viewdata = (float *) throughfocusLocker.get();

			for (int j = 0; j < imwidth * imheight; j++) {
				//viewdata[j] = acquireddata[j];
				float val = *(acquireddata++);
				*(viewdata++) = val;
				*(stackdata++) = val;
			}
		}

		else
			if (DigitalMicrograph::IsIntegerDataType(AcquiredImage.get(), 4, true)) {
				int * acquireddata = (int *) acquireLocker.get();
				float * viewdata = (float *) throughfocusLocker.get();

				for (int j = 0; j < imwidth * imheight; j++) {
					//viewdata[j] = acquireddata[j];
					int val = *(acquireddata++);
					*(viewdata++) = val;
					*(stackdata++) = val;
				}
			}

			// If its not float or integer this doesn't work.
			else {
				DigitalMicrograph::OpenAndSetProgressWindow("Image not of expected type", "", "");
				return;
			}

		acquireLocker.~ImageDataLocker();
		throughfocusLocker.~ImageDataLocker();

		view.DataChanged();
		view.SetName(t_to_string(i) + " - " + t_to_string(df) + "nm"); // Rename and Save Image.

		std::string saveloc = filedir + "// " + t_to_string(i) + " " + t_to_string(df) + "nm";
		AcquiredImage.SetName(t_to_string(i) + " " + t_to_string(df) + "nm");

		try {
			AddFocalTags(AcquiredImage, i, df, boost::math::round<int>(floatFocalStep / nmperstep)*nmperstep, intNumImages);
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::Result("Can't Add Focal Series Tags \n");
		}


		AddAberrationTags(AcquiredImage);

		// Add image and savelocation to vector for processing later
		Images.push_back(AcquiredImage);
		savelocations.push_back(saveloc);

		// Don't save here anymore, delay til after all acquisitions
		//DigitalMicrograph::SaveAsGatan(AcquiredImage,saveloc);

		// Don't show image if using continuous updating view
		//DigitalMicrograph::ShowImage(AcquiredImage);
		//AcquiredImage.GetOrCreateImageDocument().Clean();

		// TODO: Add to Stack Image (For Later)

		// Update Progress
		std::string defoc = t_to_string(df);
		std::string st3 = "@ " + defoc + "nm";
		m_progressBar.SetPos(i);
		DigitalMicrograph::OpenAndSetProgressWindow("Image", "Acquired", st3.c_str());

		// Adjust the focus
		if (focalstep != "0") {
			if (VoltageModeEnable == 0 ) {
				// Shift Focus
				DigitalMicrograph::EMChangeFocus(boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT - (intNumUnderfocus + i)*boost::math::round<int>(floatFocalStep / nmperstep));
			}
		}

		// Pause for settling time
		if (floatSettlingTime != 0) {
			DigitalMicrograph::Sleep(floatSettlingTime);
		}

	}

	// If we have been requested to stop then focus is already restored earlier
	if (GetStatus() != Stopping) {
		if (focalstep != "0") {
			if (VoltageModeEnable == 0 ) {
				// Reset Focus (one more than might be expected because it changes again after last image)
				DigitalMicrograph::EMChangeFocus(-(intNumImages - intNumUnderfocus)*boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT);
			}
		}
	}

	// Probably important...
	acqsource->FinishAcquisition();

	// Dont bother saving if we have aborted..
	if (GetStatus() != Stopping) {
		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Now Saving Images...");

		for (int i = 0; i < intNumImages; i++) {
			DigitalMicrograph::SaveAsGatan(Images[i], savelocations[i]);
		}

		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Images Saved.");
	}


	// Get rid of through focus view
	view.GetOrCreateImageDocument().Clean();
	DigitalMicrograph::CloseImage(view);

	// Close stack for editing
	stackLocker.~ImageDataLocker();

	float xscale = Images[0].GetDimensionScale(0);
	float yscale = Images[0].GetDimensionScale(1);

	DigitalMicrograph::TagGroup stacktags = stack.GetTagGroup();
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Focal Series:Adjusted focalstep", t_to_string(focalstep));
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Microscope Info:Voltage", t_to_string(operatingVoltage));

	stack.SetDimensionScale(0, xscale);
	stack.SetDimensionScale(1, yscale);

	stack.DataChanged();

	// Show stack view
	stackdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));
	// Clear vectors and other memory just to be sure

	Images.clear();
	savelocations.clear();

	// If it was low dose we need to return beam to normal position afterwards
	if (LowDoseEnable == 1) {
		DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
	}
}

void Acquisition::OnBnClickedCheckfocvol() {
	// Alternate button text between focus mode and voltage mode.
	if (VoltageModeEnable == 0) {
		VoltageModeEnable = 1;
	}

	else
		if (VoltageModeEnable == 1) {
			VoltageModeEnable = 0;
		}
}


void Acquisition::OnBnClickedChecklow() {
	// Alternate button text between focus mode and voltage mode.
	if (LowDoseEnable == 0) {
		LowDoseEnable = 1;

		std::string output = "LowDoseEnable = " + boost::lexical_cast<std::string>(LowDoseEnable) + "\n";
		//DigitalMicrograph::Result(output);
	} else
		if (LowDoseEnable == 1) {
			LowDoseEnable = 0;

			std::string output = "LowDoseEnable = " + boost::lexical_cast<std::string>(LowDoseEnable) + "\n";
			//DigitalMicrograph::Result(output);
		}
}


void Acquisition::OnBnClickedSearch() {
	// Launch Search Mode Dialog
	SearchDialog * dlg = new SearchDialog(this);
	dlg->Create(SearchDialog::IDD);
	dlg->ShowWindow(SW_SHOW);
}

void Acquisition::SingleLowDoseAcquire() {
// Retrieve Values from Dialog ready to start imaging.
	UpdateData(TRUE);

	// Cancel any currently updating live views... didn't seem to work.
	//DigitalMicrograph::StopAllAcquisitionDaemons();
	Gatan::Camera::StopCurrentCameraViewer(false); // Stop it but dont close it...

	// Get size of shift to be put into low dose mode.
	intDACShift = boost::lexical_cast<int>(DACShift);


	// First check we have a save location specified.
	DigitalMicrograph::String filepath("");
	DigitalMicrograph::TagGroup pers = DigitalMicrograph::GetPersistentTagGroup();
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:filepath", filepath);

	// Check path exists
	if (!DigitalMicrograph::DoesDirectoryExist( filepath )) {
		DigitalMicrograph::OpenAndSetProgressWindow("Root Directory Not Available", "Or Not Specified", "             :(");
		return;
	}

	// Check for Specimen Name
	DigitalMicrograph::String SpecimenName("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Specimen", SpecimenName);
	// Doesn't really matter if blank.

	// Get Area Number
	DigitalMicrograph::String AreaNumber("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:Area", AreaNumber);

	// Find out todays date
	CTime theTime;
	theTime = CTime::GetCurrentTime();
	CString date = theTime.Format( "%d %B %Y" );
	std::string datestring = ((LPCTSTR)date);



	// Set progress bar to correct size
	intNumImages = 1;
	m_progressBar.SetRange(0, intNumImages);

	intNumUnderfocus = 0;
	floatFocalStep = boost::lexical_cast<float>(focalstep);
	floatExpTime = boost::lexical_cast<float>(expTime);

	int binning(0);

	try {
		binning = boost::lexical_cast<int>(binningText);
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	if (!(binning == 1 || binning == 2 || binning == 4)) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	// Get current microscope settings
	DigitalMicrograph::String operatingVoltage;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Voltage", operatingVoltage);
	DigitalMicrograph::String indicatedMag;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Indicated Magnification", indicatedMag);

	float mag = 0.0f;

	try {
		mag = boost::lexical_cast<float>(indicatedMag.get_string());
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Indicated Magnification", "taggroup", "incorrectly set");
		return;
	}

	DigitalMicrograph::String operationMode;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Operation Mode", operationMode);

	// Check for focalstep
	DigitalMicrograph::String microscopeFStep("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:focalstep (nm)", microscopeFStep);
	std::string microscopeFStep2 = microscopeFStep;

	if (microscopeFStep2.empty()) {
		//Need to know this value
		DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:focalstep");
		return;
	}

	// Check for focalstep voltagemode
	DigitalMicrograph::String microscopeFStepVol("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:voltage focalstep (nm)", microscopeFStepVol);
	std::string microscopeFStep3 = microscopeFStepVol;

	if (VoltageModeEnable) {
		if (microscopeFStep3.empty()) {
			//Need to know this value
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:voltage focalstep");
			return;
		}
	}

	float nmperstep = 1.0f;

	if (VoltageModeEnable == 0) {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	} else {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep3);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	}

	// Setup Variables for Microscope (ARM200)
	// Should add nmperdac entry to be applicable to any microscope
	// Maybe read from tag's and create tag if none already...
	Gatan::uint32 xpixels = 0; //Get from CCD anyway
	Gatan::uint32 ypixels = 0; //Get from CCD anyway

	// Get current DM camera
	Gatan::Camera::Camera camera;

	try {
		camera = Gatan::Camera::GetCurrentCamera();
		Gatan::Camera::CCD_GetSize(camera, &xpixels, &ypixels);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("No Camera Detected", "", "");
		return;
	}

	bool inserted = false;

	try {
		inserted = Gatan::Camera::GetCameraInserted(camera);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Couldn't check camera", "status", "");
		return;
	}

	if (inserted != true) {
		DigitalMicrograph::OpenAndSetProgressWindow("Camera not inserted", "", "");
		return;
	}

	// Want gain normalized imaging unless doing post processing yourself
	Gatan::Camera::AcquisitionProcessing processing = Gatan::Camera::kGainNormalized;
	Gatan::Camera::AcquisitionParameters acqparams;

	try {
		acqparams = Gatan::Camera::CreateAcquisitionParameters(camera, processing, floatExpTime, binning, binning, 0, 0, ypixels,
		            xpixels);
		Gatan::CM::SetDoContinuousReadout(acqparams, true);
		Gatan::CM::SetQualityLevel(acqparams, 1); // Can't remember if fast or slow :D
		Gatan::Camera::Validate_AcquisitionParameters(camera, acqparams);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Problem with acquisition", "parameters", "");
		return;
	}


	// NEW BIT FOR ALTERNATE ACQUISITION
	Gatan::CM::AcquisitionPtr acq = CreateAcquisition( camera, acqparams );
	// Turn into script object for my dms function
	DigitalMicrograph::ScriptObject acqtok =
	    DigitalMicrograph::ScriptObjectProxy<Gatan::Camera::AcquisitionImp, DigitalMicrograph::DMObject>::to_object_token(acq.get());
	Gatan::CM::FrameSetInfoPtr fsi = DMScript::GetFrameSetInfoPtr(acqtok);
	Gatan::Camera::AcquisitionImageSourcePtr acqsource = Gatan::Camera::AcquisitionImageSource::New(acq, fsi, 0);

	// To add calibration information etc..
	Gatan::Camera::AcquisitionImp * Acquisition = acq.get();


	// Check / Create rest of directory structure

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring);

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string());

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string());

	// Make new folder for this focal series
	// First build directory name
	std::string magnification = FormattedMag(mag);
	std::string stdfocalstep = LPCTSTR(focalstep);

	std::string dirname = magnification;

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string() + "\\" + dirname))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string() + "\\" + dirname);

	// Count number of files already in directory
	DigitalMicrograph::TagGroup FilesInDir = DigitalMicrograph::GetFilesInDirectory(filepath.get_string() + "\\" + datestring + "\\" +
	        SpecimenName.get_string() + "\\Area " + AreaNumber.get_string() + "\\" + dirname, 2 );
	int Files = DigitalMicrograph::TagGroupCountTags(FilesInDir);
	Files++;

	std::string filedir = filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
	                      AreaNumber.get_string() + "\\" + dirname + "\\" + t_to_string<int>(Files);

	if (!DigitalMicrograph::DoesDirectoryExist(filedir)) {
		DigitalMicrograph::CreateDirectory(filedir);
	}

	// Find image size and setup new images.
	int imwidth = xpixels / binning;
	int imheight = ypixels / binning;

	// Test if its faster without displaying anything - not really
	// Images wont always be realimage - GIF is signed int..
	DigitalMicrograph::Image view = DigitalMicrograph::RealImage("Through Focus View", 4, imwidth, imheight);
	DigitalMicrograph::Image stack = DigitalMicrograph::RealImage("Defocus Stack", 4, imwidth, imheight, intNumImages);
	DigitalMicrograph::ImageDocument viewdoc = view.GetOrCreateImageDocument();
	DigitalMicrograph::ImageDocument stackdoc = stack.GetOrCreateImageDocument();
	Gatan::PlugIn::ImageDataLocker stackLocker(stack);
	float * stackdata = (float *) stackLocker.get();

	// Maybe this is wrong size, factor binning in here.
	viewdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));

	// Hold image information in vectors so can do all saving at the end (faster).
	std::vector<DigitalMicrograph::Image> Images;
	std::vector<std::string> savelocations;

	// Get HighTension info if using voltage mode.
	double startHT;

	if (VoltageModeEnable == 1 ) {
		startHT = DMScript::EMGetHighTensionOffset();
	}

	// Start Acquisition
	acqsource->BeginAcquisition();

	for (int i = 1; i <= intNumImages; i++) {
		// Establish current focus including rounding

		int pos = (-intNumUnderfocus + (i - 1)) * boost::math::round<int>(floatFocalStep / nmperstep);
		float df = pos * nmperstep;

		DigitalMicrograph::Image AcquiredImage;

		try {
			AcquiredImage = Gatan::Camera::CreateImageForAcquire(acq, "Acquired Image");
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::Result("Couldn't Create Image for Acquire\n");
			return;
		}

		// Try to add acquisition information
		try {
			Acquisition->Acquire_WriteInfo(AcquiredImage, 1);
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			std::string desc = DigitalMicrograph::GetExceptionDescription();
			DigitalMicrograph::Result(desc + "\n");
			DigitalMicrograph::Result("Couldn't Add Acquisition Tags\n");
			return;
		}

		// Try to acquire an image
		DigitalMicrograph::EMPrepareShift();

		// If low dose is enabled then we should move beam away now (ideally it would already be away from where we are now (but then we wouldnt find anything)
		// Then call the other thread to move the beam back in for the specified time intervals...
		// Remember to put beam back afterwards or people will be annoyed :)
		if (LowDoseEnable == 1) {
			DigitalMicrograph::Result("Moving Beam -" + DACShift + "\n");
			DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
		}

		bool acqprmchanged = false;

		try {
			// Not sure of the effect of max waiting time

			// Returns a bool if a read finishes within call?
			// true means restart if it was already finishing a read, this is incase parameters hadnt changed from last time yet?
			if (!acqsource->AcquireTo(AcquiredImage, true, 0.5f, acqprmchanged)) {
				// Now wait for it to finish again but dont restart if it finishes durign call....
				while (!acqsource->AcquireTo(AcquiredImage, false, 0.5f, acqprmchanged)) {
					// Waiting for read to finish
				}

				// Has read finished now??
				if (LowDoseEnable == 1) {
					// MOVE BEAM GET OUT THE WAY
					DigitalMicrograph::Result("Moving Beam +" + DACShift + "\n");
					DigitalMicrograph::EMBeamShift(intDACShift, intDACShift);
					// GET OUT THE WAY BEAM
					// GET OUT THE WAY
				}
			}
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::OpenAndSetProgressWindow("Couldn't Acquire Image", "", "");
			return;
		}

		Gatan::PlugIn::ImageDataLocker acquireLocker(AcquiredImage);
		Gatan::PlugIn::ImageDataLocker throughfocusLocker(view);

		// Check if camera records float or integer type images.
		if (DigitalMicrograph::IsFloatImage(AcquiredImage.get())) {
			float * acquireddata = (float *) acquireLocker.get();
			float * viewdata = (float *) throughfocusLocker.get();

			for (int j = 0; j < imwidth * imheight; j++) {
				//viewdata[j] = acquireddata[j];
				float val = *(acquireddata++);
				*(viewdata++) = val;
				*(stackdata++) = val;
			}
		}

		else
			if (DigitalMicrograph::IsIntegerDataType(AcquiredImage.get(), 4, true)) {
				int * acquireddata = (int *) acquireLocker.get();
				float * viewdata = (float *) throughfocusLocker.get();

				for (int j = 0; j < imwidth * imheight; j++) {
					//viewdata[j] = acquireddata[j];
					int val = *(acquireddata++);
					*(viewdata++) = val;
					*(stackdata++) = val;
				}
			}

			// If its not float or integer this doesn't work.
			else {
				DigitalMicrograph::OpenAndSetProgressWindow("Image not of expected type", "", "");
				return;
			}

		acquireLocker.~ImageDataLocker();
		throughfocusLocker.~ImageDataLocker();

		view.DataChanged();
		view.SetName(t_to_string(i) + " - " + t_to_string(df) + "nm"); // Rename and Save Image.

		std::string saveloc = filedir + "// " + t_to_string(i) + " " + t_to_string(df) + "nm";
		AcquiredImage.SetName(t_to_string(i) + " " + t_to_string(df) + "nm");

		AddAberrationTags(AcquiredImage);

		// Add image and savelocation to vector for processing later
		Images.push_back(AcquiredImage);
		savelocations.push_back(saveloc);

		// Don't save here anymore, delay til after all acquisitions
		//DigitalMicrograph::SaveAsGatan(AcquiredImage,saveloc);

		// Don't show image if using continuous updating view
		//DigitalMicrograph::ShowImage(AcquiredImage);
		//AcquiredImage.GetOrCreateImageDocument().Clean();

		// TODO: Add to Stack Image (For Later)

		// Update Progress
		std::string defoc = t_to_string(df);
		std::string st3 = "@ " + defoc + "nm";
		m_progressBar.SetPos(i);
		DigitalMicrograph::OpenAndSetProgressWindow("Image", "Acquired", st3.c_str());
	}

	// Probably important...
	acqsource->FinishAcquisition();

	// Dont bother saving if we have aborted..
	if (GetStatus() != Stopping) {
		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Now Saving Images...");

		for (int i = 0; i < intNumImages; i++) {
			DigitalMicrograph::SaveAsGatan(Images[i], savelocations[i]);
		}

		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Images Saved.");
	}


	// Get rid of through focus view
	view.GetOrCreateImageDocument().Clean();
	DigitalMicrograph::CloseImage(view);

	// Close stack for editing
	stackLocker.~ImageDataLocker();

	float xscale = Images[0].GetDimensionScale(0);
	float yscale = Images[0].GetDimensionScale(1);

	DigitalMicrograph::TagGroup stacktags = stack.GetTagGroup();
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Focal Series:Adjusted focalstep", t_to_string(focalstep));
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Microscope Info:Voltage", t_to_string(operatingVoltage));

	stack.SetDimensionScale(0, xscale);
	stack.SetDimensionScale(1, yscale);

	stack.DataChanged();

	// Show stack view
	stackdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));
	// Clear vectors and other memory just to be sure

	Images.clear();
	savelocations.clear();

	// If it was low dose we need to return beam to normal position afterwards
	if (LowDoseEnable == 1) {
		DigitalMicrograph::Result("Moving Beam -" + DACShift + "\n");
		DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
	}
}

void Acquisition::SeriesLowDoseAcquire() {
	// Retrieve Values from Dialog ready to start imaging.
	UpdateData(TRUE);

	// Cancel any currently updating live views... didn't seem to work.
	//DigitalMicrograph::StopAllAcquisitionDaemons();
	Gatan::Camera::StopCurrentCameraViewer(false); // Stop it but dont close it...

	// Get size of shift to be put into low dose mode.
	intDACShift = boost::lexical_cast<int>(DACShift);

	// First check we have a save location specified.
	DigitalMicrograph::String filepath("");
	DigitalMicrograph::TagGroup pers = DigitalMicrograph::GetPersistentTagGroup();
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:filepath", filepath);

	// Check path exists
	if (!DigitalMicrograph::DoesDirectoryExist( filepath )) {
		DigitalMicrograph::OpenAndSetProgressWindow("Root Directory Not Available", "Or Not Specified", "             :(");
		return;
	}

	// Check for Specimen Name
	DigitalMicrograph::String SpecimenName("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Specimen", SpecimenName);
	// Doesn't really matter if blank.

	// Get Area Number
	DigitalMicrograph::String AreaNumber("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Operator Panel:Area", AreaNumber);

	// Find out todays date
	CTime theTime;
	theTime = CTime::GetCurrentTime();
	CString date = theTime.Format( "%d %B %Y" );
	std::string datestring = ((LPCTSTR)date);



	// Set progress bar to correct size
	intNumImages = boost::lexical_cast<int>(numImages);
	m_progressBar.SetRange(0, intNumImages);

	intNumUnderfocus = boost::lexical_cast<int>(numUnderfocus);
	floatFocalStep = boost::lexical_cast<float>(focalstep);
	floatExpTime = boost::lexical_cast<float>(expTime);

	int binning(0);

	try {
		binning = boost::lexical_cast<int>(binningText);
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	if (!(binning == 1 || binning == 2 || binning == 4)) {
		DigitalMicrograph::OpenAndSetProgressWindow("Please Select", "Appropriate Binning", "1,2 or 4");
		return;
	}

	// Check values make sense i.e underfocus < total images
	if ((!(intNumUnderfocus <= intNumImages)) || intNumImages < 1) {
		DigitalMicrograph::OpenAndSetProgressWindow("Number of Images", "Incorrectly Set", "");
		return;
	}

	// Get current microscope settings
	DigitalMicrograph::String operatingVoltage;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Voltage", operatingVoltage);
	DigitalMicrograph::String indicatedMag;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Indicated Magnification", indicatedMag);

	float mag = 0.0f;

	try {
		mag = boost::lexical_cast<float>(indicatedMag.get_string());
	} catch ( boost::bad_lexical_cast const &) {
		DigitalMicrograph::OpenAndSetProgressWindow("Indicated Magnification", "taggroup", "incorrectly set");
		return;
	}

	DigitalMicrograph::String operationMode;
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Microscope Info:Operation Mode", operationMode);

	// Check for focalstep
	DigitalMicrograph::String microscopeFStep("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:focalstep (nm)", microscopeFStep);
	std::string microscopeFStep2 = microscopeFStep;

	if (microscopeFStep2.empty()) {
		//Need to know this value
		DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:focalstep");
		return;
	}

	// Check for focalstep voltagemode
	DigitalMicrograph::String microscopeFStepVol("");
	DigitalMicrograph::TagGroupGetTagAsString(pers, "Focal Series:voltage focalstep (nm)", microscopeFStepVol);
	std::string microscopeFStep3 = microscopeFStepVol;

	if (VoltageModeEnable) {
		if (microscopeFStep3.empty()) {
			//Need to know this value
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep not specified", "in taggroup", "Focal Series:voltage focalstep");
			return;
		}
	}

	float nmperstep = 1.0f;

	if (VoltageModeEnable == 0) {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	} else {
		try {
			nmperstep = boost::lexical_cast<float>(microscopeFStep3);
		} catch ( boost::bad_lexical_cast const &) {
			DigitalMicrograph::OpenAndSetProgressWindow("focalstep indicated", "cannot be parsed", "to a number");
			return;
		}
	}

	// Setup Variables for Microscope (ARM200)
	// Should add nmperdac entry to be applicable to any microscope
	// Maybe read from tag's and create tag if none already...
	Gatan::uint32 xpixels = 0; //Get from CCD anyway
	Gatan::uint32 ypixels = 0; //Get from CCD anyway

	// Get current DM camera
	Gatan::Camera::Camera camera;

	try {
		camera = Gatan::Camera::GetCurrentCamera();
		Gatan::Camera::CCD_GetSize(camera, &xpixels, &ypixels);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("No Camera Detected", "", "");
		return;
	}

	bool inserted = false;

	try {
		inserted = Gatan::Camera::GetCameraInserted(camera);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Couldn't check camera", "status", "");
		return;
	}

	if (inserted != true) {
		DigitalMicrograph::OpenAndSetProgressWindow("Camera not inserted", "", "");
		return;
	}

	// Want gain normalized imaging unless doing post processing yourself
	Gatan::Camera::AcquisitionProcessing processing = Gatan::Camera::kGainNormalized;
	Gatan::Camera::AcquisitionParameters acqparams;

	try {
		acqparams = Gatan::Camera::CreateAcquisitionParameters(camera, processing, floatExpTime, binning, binning, 0, 0, ypixels,
		            xpixels);
		Gatan::CM::SetDoContinuousReadout(acqparams, true);
		Gatan::CM::SetQualityLevel(acqparams, 1); // Can't remember if fast or slow :D
		Gatan::Camera::Validate_AcquisitionParameters(camera, acqparams);
	} catch (...) {
		short error;
		long context;
		DigitalMicrograph::GetException(&error, &context);
		DigitalMicrograph::ErrorDialog(error);
		DigitalMicrograph::OpenAndSetProgressWindow("Problem with acquisition", "parameters", "");
		return;
	}


	// NEW BIT FOR ALTERNATE ACQUISITION
	Gatan::CM::AcquisitionPtr acq = CreateAcquisition( camera, acqparams );
	// Turn into script object for my dms function
	DigitalMicrograph::ScriptObject acqtok =
	    DigitalMicrograph::ScriptObjectProxy<Gatan::Camera::AcquisitionImp, DigitalMicrograph::DMObject>::to_object_token(acq.get());
	Gatan::CM::FrameSetInfoPtr fsi = DMScript::GetFrameSetInfoPtr(acqtok);
	Gatan::Camera::AcquisitionImageSourcePtr acqsource = Gatan::Camera::AcquisitionImageSource::New(acq, fsi, 0);

	// To add calibration information etc..
	Gatan::Camera::AcquisitionImp * Acquisition = acq.get();


	// Check / Create rest of directory structure

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring);

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string());

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string()))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string());

	// Make new folder for this focal series
	// First build directory name
	std::string magnification = FormattedMag(mag);
	std::string stdfocalstep = LPCTSTR(focalstep);

	std::string dirname = "focalseries " + stdfocalstep + "nm " + magnification;

	if (!DigitalMicrograph::DoesDirectoryExist(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() +
	        "\\Area " + AreaNumber.get_string() + "\\" + dirname))
		DigitalMicrograph::CreateDirectory(filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
		                                   AreaNumber.get_string() + "\\" + dirname);

	// Count number of files already in directory
	DigitalMicrograph::TagGroup FilesInDir = DigitalMicrograph::GetFilesInDirectory(filepath.get_string() + "\\" + datestring + "\\" +
	        SpecimenName.get_string() + "\\Area " + AreaNumber.get_string() + "\\" + dirname, 2 );
	int Files = DigitalMicrograph::TagGroupCountTags(FilesInDir);
	Files++;

	std::string filedir = filepath.get_string() + "\\" + datestring + "\\" + SpecimenName.get_string() + "\\Area " +
	                      AreaNumber.get_string() + "\\" + dirname + "\\" + t_to_string<int>(Files);

	if (!DigitalMicrograph::DoesDirectoryExist(filedir)) {
		DigitalMicrograph::CreateDirectory(filedir);
	}

	// Find image size and setup new images.
	int imwidth = xpixels / binning;
	int imheight = ypixels / binning;

	// Test if its faster without displaying anything - not really
	// Images wont always be realimage - GIF is signed int..
	DigitalMicrograph::Image view = DigitalMicrograph::RealImage("Through Focus View", 4, imwidth, imheight);
	DigitalMicrograph::Image stack = DigitalMicrograph::RealImage("Defocus Stack", 4, imwidth, imheight, intNumImages);
	DigitalMicrograph::ImageDocument viewdoc = view.GetOrCreateImageDocument();
	DigitalMicrograph::ImageDocument stackdoc = stack.GetOrCreateImageDocument();
	Gatan::PlugIn::ImageDataLocker stackLocker(stack);
	float * stackdata = (float *) stackLocker.get();

	// Maybe this is wrong size, factor binning in here.
	viewdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));

	// Hold image information in vectors so can do all saving at the end (faster).
	std::vector<DigitalMicrograph::Image> Images;
	std::vector<std::string> savelocations;

	// Get HighTension info if using voltage mode.
	double startHT;

	if (VoltageModeEnable == 1 ) {
		startHT = DMScript::EMGetHighTensionOffset();
	}


	if (VoltageModeEnable == 0 ) {
		// Go to most underfocus
		if (focalstep != "0") {
			try {
				DigitalMicrograph::EMChangeFocus(-intNumUnderfocus * boost::math::round<int>(floatFocalStep / nmperstep));
			} catch (...) {
				short error;
				long context;
				DigitalMicrograph::GetException(&error, &context);
				DigitalMicrograph::ErrorDialog(error);
				DigitalMicrograph::OpenAndSetProgressWindow("Problem with changing", "focus", "");
				return;
			}
		}
	} else {
		// Go to most underfocus
		if (focalstep != "0") {
			try {
				DMScript::EMSetHighTensionOffset(startHT - intNumUnderfocus * boost::math::round<int>(floatFocalStep / nmperstep));
			} catch (...) {
				short error;
				long context;
				DigitalMicrograph::GetException(&error, &context);
				DigitalMicrograph::ErrorDialog(error);
				DigitalMicrograph::OpenAndSetProgressWindow("Problem with changing", "focus", "");
				return;
			}
		}
	}


	// Start Acquisition
	acqsource->BeginAcquisition();

	for (int i = 1; i <= intNumImages; i++) {
		// Check if we are requested to stop
		if (GetStatus() == Stopping) {
			if (VoltageModeEnable == 0 ) {
				// Fix by how much we have changed focus so far.
				DigitalMicrograph::EMChangeFocus((intNumUnderfocus - (i - 1))*boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT);
			}

			// Beam position will be restored at end anyway

			// Print useful message
			DigitalMicrograph::Result("Acquisition has been aborted \n");

			// Exit loop
			break;
		}

		// Establish current focus including rounding

		int pos = (-intNumUnderfocus + (i - 1)) * boost::math::round<int>(floatFocalStep / nmperstep);
		float df = pos * nmperstep;

		DigitalMicrograph::Image AcquiredImage;

		try {
			AcquiredImage = Gatan::Camera::CreateImageForAcquire(acq, "Acquired Image");
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::Result("Couldn't Create Image for Acquire\n");
			return;
		}

		// Try to add acquisition information
		try {
			Acquisition->Acquire_WriteInfo(AcquiredImage, 1);
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			std::string desc = DigitalMicrograph::GetExceptionDescription();
			DigitalMicrograph::Result(desc + "\n");
			DigitalMicrograph::Result("Couldn't Add Acquisition Tags\n");
			return;
		}

		// Try to acquire an image
		DigitalMicrograph::EMPrepareShift();

		// If low dose is enabled then we should move beam away now (ideally it would already be away from where we are now (but then we wouldnt find anything)
		// Then call the other thread to move the beam back in for the specified time intervals...
		// Remember to put beam back afterwards or people will be annoyed :)
		if (LowDoseEnable == 1) {
			DigitalMicrograph::Result("Moving Beam -" + DACShift + "\n");
			DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
		}

		bool acqprmchanged = false;

		try {
			// Not sure of the effect of max waiting time

			// Returns a bool if a read finishes within call?
			// true means restart if it was already finishing a read, this is incase parameters hadnt changed from last time yet?
			if (!acqsource->AcquireTo(AcquiredImage, true, 0.5f, acqprmchanged)) {
				// Now wait for it to finish again but dont restart if it finishes durign call....
				while (!acqsource->AcquireTo(AcquiredImage, false, 0.5f, acqprmchanged)) {
					// Waiting for read to finish
				}

				// Has read finished now??
				if (LowDoseEnable == 1) {
					// MOVE BEAM GET OUT THE WAY
					DigitalMicrograph::Result("Moving Beam +" + DACShift + "\n");
					DigitalMicrograph::EMBeamShift(intDACShift, intDACShift);
					// GET OUT THE WAY BEAM
					// GET OUT THE WAY
				}
			}
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::OpenAndSetProgressWindow("Couldn't Acquire Image", "", "");
			return;
		}

		Gatan::PlugIn::ImageDataLocker acquireLocker(AcquiredImage);
		Gatan::PlugIn::ImageDataLocker throughfocusLocker(view);

		// Check if camera records float or integer type images.
		if (DigitalMicrograph::IsFloatImage(AcquiredImage.get())) {
			float * acquireddata = (float *) acquireLocker.get();
			float * viewdata = (float *) throughfocusLocker.get();

			for (int j = 0; j < imwidth * imheight; j++) {
				//viewdata[j] = acquireddata[j];
				float val = *(acquireddata++);
				*(viewdata++) = val;
				*(stackdata++) = val;
			}
		}

		else
			if (DigitalMicrograph::IsIntegerDataType(AcquiredImage.get(), 4, true)) {
				int * acquireddata = (int *) acquireLocker.get();
				float * viewdata = (float *) throughfocusLocker.get();

				for (int j = 0; j < imwidth * imheight; j++) {
					//viewdata[j] = acquireddata[j];
					int val = *(acquireddata++);
					*(viewdata++) = val;
					*(stackdata++) = val;
				}
			}

			// If its not float or integer this doesn't work.
			else {
				DigitalMicrograph::OpenAndSetProgressWindow("Image not of expected type", "", "");
				return;
			}

		acquireLocker.~ImageDataLocker();
		throughfocusLocker.~ImageDataLocker();

		view.DataChanged();
		view.SetName(t_to_string(i) + " - " + t_to_string(df) + "nm"); // Rename and Save Image.

		std::string saveloc = filedir + "// " + t_to_string(i) + " " + t_to_string(df) + "nm";
		AcquiredImage.SetName(t_to_string(i) + " " + t_to_string(df) + "nm");

		try {
			AddFocalTags(AcquiredImage, i, df, boost::math::round<int>(floatFocalStep / nmperstep)*nmperstep, intNumImages);
		} catch (...) {
			short error;
			long context;
			DigitalMicrograph::GetException(&error, &context);
			DigitalMicrograph::ErrorDialog(error);
			DigitalMicrograph::Result("Can't Add Focal Series Tags \n");
		}


		AddAberrationTags(AcquiredImage);

		// Add image and savelocation to vector for processing later
		Images.push_back(AcquiredImage);
		savelocations.push_back(saveloc);

		// Don't save here anymore, delay til after all acquisitions
		//DigitalMicrograph::SaveAsGatan(AcquiredImage,saveloc);

		// Don't show image if using continuous updating view
		//DigitalMicrograph::ShowImage(AcquiredImage);
		//AcquiredImage.GetOrCreateImageDocument().Clean();

		// TODO: Add to Stack Image (For Later)

		// Update Progress
		std::string defoc = t_to_string(df);
		std::string st3 = "@ " + defoc + "nm";
		m_progressBar.SetPos(i);
		DigitalMicrograph::OpenAndSetProgressWindow("Image", "Acquired", st3.c_str());

		if (focalstep != "0") {
			if (VoltageModeEnable == 0 ) {
				// Shift Focus
				DigitalMicrograph::EMChangeFocus(boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT - (intNumUnderfocus + i)*boost::math::round<int>(floatFocalStep / nmperstep));
			}
		}
	}

	// If we have been requested to stop then focus is already restored earlier
	if (GetStatus() != Stopping) {
		if (focalstep != "0") {
			if (VoltageModeEnable == 0 ) {
				// Reset Focus (one more than might be expected because it changes again after last image)
				DigitalMicrograph::EMChangeFocus(-(intNumImages - intNumUnderfocus)*boost::math::round<int>(floatFocalStep / nmperstep));
			} else {
				DMScript::EMSetHighTensionOffset(startHT);
			}
		}
	}

	// Probably important...
	acqsource->FinishAcquisition();

	// Dont bother saving if we have aborted..
	if (GetStatus() != Stopping) {
		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Now Saving Images...");

		for (int i = 0; i < intNumImages; i++) {
			DigitalMicrograph::SaveAsGatan(Images[i], savelocations[i]);
		}

		DigitalMicrograph::OpenAndSetProgressWindow("Focal Series Acquisition", "Complete", "Images Saved.");
	}


	// Get rid of through focus view
	view.GetOrCreateImageDocument().Clean();
	DigitalMicrograph::CloseImage(view);

	// Close stack for editing
	stackLocker.~ImageDataLocker();

	float xscale = Images[0].GetDimensionScale(0);
	float yscale = Images[0].GetDimensionScale(1);

	DigitalMicrograph::TagGroup stacktags = stack.GetTagGroup();
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Focal Series:Adjusted focalstep", t_to_string(focalstep));
	DigitalMicrograph::TagGroupSetTagAsString(stacktags, "Microscope Info:Voltage", t_to_string(operatingVoltage));

	stack.SetDimensionScale(0, xscale);
	stack.SetDimensionScale(1, yscale);

	stack.DataChanged();

	// Show stack view
	stackdoc.ShowAtRect(30, 30, 30 + ypixels / (4), 30 + xpixels / (4));
	// Clear vectors and other memory just to be sure

	Images.clear();
	savelocations.clear();

	// If it was low dose we need to return beam to normal position afterwards
	if (LowDoseEnable == 1) {
		DigitalMicrograph::Result("Moving Beam -" + DACShift + "\n");
		DigitalMicrograph::EMBeamShift(-intDACShift, -intDACShift);
	}
}

void Acquisition::PreShiftBeam() {
	// Get size of shift to be put into low dose mode.
	intDACShift = boost::lexical_cast<int>(DACShift);

	// Call this before any beam shifting is to be done.
	DigitalMicrograph::EMPrepareShift();

	// If low dose start by moving beam out of position

	// Technically this bit should be in a seperate function, or stage movement should come afterwards, because we want to take a series somewhere we havent been before.
	if (LowDoseEnable == 1) {
		try {
			DigitalMicrograph::Result("Moving Beam +" + DACShift + "\n");
			DigitalMicrograph::EMBeamShift(intDACShift, intDACShift);
		} catch (...) {
			DigitalMicrograph::Result("Unable to shift beam position \n");
			return;
		}
	}
}