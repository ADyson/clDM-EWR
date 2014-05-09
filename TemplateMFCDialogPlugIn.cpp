// TemplateMFCDialogPlugIn.cpp : Defines the initialization routines for the DLL.
//



#include "stdafx.h"
#include "TemplateMFCDialogPlugIn.h"
#include "MyTabCtrl.h"
#include "MutualInformation.h"

#define _GATAN_USE_STL_STRING

#include "DMPluginMain.h"
extern AFX_EXTENSION_MODULE gPlugInExtensionModule = { NULL, NULL };

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

void CTemplateMFCDialogPlugInApp::Start(void)
{
	AddFunction("float MutualInformationMap(RealImagePtr,RealImagePtr,long,long)", (void *)MutualInformationMap);
	AddFunction("float MutualInformationCL(RealImagePtr,RealImagePtr,long,long)", (void *)MutualInformationCL);
	AddFunction("void XCFPCFPCPCF(RealImagePtr,RealImagePtr,long)", (void *)XCFPCFPCPCF);
	AddFunction("void Magnification(RealImagePtr,RealImagePtr)", (void *)Magnification);
	AddFunction("void MagnificationTest(RealImagePtr,RealImagePtr,float,long,float,float)", (void *)MagnificationTest);
	fMFCPaletteHandle1=RegisterMFCPalette("OCL - EWR", GetPlugInOSHandle(), RUNTIME_CLASS(CDMDialog),3);
}
void CTemplateMFCDialogPlugInApp::End(void)
{
	UnregisterMFCPalette(fMFCPaletteHandle1);
}

// The one and only CTemplateMFCDialogPlugInApp object

CTemplateMFCDialogPlugInApp theApp;

