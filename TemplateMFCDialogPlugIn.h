// TemplateMFCDialogPlugIn.h : main header file for the TemplateMFCDialogPlugIn DLL
//

#pragma once


#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

#include "DMDialog.h"
#include "SeriesFunctions.h"

// CTemplateMFCDialogPlugInApp
// See TemplateMFCDialogPlugIn.cpp for the implementation of this class
//
#pragma unmanaged
class CTemplateMFCDialogPlugInApp : public Gatan::PlugIn::PlugInMain, public Gatan::PlugIn::MFCPlugInUtility
{
private:
	ulong_ptr_t fMFCPaletteHandle1;
	ulong_ptr_t fMFCPaletteHandle2;
public:
	CTemplateMFCDialogPlugInApp()  {}
	~CTemplateMFCDialogPlugInApp() {}

	virtual void Start(void);
	virtual void End(void);
};
