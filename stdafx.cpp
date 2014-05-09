// stdafx.cpp : source file that includes just the standard includes
// TemplateMFCDialogPlugIn.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"


void Debug(std::string in)
{
	std::string out = "DEBUG: "+in+" \n";
	DigitalMicrograph::Result(out);
}