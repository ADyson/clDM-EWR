#include "stdafx.h"
#include "standardfunctions.h"

// define some useful functions
int GetWindowString(HWND hwnd, std::string &s) 
{
	char buffer[65536];
 
	int txtlen=GetWindowTextLength(hwnd);
	GetWindowText(hwnd, buffer, txtlen);
 
	s = buffer;
	return txtlen;
}

int isPowerOfTwo (unsigned int x)
{
	return ((x != 0) && !(x & (x - 1)));
}

void SetWindowString(HWND hwnd, std::string s) 
{    
	SetWindowText(hwnd, s.c_str());
}



