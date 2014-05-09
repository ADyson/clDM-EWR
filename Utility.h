
#pragma once
#include "boost/lexical_cast.hpp"
#include "CL/OpenCl.h"
#include "clKernel.h"
#include "stdafx.h"

class Utility
{
public:
	Utility(void);
	~Utility(void);

	template <class T> static bool GetDialogString(T &arg, CString variable) 
	{
		try
		{
			arg = boost::lexical_cast<T>(variable);
		}
		catch(...)
		{
			DigitalMicrograph::Result("Couldn't parse " + variable);
			return false;
		}
		return true;
	}

	static void SetProgressWindow(std::string s1, std::string s2, std::string s3)
	{
		DigitalMicrograph::OpenAndSetProgressWindow(s1.c_str(),s2.c_str(),s3.c_str());
	}

	static void SetResultWindow(std::string s1)
	{
		std::string s2 = s1 + "\n";
		DigitalMicrograph::Result(s2.c_str());
	}

	static DigitalMicrograph::Image PrintCLMemToImage(cl_mem buffer,std::string name, int width, int height, clTypes type, clQueue* clq);

	static DigitalMicrograph::Image PrintCLMemToImagePlusOne(cl_mem buffer,std::string name, int width, int height, clTypes type, clQueue* clq);

	static int LinearAccess3D(int i , int starti, int j, int startj, int k, int sizei, int sizej);

	static int LinearAccess3D(int i , int j, int k, int sizei, int sizej);

	static DigitalMicrograph::Image DisplayArray(float* buffer,std::string name, int width, int height);

	static DigitalMicrograph::Image DisplayArray(unsigned int* buffer,std::string name, int width, int height);

	static DigitalMicrograph::Image DisplayArray(int* buffer,std::string name, int width, int height);
};

