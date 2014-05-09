#include "StdAfx.h"
#include "Utility.h"


Utility::Utility(void)
{
}


Utility::~Utility(void)
{
}

DigitalMicrograph::Image Utility::PrintCLMemToImage(cl_mem buffer,std::string name, int width, int height, clTypes type, clQueue* clq)
{
	DM::Image Test;

	CString cname = name.c_str();

	if(type == clFloat)
	{
		Test = DM::RealImage(cname,4,width,height);
	}
	else if (type == clFloat2)
	{
		Test = DM::ComplexImage(cname,8,width,height);
	}

	Gatan::PlugIn::ImageDataLocker locker(Test);

	if(type == clFloat)
	{
		std::vector<float> host (width*height);
		clEnqueueReadBuffer(clq->cmdQueue,buffer,CL_TRUE,0,width*height*sizeof(cl_float),&host[0],0,NULL,NULL);

		float* data =(float*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i] = host[i];
		}

		host.clear();
		locker.~ImageDataLocker();

		Test.GetOrCreateImageDocument().Show();
	}
	else if (type == clFloat2)
	{
		std::vector<cl_float2> host (width*height);
		clEnqueueReadBuffer(clq->cmdQueue,buffer,CL_TRUE,0,width*height*sizeof(cl_float2),&host[0],0,NULL,NULL);

		cl_float2* data =(cl_float2*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i].s[0] = host[i].s[0];
			data[i].s[1] = host[i].s[1];
		}

		host.clear();
		locker.~ImageDataLocker();
		Test.GetOrCreateImageDocument().Show();
	}

	return Test;

}

DigitalMicrograph::Image Utility::PrintCLMemToImagePlusOne(cl_mem buffer,std::string name, int width, int height, clTypes type, clQueue* clq)
{
	DigitalMicrograph::Image Test;

	CString cname = name.c_str();

	if(type == clFloat)
	{
		Test = DM::RealImage(cname,4,width,height);
	}
	else if (type == clFloat2)
	{
		Test = DM::ComplexImage(cname,8,width,height);
	}

	Gatan::PlugIn::ImageDataLocker locker(Test);

	if(type == clFloat)
	{
		std::vector<float> host (width*height);
		clEnqueueReadBuffer(clq->cmdQueue,buffer,CL_TRUE,0,width*height*sizeof(cl_float),&host[0],0,NULL,NULL);

		float* data =(float*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i] = host[i];
		}

		host.clear();
		locker.~ImageDataLocker();

		Test.GetOrCreateImageDocument().Show();
	}
	else if (type == clFloat2)
	{
		std::vector<cl_float2> host (width*height);
		clEnqueueReadBuffer(clq->cmdQueue,buffer,CL_TRUE,0,width*height*sizeof(cl_float2),&host[0],0,NULL,NULL);

		cl_float2* data =(cl_float2*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i].s[0] = host[i].s[0] + 1;
			data[i].s[1] = host[i].s[1];
		}

		host.clear();
		locker.~ImageDataLocker();
		Test.GetOrCreateImageDocument().Show();
		Test.GetImageDisplay(0).SetComplexMode(5);
	}

	return Test;

}

int Utility::LinearAccess3D(int i , int starti, int j, int startj, int k, int sizei, int sizej)
{
	return i+starti + (j+startj)*sizei + k*sizei*sizej;
}

int Utility::LinearAccess3D(int i , int j, int k, int sizei, int sizej)
{
	return i + j*sizei + k*sizei*sizej;
}

DigitalMicrograph::Image Utility::DisplayArray(float* buffer,std::string name, int width, int height)
{
	DM::Image Test;

	CString cname = name.c_str();

	Test = DigitalMicrograph::RealImage(cname,4,width,height);


	Gatan::PlugIn::ImageDataLocker locker(Test);



		float* data =(float*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i] = buffer[i];
		}


		locker.~ImageDataLocker();

		Test.GetOrCreateImageDocument().Show();

	return Test;

}

DigitalMicrograph::Image Utility::DisplayArray(unsigned int* buffer,std::string name, int width, int height)
{
	DM::Image Test;

	CString cname = name.c_str();

	Test = DigitalMicrograph::IntegerImage(cname,4,false,width,height);


	Gatan::PlugIn::ImageDataLocker locker(Test);



		unsigned int* data =(unsigned int*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i] = buffer[i];
		}


		locker.~ImageDataLocker();

		Test.GetOrCreateImageDocument().Show();

	return Test;

}

DigitalMicrograph::Image Utility::DisplayArray(int* buffer,std::string name, int width, int height)
{
	DM::Image Test;

	CString cname = name.c_str();

	Test = DigitalMicrograph::IntegerImage(cname,4,true,width,height);


	Gatan::PlugIn::ImageDataLocker locker(Test);



		unsigned int* data =(unsigned int*) locker.get();

		for(int i = 0 ; i< width*height ; i++)
		{
			data[i] = buffer[i];
		}


		locker.~ImageDataLocker();

		Test.GetOrCreateImageDocument().Show();

	return Test;

}