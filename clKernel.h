#pragma once
#include "CL/OpenCl.h"
#include "DMPlugInStubs.h"

class clState;

enum clTypes
{
	clFloat = 0,
	clFloat2 = 1,
	clDouble = 2,
	clDouble2 = 3,
	clUInt = 4,
	clInt = 5
}; 

class clDevice
{
public:
	cl_uint numDevices;
	cl_device_id devices;

	clDevice::clDevice(cl_uint numDevices, cl_device_id* devices);
	clDevice::clDevice(cl_uint numDevices, cl_device_id device);
	clDevice();
	~clDevice();
};

class clQueue
{
public:
	cl_command_queue cmdQueue;
	cl_int status;

	//clQueue(cl_command_queue &cmdQueue);
	cl_int SetupQueue(cl_context &context,cl_device_id device);
	clQueue();
	~clQueue();

};


class clKernel
{
public:
	clKernel(void);
	~clKernel(void);

	const char* kernelcode;
	cl_program kernelprogram;
	//cl_int status;
	//cl_context* context;
	cl_kernel kernel;
	std::string kernelname;
	size_t log;
	//clDevice* cldev;
	//clQueue* clq;

	// Old constructor before clState gave global context... (remove when get rid of all references to it.
	clKernel(const char* codestring, cl_context &context, clDevice* cldev, std::string kernelname,clQueue* Queue);

	void SetCodeAndName(const char* codestring, std::string kernelname);
	void BuildKernel();
	cl_int StatusOK();
	void Enqueue(size_t* globalWorkSize);
	void Enqueue3D(size_t* globalWorkSize);
	void Enqueue3D(size_t* globalWorkSize, size_t* localWorkSize);

	// Dont really need these with badass template :D
	//void SetArg(int position,cl_mem argument);
	//void SetArg(int position,int argument);
	//void SetArg(int position,float argument);

	// Function definition has to be in header for templates...
	// Sets arguments for clKernel
	template <class T> void SetArgT(int position, T &arg) 
	{
		clSetKernelArg(kernel,position,sizeof(T),&arg);
	}

	void SetArgLocalMemory(int position, int size, clTypes type) 
	{
		switch (type)
		{
		case clFloat:
			clSetKernelArg(kernel,position,size*sizeof(cl_float),NULL);
			break;
		case clFloat2:
			clSetKernelArg(kernel,position,size*sizeof(cl_float2),NULL);
			break;
		case clDouble:
			clSetKernelArg(kernel,position,size*sizeof(cl_double),NULL);
			break;
		case clDouble2:
			clSetKernelArg(kernel,position,size*sizeof(cl_double2),NULL);
			break;
		case clUInt:
			clSetKernelArg(kernel,position,size*sizeof(cl_uint),NULL);
			break;
		case clInt:
			clSetKernelArg(kernel,position,size*sizeof(cl_int),NULL);
			break;
		}
		
	}


};

