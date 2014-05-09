#pragma once
#include "CL\cl.h"
#include "clKernel.h"

class clState
{
private:
	clState(void);
	~clState(void);
public:
	static cl_int status;
	static cl_context context;
	static clDevice* cldev;
	static clQueue* clq;
	static bool OpenCLAvailable;

	static void Setup();
	static cl_int GetStatus();
	static cl_context GetContext();
	static clDevice* GetDevicePtr();
	static clQueue* GetQueuePtr();
};

