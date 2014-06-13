#include "StdAfx.h"
#include "clKernel.h"
#include "standardfunctions.h"
#include "clState.h"


clKernel::clKernel(void) {
}


clKernel::~clKernel(void) {
	clReleaseProgram(kernelprogram);
	clReleaseKernel(kernel);
}

// Constructor that sets all Command Queues and Contexts etc..
clKernel::clKernel(const char * codestring, cl_context & context, clDevice * cldev, std::string kernelname, clQueue * Queue) {
	//this->cldev = cldev;
	//this->context = &context;
	//this->numDevices = numdevices;
	//this->devices = devices;
	this->kernelcode = codestring;
	this->kernelname = kernelname;
	//this->clq = Queue;
}

void clKernel::SetCodeAndName(const char * codestring, std::string kernelname) {
	this->kernelcode = codestring;
	this->kernelname = kernelname;
}

void clKernel::BuildKernel() {
	this->kernelprogram = clCreateProgramWithSource(clState::context, 1, &kernelcode, NULL, &clState::status);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " source " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		return;
	}

	clState::status = clBuildProgram(kernelprogram, clState::cldev->numDevices, &clState::cldev->devices, NULL, NULL, NULL);

	clGetProgramBuildInfo(kernelprogram, clState::cldev->devices, CL_PROGRAM_BUILD_LOG, 0, NULL, &log);
	char * buildlog = (char *)malloc(log * sizeof(char));
	clGetProgramBuildInfo(kernelprogram, clState::cldev->devices, CL_PROGRAM_BUILD_LOG, log, buildlog, NULL);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " build " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		DigitalMicrograph::Result(buildlog);
		DigitalMicrograph::Result("\n");
		return;
	}

	free(buildlog);

	this->kernel = clCreateKernel(kernelprogram, kernelname.c_str(), &clState::status);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " kernel " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		return;
	}

}

// 0 is success
cl_int clKernel::StatusOK() {
	if (clState::status != CL_SUCCESS) {
		return clState::status;
	} else return clState::status;
}

void clKernel::Enqueue(size_t * globalWorkSize ) {
	clState::status = clEnqueueNDRangeKernel(clState::clq->cmdQueue, kernel, 2, NULL, globalWorkSize, NULL, 0, NULL, NULL);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " enqueue " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		return;
	}
}

void clKernel::Enqueue3D(size_t * globalWorkSize, size_t * localWorkSize ) {
	clState::status = clEnqueueNDRangeKernel(clState::clq->cmdQueue, kernel, 3, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " enqueue " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		return;
	}
}
void clKernel::Enqueue3D(size_t * globalWorkSize ) {
	clState::status = clEnqueueNDRangeKernel(clState::clq->cmdQueue, kernel, 3, NULL, globalWorkSize, NULL, 0, NULL, NULL);

	if (!clState::status == 0) {
		DigitalMicrograph::Result("Problem with " + kernelname + " enqueue " + t_to_string(clState::status));
		DigitalMicrograph::Result("\n");
		return;
	}
}


clDevice::clDevice(cl_uint numDevices, cl_device_id * devices) {
	this->numDevices = numDevices;
	this->devices = devices[0];
}

clDevice::clDevice(cl_uint numDevices, cl_device_id  device) {
	this->numDevices = numDevices;
	this->devices = device;
}

cl_int clQueue::SetupQueue(cl_context & context, cl_device_id device) {
	this->cmdQueue = clCreateCommandQueue(context, device, 0, &status);

	if (!status == CL_SUCCESS) {
		DigitalMicrograph::Result("Problem with Command Queue generation \n");
		return status;
	}


	return status;
}

clQueue::clQueue(void) {
}