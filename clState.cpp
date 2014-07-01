#include "StdAfx.h"
#include "clState.h"
#include "boost\lexical_cast.hpp"


clState::clState(void)
{
}


clState::~clState(void)
{
}

// initialise statics?
cl_int clState::status = 0;
cl_context clState::context = NULL;
clDevice* clState::cldev = NULL;
clQueue* clState::clq = NULL;
bool clState::OpenCLAvailable = false;

// Call this once somewhere during plugin load.
void clState::Setup()
{
	
	// TODO: Check for persistent tags to get platform and devnumber
	DigitalMicrograph::TagGroup PersistentTags = DigitalMicrograph::GetPersistentTagGroup();
	DigitalMicrograph::String clPlatformTag;
	DigitalMicrograph::String clDeviceTag;


	DigitalMicrograph::TagGroupGetTagAsString(PersistentTags,"OpenCL:Platform",clPlatformTag);
	DigitalMicrograph::TagGroupGetTagAsString(PersistentTags,"OpenCL:Device",clDeviceTag);

	int platformnumber;
	int devicenumber;

	try
	{
		platformnumber = boost::lexical_cast<int>(clPlatformTag);
		devicenumber = boost::lexical_cast<int>(clDeviceTag);
	}
	catch(boost::bad_lexical_cast e)
	{
		DigitalMicrograph::Result("Open CL Device and Platform incorrectly set in global tags \n");
		DigitalMicrograph::Result("Set Platform and Device 0,1,2 etc in tags OpenCl:Platform and OpenCL:Device \n");
		DigitalMicrograph::Result("Run clinfo at command prompt to find available devices and platforms");
		return;
	}

	cl_uint numDevices;
	cl_device_id* devices;

	//Setup OpenCL
	context = NULL;
	numDevices = 0;
	devices = NULL;

	// Maybe Can Do OpenCL setup and device registering here - Print to Ouput with device data?
	// Discover and initialize available platforms
	cl_uint numPlatforms = 0;
	cl_platform_id * platforms = NULL;

	// Use clGetPlatformIds() to retrieve the number of platforms
	status = clGetPlatformIDs(0,NULL,&numPlatforms);

	// Allocate enough space for each platform
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));

	// Fill in platforms with clGetPlatformIDs()
	status = clGetPlatformIDs(numPlatforms,platforms,NULL);

	// Discover and initialize available devices	
	// use clGetDeviceIDs() to retrieve number of devices present
	status = clGetDeviceIDs(platforms[platformnumber],CL_DEVICE_TYPE_ALL,0,NULL,&numDevices);

	// Allocate enough space for each device
	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));

	// Fill in devices with clGetDeviceIDs()
	status = clGetDeviceIDs(platforms[platformnumber],CL_DEVICE_TYPE_ALL,numDevices,devices,NULL);

	// Most of initialisation is done, would be nice to print device information...
	//Getting the device name
	size_t deviceNameLength = 4096;
	size_t actualSize;
	char* tempDeviceName = (char*)malloc(4096);
	char* deviceName;
	status |= clGetDeviceInfo(devices[devicenumber], CL_DEVICE_NAME, deviceNameLength, tempDeviceName, &actualSize);

	if(status == CL_SUCCESS)
	{
		deviceName = (char*)malloc(actualSize);
		memcpy(deviceName, tempDeviceName, actualSize);
		free(tempDeviceName);

		std::string devName(deviceName);
		DigitalMicrograph::Result("Using OpenCL on device "+devName+" - OCL - EWR\n");
		DigitalMicrograph::Result("To change edit the Global Tags OpenCL:Platform and OpenCL:Device then restart DM\n");
		OpenCLAvailable = true;

		context = clCreateContext(NULL,1,&devices[devicenumber],NULL,NULL,&status);
		clq = new clQueue();
		clq->SetupQueue(context,devices[devicenumber]);
		cldev = new clDevice(1,&devices[devicenumber]);
	}
	if(status!=CL_SUCCESS)
	{
		DigitalMicrograph::Result("Could not setup OpenCL on this computer, run clinfo to check availability\n");
	}
}

cl_context clState::GetContext()
{
	return context;
}

clDevice* clState::GetDevicePtr()
{
	return cldev;
}

clQueue* clState::GetQueuePtr()
{
	return clq;
}

cl_int clState::GetStatus()
{
	return status;
}