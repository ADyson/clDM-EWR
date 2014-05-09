#include "stdafx.h"
#include "WorkerClass.h"

//	Destructor- check to make sure we're not still running.
WorkerClass::~WorkerClass()
{
	Stop();
}

//	Start the work.
bool WorkerClass::Start()
{
	//	Check for the trivial and error case.
	if(m_status == Working)
		return true;
	if(m_status == Stopping)
		return false;

	//	The current status is 'Stopped' so we can start the collection thread.
	m_hThread = ::CreateThread(NULL, 0, reinterpret_cast<LPTHREAD_START_ROUTINE>(ThreadProcess), 
		this, NULL, &m_dwThreadID);
	if(m_hThread == NULL)
		return false;

	//	The thread has been created- we must wait until the worker class
	//	has been told explicitly it has started until we return.
	while(m_status != Working)
		Sleep(1);

	//	The status is now working, the thread is running- we're done.
	return true;
}

//	Stop the work.
bool WorkerClass::Stop()
{
	//	Check for the trivial case.
	if(m_status == Stopped)
		return true;

	//	Set the status to stopping.
	m_status = Stopping;

	//	Wait for the work to stop.
	while(m_status != Stopped)
		Sleep(1);

	//	We may as well zero the thread id- we can't use it now.
	m_dwThreadID = 0;
	m_hThread = NULL;

	//	Finally the work has stopped- so we can return.
	return true;
}

//	Terminate the thread.
bool WorkerClass::Terminate()
{
	//	Check for the trivial case.
	if(m_status == Stopped)
		return true;

	//	Set the status to stopping.
	m_status = Stopping;

	//	Kill the thread.
	::TerminateThread(m_hThread, 0);

	//	Zero the thread and return.
	m_status = Stopped;
	m_dwThreadID = 0;
	m_hThread = NULL;

	return true;
}

//	Pause the thread.
bool WorkerClass::Pause()
{
	//	Check for the trivial case.
	if(m_status == Stopped)
		return true;

	//	If we are currently not paused...
	if(m_status == Working)
	{
		//	...then pause the thread.
		::SuspendThread(m_hThread);
		m_status = Paused;
	}
	else if(m_status == Paused)
	{
		//	Or if we're paused, resume the thread.
		::ResumeThread(m_hThread);
		m_status = Working;
	}

	//	Et voila.
	return true;
}

//	The main thread process to start collection.
DWORD WorkerClass::ThreadProcess(LPVOID lpParameter)
{
	//	Cast the parameter.
	WorkerClass* pThis = reinterpret_cast<WorkerClass*>(lpParameter);

	//	Set the 'Working' flag.
	pThis->m_status = Working;

	//	Start the 'DoWork' function. The clients implement this.
	pThis->DoWork();

	//	If this function stops, then we're done working- we're stopped.
	pThis->m_status = Stopped;

	//	We're done, return success.
	return 0;
}