#pragma once

//	Copyright (C) Dave Kerr, 2009
//	dopedcoder@googlemail.com

//	This class makes implementing a Worker Class (a class with the 'DoWork' function
//	that is threaded) trivial.

//	Usage:
//
//	Start a WorkerClass by calling 'Start'.
//	Stop a WorkerClass by calling 'Stop'.
//	Pause/Unpause a WorkerClass by calling 'Pause'.
//	Make Worker Classes by deriving, and override 'DoWork'.

class WorkerClass
{
public:

	//	Represents the status of a collector.
	enum Status
	{
		Stopped,	//	The worker class isn't working.
		Stopping,	//	The worker class has been told to stop working, and is shutting down the thread.
		Paused,		//	The worker class is currently paused.
		Working		//	The worker class is working.
	};


	//	Constructor / Destructor.
	WorkerClass() : m_status(Stopped), m_dwThreadID(0), m_hThread(NULL) {}
	virtual ~WorkerClass();

	//	Stop collection. This function returns when collection has 
	//	genuinely stopped- i.e the threads stopped.
	virtual bool Stop();

	//	Kill the working function- this is not recommended but sometimes needed.
	virtual bool Terminate();

	//	Pause (or unpause!) the working function.
	virtual bool Pause();

	//	Accessor for the status.
	virtual Status GetStatus() const {return m_status;}

protected:

	//	Start collection. Derived classes should make a public implementation of
	//	start that accepts any parameters required and then calls this.
	virtual bool Start();

private:
    
	//	To implement a worker class, implement the function below.
	virtual void DoWork() = 0;

	//	This class ALWAYS handles the status- the only thing child 
	//	classes can do is check it (typically in DoWork for 'Stopping').
	volatile Status m_status;

	//	This function controls the status and makes sure that the DoWork function
	//	is properly called.
	static DWORD ThreadProcess(LPVOID lpParameter);

	//	This is the ID of the working thread (if one exists).
	DWORD m_dwThreadID;
	HANDLE m_hThread;
};