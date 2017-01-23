#ifndef SEMAPHORE_H
#define SEMAPHORE_H
// Compile with g++ thread.cpp -o thread -pthread -std=c++11

#include <mutex>
#include <condition_variable>

using namespace std;

/********************************************************************************************************************
 * Semaphore
 * This is a semaphore class that uses C++11 features.
 *******************************************************************************************************************/
class Semaphore {
public:
	Semaphore() = delete;
	Semaphore(Semaphore&& s) = delete;
	Semaphore(int c=1) : counter(c) {}
	void wait()		{unique_lock<mutex> lck(counter_mtx); while(!counter) cv.wait(lck); counter--;}
	void notify()		{unique_lock<mutex> lck(counter_mtx); counter++; cv.notify_one();}
protected:
	condition_variable	cv;
	int			counter;
	mutex			counter_mtx;
};

#endif	// SEMAPHORE_H
