#include <functional>
#include <thread>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <mutex> // to enable multiple threads to access RHS one at the time!

#include "AadDataNew.hpp"
#include "OdeDriverNew.hpp"

class ThreadPool
{
public:
    void Start()
    {
        const uint32_t num_threads = std::thread::hardware_concurrency(); // Max # of threads the system supports
        threads.resize(num_threads);
        for (uint32_t i = 0; i < num_threads; i++)
        {
            std::thread th(&ThreadPool::ThreadLoop, this);
            threads.push_back(std::move(th));
        }
    }
    void QueueJob(const std::function<void()> &job)
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            jobs.push(job);
        }
        mutex_condition.notify_one();
    }
    void Stop()
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            should_terminate = true;
        }
        mutex_condition.notify_all();
        for (std::thread &active_thread : threads)
        {
            active_thread.join();
        }
        threads.clear();
    }
    // Busy:
    bool Busy()
    {
        bool poolbusy;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            poolbusy = jobs.empty();
        }
        return poolbusy;
    }

private:
    // ThreadLoop:
    void ThreadLoop()
    {
        while (true)
        {
            std::function<void()> job;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                mutex_condition.wait(lock, [this]
                                     { return !jobs.empty() || should_terminate; });
                if (should_terminate)
                {
                    return;
                }
                job = jobs.front();
                jobs.pop();
            }
            job();
        }
    }

    bool should_terminate = false;           // Tells threads to stop looking for jobs
    std::mutex queue_mutex;                  // Prevents data races to the job queue
    std::condition_variable mutex_condition; // Allows threads to wait on new jobs or termination
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> jobs;
};
