#pragma once
#include <mutex>
#include <condition_variable>

namespace aadc {
class AADCFunctionSynchronizer {
public:

    AADCFunctionSynchronizer(int num_threads)
        : _num_threads(num_threads)
        , _counter(0)
        , _num_threads_at_this_sync(0)
    {}

    void syncThreads() {
        int waiting_at = _counter;
        {
            std::lock_guard<std::mutex> lock(mutex);
            ++_num_threads_at_this_sync;
            if (_num_threads_at_this_sync >= _num_threads) {
                _counter++;
                _num_threads_at_this_sync = 0;
                cv.notify_all();
            }
        }

        std::unique_lock<std::mutex> lk(mutex);
        cv.wait(lk, [&]{return _counter > waiting_at;});
    }

private:
    const int _num_threads;
    int _counter;
    int _num_threads_at_this_sync;
    std::condition_variable cv;
    std::mutex mutex;
};

};// namespace aadc