#pragma once
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <future>
#include <atomic>
#include "types.hpp"

namespace sov {

// Layer L1 Core: Native C++17 ThreadPool for parallel portfolio, strong branching, and PDHG kernels
class ThreadPool {
public:
    explicit ThreadPool(int num_threads = 1) : stop_(false) {
        int n = std::max(1, num_threads);
        for (int i = 0; i < n; ++i) {
            workers_.emplace_back([this]() {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(mu_);
                        cv_.wait(lock, [this]() { return stop_ || !tasks_.empty(); });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto& w : workers_) {
            if (w.joinable()) w.join();
        }
    }

    template <typename F>
    auto submit(F&& f) -> std::future<decltype(f())> {
        using Ret = decltype(f());
        auto task = std::make_shared<std::packaged_task<Ret()>>(std::forward<F>(f));
        std::future<Ret> res = task->get_future();
        {
            std::lock_guard<std::mutex> lock(mu_);
            tasks_.emplace([task]() { (*task)(); });
        }
        cv_.notify_one();
        return res;
    }

    void parallel_for(Idx begin, Idx end, const std::function<void(Idx)>& body) {
        Idx total = end - begin;
        if (total <= 0) return;
        int nw = static_cast<int>(workers_.size());
        if (nw <= 1 || total < 64) {
            for (Idx i = begin; i < end; ++i) body(i);
            return;
        }
        Idx chunk = (total + nw - 1) / nw;
        std::vector<std::future<void>> futs;
        for (int t = 0; t < nw; ++t) {
            Idx b = begin + t * chunk;
            Idx e = std::min(end, b + chunk);
            if (b >= e) break;
            futs.push_back(submit([b, e, &body]() {
                for (Idx i = b; i < e; ++i) body(i);
            }));
        }
        for (auto& f : futs) f.get();
    }

    int size() const { return static_cast<int>(workers_.size()); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool stop_;
};

} // namespace sov
