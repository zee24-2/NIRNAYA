#pragma once
#include <chrono>
#include <atomic>
#include <cstdint>
#include "types.hpp"

namespace sov {

class WorkCounter {
public:
    void add(int64_t k) { count_ += k; }
    int64_t count() const { return count_; }
    void reset() { count_ = 0; }
private:
    int64_t count_ = 0;
};

class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    void reset() {
        start_ = std::chrono::steady_clock::now();
    }

    double elapsed_sec() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - start_).count();
    }
private:
    std::chrono::steady_clock::time_point start_;
};

struct Deadline {
    Timer timer;
    double wall_limit_sec = kInf;
    double work_limit     = kInf;
    WorkCounter work;
    std::atomic<bool>* stop_flag = nullptr;

    bool expired() const {
        if (stop_flag && stop_flag->load(std::memory_order_relaxed)) return true;
        if (timer.elapsed_sec() >= wall_limit_sec) return true;
        if (static_cast<double>(work.count()) >= work_limit) return true;
        return false;
    }

    bool is_interrupted() const {
        return stop_flag && stop_flag->load(std::memory_order_relaxed);
    }

    bool is_work_expired() const {
        return static_cast<double>(work.count()) >= work_limit;
    }
};

} // namespace sov
