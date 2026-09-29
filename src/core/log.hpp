#pragma once
#include <cstdio>
#include <cstdarg>
#include <mutex>

namespace sov {

class Logger {
public:
    explicit Logger(int level = 1) : level_(level) {}

    void set_level(int l) { level_ = l; }
    int level() const { return level_; }

    void log(int min_level, const char* fmt, ...) const {
        if (level_ < min_level) return;
        std::lock_guard<std::mutex> lock(mu_);
        va_list args;
        va_start(args, fmt);
        std::vprintf(fmt, args);
        va_end(args);
        std::fflush(stdout);
    }

private:
    int level_ = 1;
    mutable std::mutex mu_;
};

} // namespace sov
