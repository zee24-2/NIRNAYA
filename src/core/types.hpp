#pragma once
#include <cstdint>
#include <limits>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

namespace sov {

using Idx = int32_t;
using Off = int64_t;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kInfInput = 1e20;
constexpr Idx kNone = -1;
constexpr double kTinySentinel = 1e-50;

enum class VarType : int8_t {
    Continuous = 0,
    Integer    = 1,
    Binary     = 2,
    SemiCont   = 3
};

enum class Sense : int8_t {
    Minimize = 1,
    Maximize = -1
};

enum class RowType : int8_t {
    N = 0, // Free
    L = 1, // <=
    G = 2, // >=
    E = 3, // ==
    Ranged = 4
};

enum class VarStatus : int8_t {
    Basic      = 0,
    AtLower    = 1,
    AtUpper    = 2,
    Fixed      = 3,
    FreeZero   = 4,
    SuperBasic = 5
};

inline bool is_finite_bound(double b) {
    return std::isfinite(b) && std::abs(b) < kInfInput;
}

inline double sanitize_bound_lo(double b) {
    if (!std::isfinite(b) || b <= -kInfInput) return -kInf;
    return b;
}

inline double sanitize_bound_up(double b) {
    if (!std::isfinite(b) || b >= kInfInput) return kInf;
    return b;
}

} // namespace sov
