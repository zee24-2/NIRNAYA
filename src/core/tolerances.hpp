#pragma once
#include <cmath>
#include <algorithm>
#include "types.hpp"

namespace sov {

struct Tolerances {
    double primal_feas_tol       = 1e-7;
    double dual_feas_tol         = 1e-7;
    double int_tol               = 1e-6;
    double zero_tol              = 1e-12;
    double drop_tol              = 1e-14;
    double pivot_abs_tol         = 1e-9;
    double pivot_rel_tol         = 1e-7;
    double lu_threshold_u        = 0.1;
    double lu_singular_abs_tol   = 1e-9;
    double ft_stability_tol      = 1e-8;
    double dual_shift_limit      = 1e-6;
    double perturb_base          = 5e-7;
    double dse_min_weight        = 1e-4;
    double ipm_tol               = 1e-8;
    double ipm_step_eta          = 0.9995;
    double ipm_reg_primal        = 1e-8;
    double ipm_reg_dual          = 1e-8;
    double cholesky_pivot_tol    = 1e-13;
    double mip_gap               = 1e-4;
    double mip_abs_gap           = 1e-6;
    double cut_min_efficacy_root = 1e-4;
    double cut_min_efficacy_tree = 1e-3;
    double cut_max_dynamism      = 1e6;
    double cut_parallelism_limit = 0.95;
    double checker_tol           = 1e-6;
};

constexpr double kZeroTol = 1e-12;
constexpr double kIntTol  = 1e-6;
constexpr double kDropTol = 1e-14;

inline bool is_zero(double v, double tol = kZeroTol) {
    return std::abs(v) <= tol;
}

inline bool feas_le(double a, double b, double tol) {
    return a <= b + tol;
}

inline bool feas_ge(double a, double b, double tol) {
    return a >= b - tol;
}

inline bool rel_le(double a, double b, double tol) {
    return a <= b + tol * (1.0 + std::abs(b));
}

inline bool rel_ge(double a, double b, double tol) {
    return a >= b - tol * (1.0 + std::abs(b));
}

inline bool is_integral(double v, double tol = kIntTol) {
    return std::abs(v - std::round(v)) <= tol;
}

inline double safe_div(double a, double b) {
    if (b == 0.0) {
        if (a == 0.0) return 0.0;
        return (a > 0.0) ? kInf : -kInf;
    }
    return a / b;
}

inline double frac(double v) {
    return v - std::floor(v);
}

} // namespace sov
