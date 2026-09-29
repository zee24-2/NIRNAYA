#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "model.hpp"

namespace sov {

struct ScalingInfo {
    bool active = false;
    std::vector<double> R; // Row scales (m)
    std::vector<double> C; // Col scales (n)

    static double round_pow2(double v) {
        if (v <= 0.0 || !std::isfinite(v)) return 1.0;
        int exp = static_cast<int>(std::round(std::log2(v)));
        exp = std::max(-30, std::min(30, exp));
        return std::ldexp(1.0, exp);
    }

    void unscale_primal(const std::vector<double>& x_scaled, std::vector<double>& x_orig) const {
        x_orig.resize(x_scaled.size());
        if (!active) {
            x_orig = x_scaled;
            return;
        }
        for (size_t j = 0; j < x_scaled.size(); ++j) {
            x_orig[j] = x_scaled[j] * C[j];
        }
    }

    void unscale_row_duals(const std::vector<double>& y_scaled, std::vector<double>& y_orig) const {
        y_orig.resize(y_scaled.size());
        if (!active) {
            y_orig = y_scaled;
            return;
        }
        for (size_t i = 0; i < y_scaled.size(); ++i) {
            y_orig[i] = y_scaled[i] * R[i];
        }
    }

    void unscale_reduced_costs(const std::vector<double>& d_scaled, std::vector<double>& d_orig) const {
        d_orig.resize(d_scaled.size());
        if (!active) {
            d_orig = d_scaled;
            return;
        }
        for (size_t j = 0; j < d_scaled.size(); ++j) {
            d_orig[j] = d_scaled[j] / C[j];
        }
    }
};

// Iterated geometric mean + equilibration scaling rounded to powers of two (Part 4)
inline Model scale_model(const Model& orig, bool enable_scaling, ScalingInfo& info) {
    Idx m = orig.nrow;
    Idx n = orig.ncol;
    info.R.assign(m, 1.0);
    info.C.assign(n, 1.0);
    info.active = false;

    if (!enable_scaling || orig.A.nnz() == 0) {
        return orig;
    }

    double max_a = 0.0, min_a = kInf;
    for (double v : orig.A.value) {
        double av = std::abs(v);
        if (av > 0.0) {
            max_a = std::max(max_a, av);
            min_a = std::min(min_a, av);
        }
    }
    if (min_a >= kInf || (max_a / min_a) <= 10.0) {
        return orig;
    }

    info.active = true;
    int passes = ((max_a / min_a) < 1e3) ? 2 : 6;
    std::vector<double> R(m, 1.0), C(n, 1.0);

    for (int pass = 0; pass < passes; ++pass) {
        double max_change = 1.0;
        std::vector<double> row_min(m, kInf), row_max(m, 0.0);
        for (Idx j = 0; j < n; ++j) {
            for (Off k = orig.A.start[j]; k < orig.A.start[j + 1]; ++k) {
                Idx i = orig.A.index[k];
                double v = std::abs(orig.A.value[k]) * R[i] * C[j];
                if (v > 0.0) {
                    row_min[i] = std::min(row_min[i], v);
                    row_max[i] = std::max(row_max[i], v);
                }
            }
        }
        for (Idx i = 0; i < m; ++i) {
            if (row_max[i] > 0.0 && row_min[i] < kInf) {
                double f = 1.0 / std::sqrt(row_min[i] * row_max[i]);
                max_change = std::max(max_change, std::max(f, 1.0 / f));
                R[i] *= f;
            }
        }

        std::vector<double> col_min(n, kInf), col_max(n, 0.0);
        for (Idx j = 0; j < n; ++j) {
            if (orig.vartype[j] != VarType::Continuous) continue; // Keep C_j = 1 for integers!
            for (Off k = orig.A.start[j]; k < orig.A.start[j + 1]; ++k) {
                Idx i = orig.A.index[k];
                double v = std::abs(orig.A.value[k]) * R[i] * C[j];
                if (v > 0.0) {
                    col_min[j] = std::min(col_min[j], v);
                    col_max[j] = std::max(col_max[j], v);
                }
            }
            if (col_max[j] > 0.0 && col_min[j] < kInf) {
                double f = 1.0 / std::sqrt(col_min[j] * col_max[j]);
                max_change = std::max(max_change, std::max(f, 1.0 / f));
                C[j] *= f;
            }
        }
        if (max_change < 1.05) break;
    }

    for (Idx i = 0; i < m; ++i) info.R[i] = ScalingInfo::round_pow2(R[i]);
    for (Idx j = 0; j < n; ++j) {
        info.C[j] = (orig.vartype[j] == VarType::Continuous) ? ScalingInfo::round_pow2(C[j]) : 1.0;
    }

    Model scaled = orig;
    for (Idx j = 0; j < n; ++j) {
        double cj = info.C[j];
        scaled.cost[j] = orig.cost[j] * cj;
        if (is_finite_bound(orig.collo[j])) scaled.collo[j] = orig.collo[j] / cj;
        if (is_finite_bound(orig.colup[j])) scaled.colup[j] = orig.colup[j] / cj;
        for (Off k = scaled.A.start[j]; k < scaled.A.start[j + 1]; ++k) {
            Idx i = scaled.A.index[k];
            scaled.A.value[k] = orig.A.value[k] * info.R[i] * cj;
        }
    }
    for (Idx i = 0; i < m; ++i) {
        double ri = info.R[i];
        if (is_finite_bound(orig.rowlo[i])) scaled.rowlo[i] = orig.rowlo[i] * ri;
        if (is_finite_bound(orig.rowup[i])) scaled.rowup[i] = orig.rowup[i] * ri;
    }
    if (scaled.Q.nnz() > 0) {
        for (Idx j = 0; j < n; ++j) {
            for (Off k = scaled.Q.start[j]; k < scaled.Q.start[j + 1]; ++k) {
                Idx i = scaled.Q.index[k];
                scaled.Q.value[k] = orig.Q.value[k] * info.C[i] * info.C[j];
            }
        }
    }
    return scaled;
}

} // namespace sov
