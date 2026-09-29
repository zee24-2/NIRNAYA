#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/tolerances.hpp"

namespace sov {

// Regularized Cholesky / LDL^T solver for IPM normal equations and QP KKT systems (Part 5.7 & 5.8)
class CholeskySolver {
public:
    Idx n = 0;
    std::vector<double> L; // Lower triangular (column-major n x n)
    std::vector<double> D; // Diagonal (n)

    bool factor_spd(Idx dim, const std::vector<double>& S, double reg = 1e-10) {
        n = dim;
        L.assign(static_cast<size_t>(n) * n, 0.0);
        D.assign(n, 1.0);
        if (n == 0) return true;

        double max_diag = 0.0;
        for (Idx i = 0; i < n; ++i) {
            max_diag = std::max(max_diag, std::abs(S[static_cast<size_t>(i) * n + i]));
        }
        double delta = std::max(1e-12, reg * std::max(1.0, max_diag));

        for (Idx j = 0; j < n; ++j) {
            double d = S[static_cast<size_t>(j) * n + j] + delta;
            for (Idx k = 0; k < j; ++k) {
                double ljk = L[static_cast<size_t>(k) * n + j];
                d -= ljk * ljk * D[k];
            }
            if (d <= 1e-13 * std::max(1.0, max_diag)) {
                d = 1e-10 * std::max(1.0, max_diag);
            }
            D[j] = d;
            L[static_cast<size_t>(j) * n + j] = 1.0;

            for (Idx i = j + 1; i < n; ++i) {
                double sum = S[static_cast<size_t>(j) * n + i];
                for (Idx k = 0; k < j; ++k) {
                    sum -= L[static_cast<size_t>(k) * n + i] * L[static_cast<size_t>(k) * n + j] * D[k];
                }
                L[static_cast<size_t>(j) * n + i] = sum / d;
            }
        }
        return true;
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x) const {
        x = rhs;
        if (n == 0) return;
        // Forward substitution: L z = rhs
        for (Idx j = 0; j < n; ++j) {
            double zj = x[j];
            for (Idx i = j + 1; i < n; ++i) {
                x[i] -= L[static_cast<size_t>(j) * n + i] * zj;
            }
        }
        // Diagonal scaling: w = D^-1 z
        for (Idx j = 0; j < n; ++j) {
            x[j] /= D[j];
        }
        // Backward substitution: L^T x = w
        for (Idx j = n - 1; j >= 0; --j) {
            double xj = x[j];
            for (Idx i = 0; i < j; ++i) {
                x[i] -= L[static_cast<size_t>(i) * n + j] * xj;
            }
        }
    }
};

} // namespace sov
