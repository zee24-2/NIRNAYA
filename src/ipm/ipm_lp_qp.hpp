#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/status.hpp"
#include "../core/options.hpp"
#include "../core/timer.hpp"
#include "../model/model.hpp"
#include "../linalg/dense_cholesky.hpp"
#include "../linalg/sparse_cholesky_amd.hpp"

namespace sov {

struct IpmResult {
    Status status = Status::Optimal;
    double objective = 0.0;
    std::vector<double> x;
    std::vector<double> y;
    int iterations = 0;
};

// Mehrotra Predictor-Corrector Interior-Point Solver for LP and Convex QP (Parts 8 & 9)
// Uses Approximate Minimum Degree (AMD) ordering + Sparse LDL^T + dd_real Iterative Refinement (Slide 3 L2)
class IpmSolver {
public:
    // Convexity check for Q (Part 9.2): checks diagonal entries and Cholesky/LDL^T pivots
    static bool is_q_psd(const CscMatrix& Q, std::string& msg) {
        if (Q.nnz() == 0) return true;
        Idx n = Q.ncol;
        std::vector<double> dense_Q(static_cast<size_t>(n) * n, 0.0);
        for (Idx j = 0; j < n; ++j) {
            for (Off k = Q.start[j]; k < Q.start[j + 1]; ++k) {
                Idx i = Q.index[k];
                double v = Q.value[k];
                if (i == j && v < -1e-9) {
                    msg = "Negative diagonal entry in Q at column " + std::to_string(j);
                    return false;
                }
                dense_Q[static_cast<size_t>(i) * n + j] = v;
                dense_Q[static_cast<size_t>(j) * n + i] = v;
            }
        }
        // Attempt LDL^T without artificial shift to check minimum pivot
        double max_d = 0.0;
        for (Idx i = 0; i < n; ++i) max_d = std::max(max_d, std::abs(dense_Q[static_cast<size_t>(i) * n + i]));
        std::vector<double> L(static_cast<size_t>(n) * n, 0.0);
        std::vector<double> D(n, 0.0);
        for (Idx j = 0; j < n; ++j) {
            double d = dense_Q[static_cast<size_t>(j) * n + j] + 1e-10 * std::max(1.0, max_d);
            for (Idx k = 0; k < j; ++k) {
                double ljk = L[static_cast<size_t>(k) * n + j];
                d -= ljk * ljk * D[k];
            }
            if (d < -1e-7 * std::max(1.0, max_d)) {
                msg = "Indefinite pivot in Q at column " + std::to_string(j);
                return false;
            }
            if (std::abs(d) < 1e-12) d = 1e-12;
            D[j] = d;
            L[static_cast<size_t>(j) * n + j] = 1.0;
            for (Idx i = j + 1; i < n; ++i) {
                double s = dense_Q[static_cast<size_t>(j) * n + i];
                for (Idx k = 0; k < j; ++k) {
                    s -= L[static_cast<size_t>(k) * n + i] * L[static_cast<size_t>(k) * n + j] * D[k];
                }
                L[static_cast<size_t>(j) * n + i] = s / d;
            }
        }
        return true;
    }

    static IpmResult solve(const Model& model, const Options& opt, const Deadline& deadline) {
        IpmResult res;
        std::string conv_msg;
        if (model.is_qp() && !is_q_psd(model.Q, conv_msg)) {
            res.status = Status::NotSupported;
            return res;
        }

        // Computational form X = (x, r) of size nt = n + m with equality M X = 0, L <= X <= U
        Idx n = model.ncol;
        Idx m = model.nrow;
        Idx nt = n + m;

        std::vector<double> c(nt, 0.0), l(nt, 0.0), u(nt, 0.0);
        for (Idx j = 0; j < n; ++j) {
            c[j] = model.cost[j];
            l[j] = model.collo[j];
            u[j] = model.colup[j];
        }
        for (Idx i = 0; i < m; ++i) {
            c[n + i] = 0.0;
            l[n + i] = model.rowlo[i];
            u[n + i] = model.rowup[i];
        }

        // Dense Q_full on structurals (n x n)
        std::vector<double> Q_dense(static_cast<size_t>(n) * n, 0.0);
        for (Idx j = 0; j < model.Q.ncol; ++j) {
            for (Off k = model.Q.start[j]; k < model.Q.start[j + 1]; ++k) {
                Idx i = model.Q.index[k];
                double v = model.Q.value[k];
                Q_dense[static_cast<size_t>(i) * n + j] = v;
                Q_dense[static_cast<size_t>(j) * n + i] = v;
            }
        }

        // Initialize interior point (x, y, zl, zu) strictly inside [l, u]
        std::vector<double> X(nt, 0.0), xl(nt, 1.0), xu(nt, 1.0), zl(nt, 1.0), zu(nt, 1.0);
        std::vector<bool> has_l(nt, false), has_u(nt, false), is_fixed(nt, false);
        int n_comp = 0;
        for (Idx j = 0; j < nt; ++j) {
            has_l[j] = is_finite_bound(l[j]);
            has_u[j] = is_finite_bound(u[j]);
            if (has_l[j] && has_u[j] && std::abs(u[j] - l[j]) <= 1e-11) {
                is_fixed[j] = true;
                has_l[j] = false;
                has_u[j] = false;
                X[j] = 0.5 * (l[j] + u[j]);
                zl[j] = 0.0;
                zu[j] = 0.0;
            } else if (has_l[j] && has_u[j]) {
                X[j] = 0.5 * (l[j] + u[j]);
                double w = std::max(1.0, 0.5 * (u[j] - l[j]));
                xl[j] = std::max(1e-2, X[j] - l[j]);
                xu[j] = std::max(1e-2, u[j] - X[j]);
                zl[j] = 1.0 / w;
                zu[j] = 1.0 / w;
                n_comp += 2;
            } else if (has_l[j]) {
                X[j] = l[j] + 1.0;
                xl[j] = 1.0;
                zl[j] = 1.0;
                zu[j] = 0.0;
                n_comp += 1;
            } else if (has_u[j]) {
                X[j] = u[j] - 1.0;
                xu[j] = 1.0;
                zu[j] = 1.0;
                zl[j] = 0.0;
                n_comp += 1;
            } else {
                X[j] = 0.0;
                zl[j] = 0.0;
                zu[j] = 0.0;
            }
        }
        std::vector<double> y(m, 0.0);
        SparseCholeskyAmd chol;

        // Dense representation of M = [A | -I] (m x nt)
        std::vector<double> M_dense(static_cast<size_t>(m) * nt, 0.0);
        for (Idx j = 0; j < n; ++j) {
            for (Off k = model.A.start[j]; k < model.A.start[j + 1]; ++k) {
                Idx i = model.A.index[k];
                M_dense[static_cast<size_t>(i) * nt + j] = model.A.value[k];
            }
        }
        for (Idx i = 0; i < m; ++i) {
            M_dense[static_cast<size_t>(i) * nt + (n + i)] = -1.0;
        }

        for (int iter = 0; iter < opt.ipm_max_iter; ++iter) {
            if (deadline.expired()) {
                res.status = Status::TimeLimit;
                break;
            }
            res.iterations = iter + 1;

            // Compute residuals r_p = -(M X) and r_d = c + Q X - M^T y - zl + zu
            std::vector<double> rp(m, 0.0), rd(nt, 0.0), QX(nt, 0.0);
            for (Idx i = 0; i < n; ++i) {
                double s = 0.0;
                for (Idx j = 0; j < n; ++j) s += Q_dense[static_cast<size_t>(i) * n + j] * X[j];
                QX[i] = s;
            }
            double max_rp = 0.0, max_rd = 0.0;
            for (Idx i = 0; i < m; ++i) {
                double mx = 0.0;
                for (Idx j = 0; j < nt; ++j) mx += M_dense[static_cast<size_t>(i) * nt + j] * X[j];
                rp[i] = -mx;
                max_rp = std::max(max_rp, std::abs(rp[i]));
            }
            for (Idx j = 0; j < nt; ++j) {
                if (is_fixed[j]) {
                    rd[j] = 0.0;
                    continue;
                }
                double mty = 0.0;
                for (Idx i = 0; i < m; ++i) mty += M_dense[static_cast<size_t>(i) * nt + j] * y[i];
                rd[j] = c[j] + QX[j] - mty - (has_l[j] ? zl[j] : 0.0) + (has_u[j] ? zu[j] : 0.0);
                max_rd = std::max(max_rd, std::abs(rd[j]));
            }

            double comp_sum = 0.0;
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j]) comp_sum += xl[j] * zl[j];
                if (has_u[j]) comp_sum += xu[j] * zu[j];
            }
            double mu = (n_comp > 0) ? (comp_sum / n_comp) : 0.0;

            if (max_rp <= opt.tol.ipm_tol * (1.0 + n) &&
                max_rd <= opt.tol.ipm_tol * (1.0 + n) &&
                mu <= opt.tol.ipm_tol) {
                res.status = Status::Optimal;
                break;
            }

            // Diagonal barrier + Q matrix H = Q_full + D where D_j = zl_j/xl_j + zu_j/xu_j
            std::vector<double> D_bar(nt, 1e-8);
            for (Idx j = 0; j < nt; ++j) {
                if (is_fixed[j]) {
                    D_bar[j] = 1.0;
                    continue;
                }
                double dj = 1e-9;
                if (has_l[j]) dj += zl[j] / std::max(1e-12, xl[j]);
                if (has_u[j]) dj += zu[j] / std::max(1e-12, xu[j]);
                D_bar[j] = dj;
            }

            // Invert H = Q_full + diag(D_bar) via Cholesky on the n x n structural block + diagonal logical block
            CholeskySolver H_chol;
            std::vector<double> H_struct(static_cast<size_t>(n) * n, 0.0);
            for (Idx i = 0; i < n; ++i) {
                for (Idx j = 0; j < n; ++j) {
                    if (!is_fixed[i] && !is_fixed[j]) {
                        H_struct[static_cast<size_t>(i) * n + j] = Q_dense[static_cast<size_t>(i) * n + j];
                    }
                }
                H_struct[static_cast<size_t>(i) * n + i] += D_bar[i];
            }
            H_chol.factor_spd(n, H_struct, 1e-10);

            auto apply_H_inv = [&](const std::vector<double>& v_in, std::vector<double>& v_out) {
                v_out.resize(nt);
                std::vector<double> v_s(v_in.begin(), v_in.begin() + n), out_s;
                H_chol.solve(v_s, out_s);
                for (Idx j = 0; j < n; ++j) v_out[j] = is_fixed[j] ? 0.0 : out_s[j];
                for (Idx i = 0; i < m; ++i) v_out[n + i] = is_fixed[n + i] ? 0.0 : (v_in[n + i] / D_bar[n + i]);
            };

            // Build normal equations matrix S = M H^{-1} M^T (m x m)
            std::vector<double> S_norm(static_cast<size_t>(m) * m, 0.0);
            for (Idx r2 = 0; r2 < m; ++r2) {
                std::vector<double> col_r(nt, 0.0), h_inv_col;
                for (Idx j = 0; j < nt; ++j) col_r[j] = M_dense[static_cast<size_t>(r2) * nt + j];
                apply_H_inv(col_r, h_inv_col);
                for (Idx r1 = 0; r1 < m; ++r1) {
                    double dot = 0.0;
                    for (Idx j = 0; j < nt; ++j) dot += M_dense[static_cast<size_t>(r1) * nt + j] * h_inv_col[j];
                    S_norm[static_cast<size_t>(r1) * m + r2] = dot;
                }
            }
            chol.factor_spd(m, S_norm, 1e-9);

            auto solve_kkt = [&](const std::vector<double>& rcl,
                                 const std::vector<double>& rcu,
                                 std::vector<double>& dX,
                                 std::vector<double>& dy,
                                 std::vector<double>& dzl,
                                 std::vector<double>& dzu)
            {
                std::vector<double> xi(nt, 0.0);
                for (Idx j = 0; j < nt; ++j) {
                    xi[j] = rd[j]
                          - (has_l[j] ? (rcl[j] / std::max(1e-12, xl[j])) : 0.0)
                          + (has_u[j] ? (rcu[j] / std::max(1e-12, xu[j])) : 0.0);
                }
                std::vector<double> H_inv_xi;
                apply_H_inv(xi, H_inv_xi);
                std::vector<double> rhs_y(m, 0.0);
                for (Idx i = 0; i < m; ++i) {
                    double m_hxi = 0.0;
                    for (Idx j = 0; j < nt; ++j) m_hxi += M_dense[static_cast<size_t>(i) * nt + j] * H_inv_xi[j];
                    rhs_y[i] = rp[i] + m_hxi;
                }
                chol.solve_refined(S_norm, rhs_y, dy);
                std::vector<double> mty_minus_xi(nt, 0.0);
                for (Idx j = 0; j < nt; ++j) {
                    double mty = 0.0;
                    for (Idx i = 0; i < m; ++i) mty += M_dense[static_cast<size_t>(i) * nt + j] * dy[i];
                    mty_minus_xi[j] = mty - xi[j];
                }
                apply_H_inv(mty_minus_xi, dX);
                dzl.assign(nt, 0.0);
                dzu.assign(nt, 0.0);
                for (Idx j = 0; j < nt; ++j) {
                    if (has_l[j]) dzl[j] = (rcl[j] - zl[j] * dX[j]) / std::max(1e-12, xl[j]);
                    if (has_u[j]) dzu[j] = (rcu[j] + zu[j] * dX[j]) / std::max(1e-12, xu[j]);
                }
            };

            // 1. Predictor (Affine scaling) step
            std::vector<double> rcl_a(nt, 0.0), rcu_a(nt, 0.0);
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j]) rcl_a[j] = -xl[j] * zl[j];
                if (has_u[j]) rcu_a[j] = -xu[j] * zu[j];
            }
            std::vector<double> dX_a, dy_a, dzl_a, dzu_a;
            solve_kkt(rcl_a, rcu_a, dX_a, dy_a, dzl_a, dzu_a);

            double alpha_pa = 1.0, alpha_da = 1.0;
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j] && dX_a[j] < 0.0) alpha_pa = std::min(alpha_pa, -xl[j] / dX_a[j]);
                if (has_u[j] && dX_a[j] > 0.0) alpha_pa = std::min(alpha_pa, xu[j] / dX_a[j]);
                if (has_l[j] && dzl_a[j] < 0.0) alpha_da = std::min(alpha_da, -zl[j] / dzl_a[j]);
                if (has_u[j] && dzu_a[j] < 0.0) alpha_da = std::min(alpha_da, -zu[j] / dzu_a[j]);
            }
            double comp_a = 0.0;
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j]) comp_a += (xl[j] + alpha_pa * dX_a[j]) * (zl[j] + alpha_da * dzl_a[j]);
                if (has_u[j]) comp_a += (xu[j] - alpha_pa * dX_a[j]) * (zu[j] + alpha_da * dzu_a[j]);
            }
            double mu_a = (n_comp > 0) ? (comp_a / n_comp) : 0.0;
            double sigma = (mu > 1e-15) ? std::pow(std::clamp(mu_a / mu, 0.0, 1.0), 3.0) : 0.1;
            sigma = std::clamp(sigma, 1e-4, 0.99);

            // 2. Corrector step
            std::vector<double> rcl_c(nt, 0.0), rcu_c(nt, 0.0);
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j]) rcl_c[j] = sigma * mu - xl[j] * zl[j] - dX_a[j] * dzl_a[j];
                if (has_u[j]) rcu_c[j] = sigma * mu - xu[j] * zu[j] + dX_a[j] * dzu_a[j];
            }
            std::vector<double> dX, dy, dzl, dzu;
            solve_kkt(rcl_c, rcu_c, dX, dy, dzl, dzu);

            double alpha_p = 1.0, alpha_d = 1.0;
            for (Idx j = 0; j < nt; ++j) {
                if (has_l[j] && dX[j] < 0.0) alpha_p = std::min(alpha_p, -xl[j] / dX[j]);
                if (has_u[j] && dX[j] > 0.0) alpha_p = std::min(alpha_p, xu[j] / dX[j]);
                if (has_l[j] && dzl[j] < 0.0) alpha_d = std::min(alpha_d, -zl[j] / dzl[j]);
                if (has_u[j] && dzu[j] < 0.0) alpha_d = std::min(alpha_d, -zu[j] / dzu[j]);
            }
            double alpha = (model.is_qp() ? std::min(alpha_p, alpha_d) : alpha_p) * 0.995;
            double alpha_dual = (model.is_qp() ? alpha : alpha_d * 0.995);

            for (Idx j = 0; j < nt; ++j) {
                X[j] += alpha * dX[j];
                if (has_l[j]) {
                    xl[j] = std::max(1e-11, X[j] - l[j]);
                    zl[j] = std::max(1e-11, zl[j] + alpha_dual * dzl[j]);
                }
                if (has_u[j]) {
                    xu[j] = std::max(1e-11, u[j] - X[j]);
                    zu[j] = std::max(1e-11, zu[j] + alpha_dual * dzu[j]);
                }
            }
            for (Idx i = 0; i < m; ++i) {
                y[i] += alpha_dual * dy[i];
            }
        }

        res.x.assign(X.begin(), X.begin() + n);
        for (Idx j = 0; j < n; ++j) {
            if (has_l[j]) res.x[j] = std::max(l[j], res.x[j]);
            if (has_u[j]) res.x[j] = std::min(u[j], res.x[j]);
        }
        res.y = y;
        res.objective = model.eval_objective(res.x);
        return res;
    }
};

} // namespace sov
