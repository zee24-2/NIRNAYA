#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/status.hpp"
#include "../core/tolerances.hpp"
#include "../core/options.hpp"
#include "../core/timer.hpp"
#include "../core/rng.hpp"
#include "../core/hvector.hpp"
#include "../core/csc_matrix.hpp"
#include "../model/model.hpp"
#include "../linalg/lu_factor.hpp"

namespace sov {

struct BasisSnapshot {
    std::vector<VarStatus> status;
    std::vector<Idx> basic_index;
};

struct SparseRow {
    double lo = -kInf;
    double up = kInf;
    std::vector<Idx> idx;
    std::vector<double> val;
};

struct LpStats {
    Idx iterations = 0;
    Idx phase1_iterations = 0;
    Idx phase2_iterations = 0;
    Idx refactors = 0;
    Idx bound_flips = 0;
};

// Full Dual & Primal Simplex Engine with Sparse LU, DSE, Harris, BFRT, and Warm-Start API (Part 6)
class LpSolver {
public:
    Idx m = 0;
    Idx n = 0;
    Idx nt = 0;
    MatrixPair A;
    std::vector<double> cost;
    std::vector<double> cost_orig;
    std::vector<double> cost_shift;
    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<double> value;
    std::vector<double> dual;      // Reduced costs d_j (0 for basic)
    std::vector<double> row_duals; // y = B^{-T} c_B
    std::vector<VarStatus> status;
    std::vector<Idx> basic_index;
    std::vector<Idx> pos_in_basis;
    std::vector<double> primal_infeas;
    std::vector<double> dse_weight;
    std::vector<double> devex_weight;

    double obj_offset = 0.0;
    double objective = 0.0;
    double objective_cutoff = kInf;
    BasisFactor factor;
    Options opt;
    LpStats stats;
    Rng rng;

    HVector rho, alpha_row, alpha_col, spike, tau_vec, flip_rhs;
    std::vector<double> farkas_ray;
    std::vector<double> primal_ray;

    explicit LpSolver(const Options& options = Options())
        : opt(options), rng(options.seed) {}

    void load(const Model& model) {
        m = model.nrow;
        n = model.ncol;
        nt = n + m;
        obj_offset = model.obj_offset;
        A.build(model.A);

        cost.assign(nt, 0.0);
        lower.assign(nt, 0.0);
        upper.assign(nt, 0.0);
        for (Idx j = 0; j < n; ++j) {
            cost[j] = model.cost[j];
            lower[j] = model.collo[j];
            upper[j] = model.colup[j];
        }
        for (Idx i = 0; i < m; ++i) {
            cost[n + i] = 0.0;
            lower[n + i] = model.rowlo[i];
            upper[n + i] = model.rowup[i];
        }
        cost_orig = cost;
        cost_shift.assign(nt, 0.0);
        value.assign(nt, 0.0);
        dual.assign(nt, 0.0);
        row_duals.assign(m, 0.0);
        status.assign(nt, VarStatus::AtLower);
        basic_index.resize(m);
        pos_in_basis.assign(nt, kNone);
        primal_infeas.assign(m, 0.0);
        dse_weight.assign(m, 1.0);
        devex_weight.assign(nt, 1.0);

        rho.resize(m);
        alpha_row.resize(nt);
        alpha_col.resize(m);
        spike.resize(m);
        tau_vec.resize(m);
        flip_rhs.resize(m);

        factor.threshold_u = opt.tol.lu_threshold_u;
        factor.max_updates = opt.refactor_interval;

        init_slack_basis();
    }

    void init_slack_basis() {
        std::fill(pos_in_basis.begin(), pos_in_basis.end(), kNone);
        for (Idx j = 0; j < n; ++j) {
            if (is_finite_bound(lower[j]) && is_finite_bound(upper[j]) && lower[j] == upper[j]) {
                status[j] = VarStatus::Fixed;
                value[j] = lower[j];
            } else if (is_finite_bound(lower[j]) && is_finite_bound(upper[j])) {
                status[j] = (std::abs(lower[j]) <= std::abs(upper[j])) ? VarStatus::AtLower : VarStatus::AtUpper;
                value[j] = (status[j] == VarStatus::AtLower) ? lower[j] : upper[j];
            } else if (is_finite_bound(lower[j])) {
                status[j] = VarStatus::AtLower;
                value[j] = lower[j];
            } else if (is_finite_bound(upper[j])) {
                status[j] = VarStatus::AtUpper;
                value[j] = upper[j];
            } else {
                status[j] = VarStatus::FreeZero;
                value[j] = 0.0;
            }
        }
        for (Idx i = 0; i < m; ++i) {
            Idx log_var = n + i;
            status[log_var] = VarStatus::Basic;
            basic_index[i] = log_var;
            pos_in_basis[log_var] = i;
            dse_weight[i] = 1.0;
        }
    }

    void extract_col_of_M(Idx var, HVector& col) const {
        col.clear();
        if (var < n) {
            for (Off k = A.csc.start[var]; k < A.csc.start[var + 1]; ++k) {
                col.add(A.csc.index[k], A.csc.value[k]);
            }
        } else {
            col.add(var - n, -1.0);
        }
        col.tighten(0.0);
    }

    void sync_nonbasic_values() {
        for (Idx j = 0; j < nt; ++j) {
            if (status[j] == VarStatus::Basic) continue;
            if (lower[j] == upper[j] && is_finite_bound(lower[j])) {
                status[j] = VarStatus::Fixed;
                value[j] = lower[j];
            } else if (status[j] == VarStatus::AtLower) {
                if (is_finite_bound(lower[j])) value[j] = lower[j];
                else if (is_finite_bound(upper[j])) { status[j] = VarStatus::AtUpper; value[j] = upper[j]; }
                else { status[j] = VarStatus::FreeZero; value[j] = 0.0; }
            } else if (status[j] == VarStatus::AtUpper) {
                if (is_finite_bound(upper[j])) value[j] = upper[j];
                else if (is_finite_bound(lower[j])) { status[j] = VarStatus::AtLower; value[j] = lower[j]; }
                else { status[j] = VarStatus::FreeZero; value[j] = 0.0; }
            } else if (status[j] == VarStatus::Fixed) {
                value[j] = lower[j];
            } else {
                value[j] = 0.0;
            }
        }
    }

    void refactor_and_recompute() {
        auto fres = factor.factor(A.csc, basic_index);
        ++stats.refactors;
        if (!fres.repairs.empty()) {
            for (const auto& rep : fres.repairs) {
                Idx out_v = rep.out_var;
                Idx in_v = n + rep.in_logical_row;
                pos_in_basis[out_v] = kNone;
                if (is_finite_bound(lower[out_v]) && is_finite_bound(upper[out_v]) && lower[out_v] == upper[out_v]) {
                    status[out_v] = VarStatus::Fixed;
                    value[out_v] = lower[out_v];
                } else if (is_finite_bound(lower[out_v])) {
                    status[out_v] = VarStatus::AtLower;
                    value[out_v] = lower[out_v];
                } else if (is_finite_bound(upper[out_v])) {
                    status[out_v] = VarStatus::AtUpper;
                    value[out_v] = upper[out_v];
                } else {
                    status[out_v] = VarStatus::FreeZero;
                    value[out_v] = 0.0;
                }
                basic_index[rep.basis_pos] = in_v;
                pos_in_basis[in_v] = rep.basis_pos;
                status[in_v] = VarStatus::Basic;
            }
            factor.factor(A.csc, basic_index);
        }
        compute_primal();
        compute_dual();
    }

    // Compute X_B = -B^{-1} N X_N (6.4.1)
    void compute_primal() {
        sync_nonbasic_values();
        HVector rhs(m);
        for (Idx j = 0; j < n; ++j) {
            if (status[j] == VarStatus::Basic) continue;
            double xj = value[j];
            if (xj == 0.0) continue;
            for (Off k = A.csc.start[j]; k < A.csc.start[j + 1]; ++k) {
                rhs.add(A.csc.index[k], -A.csc.value[k] * xj);
            }
        }
        for (Idx i = 0; i < m; ++i) {
            Idx log_v = n + i;
            if (status[log_v] == VarStatus::Basic) continue;
            double r_val = value[log_v];
            if (r_val != 0.0) {
                // Column of logical is -e_i, so -a_j * r_val = +r_val * e_i
                rhs.add(i, r_val);
            }
        }
        rhs.tighten();
        factor.ftran(rhs);
        for (Idx k = 0; k < m; ++k) {
            Idx var = basic_index[k];
            value[var] = rhs.val[k];
        }
        update_primal_infeas_all();
        recompute_objective();
    }

    void recompute_objective() {
        double sum = obj_offset;
        for (Idx j = 0; j < n; ++j) {
            sum += cost_orig[j] * value[j];
        }
        objective = sum;
    }

    void update_primal_infeas_pos(Idx k) {
        Idx var = basic_index[k];
        double x = value[var];
        double eps_p = opt.tol.primal_feas_tol;
        double inf_v = 0.0;
        if (x < lower[var] - eps_p) {
            inf_v = lower[var] - x;
        } else if (x > upper[var] + eps_p) {
            inf_v = x - upper[var];
        }
        primal_infeas[k] = inf_v * inf_v;
    }

    void update_primal_infeas_all() {
        for (Idx k = 0; k < m; ++k) {
            update_primal_infeas_pos(k);
        }
    }

    // Compute y = B^{-T} c_B and d_j = c_j - a_j^T y (6.4.2)
    void compute_dual() {
        HVector cb(m);
        for (Idx k = 0; k < m; ++k) {
            double c = cost[basic_index[k]];
            if (c != 0.0) cb.add(k, c);
        }
        cb.tighten();
        factor.btran(cb);
        row_duals.assign(m, 0.0);
        for (Idx k = 0; k < cb.count; ++k) {
            Idx i = cb.idx[k];
            row_duals[i] = cb.val[i];
        }
        for (Idx j = 0; j < n; ++j) {
            if (status[j] == VarStatus::Basic) {
                dual[j] = 0.0;
            } else {
                double ay = 0.0;
                for (Off k = A.csc.start[j]; k < A.csc.start[j + 1]; ++k) {
                    ay += A.csc.value[k] * row_duals[A.csc.index[k]];
                }
                dual[j] = cost[j] - ay;
            }
        }
        for (Idx i = 0; i < m; ++i) {
            Idx log_v = n + i;
            if (status[log_v] == VarStatus::Basic) {
                dual[log_v] = 0.0;
            } else {
                // a_{n+i} = -e_i => d_{n+i} = cost_{n+i} - (-y_i) = cost_{n+i} + y_i
                dual[log_v] = cost[log_v] + row_duals[i];
            }
        }
    }

    // Place boxed nonbasics to satisfy dual feasibility and check remaining dual infeasibilities (6.4.3)
    Idx make_boxed_dual_feasible() {
        double eps_d = opt.tol.dual_feas_tol;
        bool changed_primal = false;
        Idx num_dual_infeas = 0;

        for (Idx j = 0; j < nt; ++j) {
            if (status[j] == VarStatus::Basic || status[j] == VarStatus::Fixed) continue;
            bool has_lo = is_finite_bound(lower[j]);
            bool has_up = is_finite_bound(upper[j]);
            double d = dual[j];

            if (has_lo && has_up) {
                if (status[j] == VarStatus::AtLower && d < -eps_d) {
                    status[j] = VarStatus::AtUpper;
                    value[j] = upper[j];
                    changed_primal = true;
                } else if (status[j] == VarStatus::AtUpper && d > eps_d) {
                    status[j] = VarStatus::AtLower;
                    value[j] = lower[j];
                    changed_primal = true;
                }
            } else if (has_lo && !has_up) {
                status[j] = VarStatus::AtLower;
                value[j] = lower[j];
                if (d < -eps_d) ++num_dual_infeas;
            } else if (!has_lo && has_up) {
                status[j] = VarStatus::AtUpper;
                value[j] = upper[j];
                if (d > eps_d) ++num_dual_infeas;
            } else {
                status[j] = VarStatus::FreeZero;
                value[j] = 0.0;
                if (std::abs(d) > eps_d) ++num_dual_infeas;
            }
        }
        if (changed_primal) {
            compute_primal();
        }
        return num_dual_infeas;
    }

    // Compute rho^T a_j for all nonbasic j (STEP 3 PRICE)
    void price_pivot_row(const HVector& rho_vec, HVector& out_alpha_row) const {
        out_alpha_row.clear();
        for (Idx k = 0; k < rho_vec.count; ++k) {
            Idx i = rho_vec.idx[k];
            double ri = rho_vec.val[i];
            if (ri == 0.0) continue;
            for (Off p = A.csr.start[i]; p < A.csr.start[i + 1]; ++p) {
                Idx j = A.csr.index[p];
                if (status[j] != VarStatus::Basic) {
                    out_alpha_row.add(j, ri * A.csr.value[p]);
                }
            }
            Idx log_v = n + i;
            if (status[log_v] != VarStatus::Basic) {
                out_alpha_row.add(log_v, -ri);
            }
        }
        out_alpha_row.tighten(1e-13);
    }

    // Core Dual Simplex Loop (Part 6.6)
    Status run_dual_phase2(const Deadline& deadline, Idx max_iters, bool enforce_cutoff) {
        double eps_p = opt.tol.primal_feas_tol;
        double eps_d = opt.tol.dual_feas_tol;
        double piv_tol = opt.tol.pivot_abs_tol;

        for (Idx iter = 0; iter < max_iters; ++iter) {
            if ((iter & 63) == 0 && deadline.expired()) {
                return deadline.is_interrupted() ? Status::Interrupted : Status::TimeLimit;
            }
            if (factor.need_refactor()) {
                refactor_and_recompute();
            }

            // Cutoff check (for MIP node pruning, 6.12)
            if (enforce_cutoff && objective_cutoff < kInf) {
                if (objective >= objective_cutoff + 1e-7 * (1.0 + std::abs(objective_cutoff))) {
                    return Status::CutoffReached;
                }
            }

            // STEP 1: CHUZR (Choose leaving basic variable with maximum DSE score)
            Idx r = kNone;
            double best_score = 0.0;
            for (Idx k = 0; k < m; ++k) {
                double inf_sq = primal_infeas[k];
                if (inf_sq > eps_p * eps_p) {
                    double score = inf_sq / std::max(1e-4, dse_weight[k]);
                    if (score > best_score) {
                        best_score = score;
                        r = k;
                    }
                }
            }

            if (r == kNone) {
                // Recompute from scratch to be 100% certain of primal/dual feasibility
                refactor_and_recompute();
                double max_p_viol = 0.0;
                for (Idx k = 0; k < m; ++k) {
                    max_p_viol = std::max(max_p_viol, std::sqrt(primal_infeas[k]));
                }
                if (max_p_viol <= eps_p * 2.0) {
                    return Status::Optimal;
                }
                continue;
            }

            Idx p = basic_index[r];
            double xp = value[p];
            double s = (xp < lower[p]) ? +1.0 : -1.0;
            double target = (xp < lower[p]) ? lower[p] : upper[p];
            double delta = (xp < lower[p]) ? (lower[p] - xp) : (xp - upper[p]);

            // STEP 2: BTRAN: rho = B^{-T} e_r
            rho.clear();
            rho.add(r, 1.0);
            factor.btran(rho);
            double exact_wr = rho.norm2_sq();
            if (exact_wr >= 1e-4) dse_weight[r] = exact_wr;

            // STEP 3: PRICE: alpha_j = rho^T a_j for nonbasic j
            price_pivot_row(rho, alpha_row);

            // STEP 4: CHUZC (Harris two-pass ratio test + Bound Flipping)
            struct Cand {
                Idx j;
                double tilde_alpha;
                double ratio;
                double harris_ratio;
            };
            std::vector<Cand> cands;
            double max_abs_alpha = 0.0;
            for (Idx k = 0; k < alpha_row.count; ++k) {
                Idx j = alpha_row.idx[k];
                if (status[j] == VarStatus::Basic || status[j] == VarStatus::Fixed) continue;
                double ta = s * alpha_row.val[j];
                max_abs_alpha = std::max(max_abs_alpha, std::abs(ta));
                if (status[j] == VarStatus::AtLower && ta < -piv_tol) {
                    double dj = std::max(0.0, dual[j]);
                    cands.push_back({j, ta, dj / (-ta), (dj + eps_d) / (-ta)});
                } else if (status[j] == VarStatus::AtUpper && ta > piv_tol) {
                    double dj = std::max(0.0, -dual[j]);
                    cands.push_back({j, ta, dj / ta, (dj + eps_d) / ta});
                } else if (status[j] == VarStatus::FreeZero && std::abs(ta) > piv_tol) {
                    cands.push_back({j, ta, 0.0, 0.0});
                }
            }

            if (cands.empty()) {
                if (factor.num_updates() > 0) {
                    refactor_and_recompute();
                    continue;
                }
                farkas_ray.assign(m, 0.0);
                for (Idx k = 0; k < rho.count; ++k) farkas_ray[rho.idx[k]] = s * rho.val[rho.idx[k]];
                return Status::Infeasible;
            }

            // Harris Pass 1: find tau_max
            double tau_max = kInf;
            for (const auto& c : cands) {
                tau_max = std::min(tau_max, c.harris_ratio);
            }

            // Harris Pass 2: among candidates with ratio <= tau_max, pick largest |tilde_alpha|
            Idx q = kNone;
            double best_piv_abs = -1.0;
            double tau = 0.0;
            for (const auto& c : cands) {
                if (c.ratio <= tau_max + 1e-13) {
                    if (std::abs(c.tilde_alpha) > best_piv_abs) {
                        best_piv_abs = std::abs(c.tilde_alpha);
                        q = c.j;
                        tau = c.ratio;
                    }
                }
            }
            if (q == kNone) {
                q = cands[0].j;
                tau = cands[0].ratio;
            }

            // STEP 5: FTRAN entering column a_q
            extract_col_of_M(q, alpha_col);
            factor.ftran(alpha_col, &spike);
            double alpha_rq = alpha_col.val[r];
            double alpha_rq_row = alpha_row.val[q];

            if (std::abs(alpha_rq) < 1e-9 ||
                std::abs(alpha_rq - alpha_rq_row) > 1e-5 * (1.0 + std::abs(alpha_rq_row))) {
                if (factor.num_updates() > 0) {
                    refactor_and_recompute();
                    continue;
                }
            }

            // STEP 7: DSE weight update (Forrest-Goldfarb)
            tau_vec.copy_from(rho);
            factor.ftran(tau_vec);
            double wr_old = dse_weight[r];
            dse_weight[r] = std::max(1e-4, wr_old / (alpha_rq * alpha_rq));
            for (Idx k = 0; k < alpha_col.count; ++k) {
                Idx i = alpha_col.idx[k];
                if (i == r) continue;
                double ratio_i = alpha_col.val[i] / alpha_rq;
                double new_w = dse_weight[i] - 2.0 * ratio_i * tau_vec.val[i] + ratio_i * ratio_i * wr_old;
                dse_weight[i] = std::max(1e-4, new_w);
            }

            // STEP 8: Primal update
            double theta_P = (xp - target) / alpha_rq;
            for (Idx k = 0; k < alpha_col.count; ++k) {
                Idx i = alpha_col.idx[k];
                value[basic_index[i]] -= theta_P * alpha_col.val[i];
            }
            value[q] += theta_P;
            value[p] = target;

            // STEP 9: Dual update with cost-shift consistency (Part 6.6 STEP 9)
            for (Idx k = 0; k < alpha_row.count; ++k) {
                Idx j = alpha_row.idx[k];
                if (status[j] == VarStatus::Basic) continue;
                dual[j] += tau * (s * alpha_row.val[j]);
                if (status[j] == VarStatus::AtLower && dual[j] < 0.0) {
                    cost[j] -= dual[j];
                    cost_shift[j] -= dual[j];
                    dual[j] = 0.0;
                } else if (status[j] == VarStatus::AtUpper && dual[j] > 0.0) {
                    cost[j] -= dual[j];
                    cost_shift[j] -= dual[j];
                    dual[j] = 0.0;
                }
            }
            dual[q] = 0.0;
            dual[p] = s * tau;

            // STEP 10: Basis change & LU update
            basic_index[r] = q;
            pos_in_basis[q] = r;
            status[q] = VarStatus::Basic;
            pos_in_basis[p] = kNone;
            if (lower[p] == upper[p] && is_finite_bound(lower[p])) {
                status[p] = VarStatus::Fixed;
            } else {
                status[p] = (s > 0.0) ? VarStatus::AtLower : VarStatus::AtUpper;
            }

            auto ustat = factor.update(r, alpha_col, &spike);
            if (ustat == UpdateStatus::Singular) {
                refactor_and_recompute();
            } else {
                for (Idx k = 0; k < alpha_col.count; ++k) {
                    update_primal_infeas_pos(alpha_col.idx[k]);
                }
                update_primal_infeas_pos(r);
                objective += tau * delta;
            }
            ++stats.iterations;
            ++stats.phase2_iterations;
        }
        return Status::IterationLimit;
    }

    // Primal Simplex (Composite Phase 1 + Phase 2) for unbounded detection, Phase 1 fallback, and cleanup (Part 6.10)
    Status run_primal_simplex(const Deadline& deadline, Idx max_iters) {
        double eps_p = opt.tol.primal_feas_tol;
        double eps_d = opt.tol.dual_feas_tol;
        double piv_tol = opt.tol.pivot_abs_tol;

        refactor_and_recompute();
        Idx zero_step_streak = 0;

        for (Idx iter = 0; iter < max_iters; ++iter) {
            if ((iter & 63) == 0 && deadline.expired()) {
                return deadline.is_interrupted() ? Status::Interrupted : Status::TimeLimit;
            }
            if (factor.need_refactor()) {
                refactor_and_recompute();
            }

            // Check whether current basic solution is primal feasible
            double sum_p_infeas = 0.0;
            HVector c1(m);
            for (Idx k = 0; k < m; ++k) {
                Idx var = basic_index[k];
                double x = value[var];
                if (x < lower[var] - eps_p) {
                    sum_p_infeas += (lower[var] - x);
                    c1.add(k, -1.0);
                } else if (x > upper[var] + eps_p) {
                    sum_p_infeas += (x - upper[var]);
                    c1.add(k, +1.0);
                }
            }
            bool in_phase1 = (sum_p_infeas > eps_p);

            // Compute pricing reduced costs (either Phase 1 or Phase 2)
            std::vector<double> d_price(nt, 0.0);
            if (in_phase1) {
                c1.tighten();
                factor.btran(c1);
                std::vector<double> y1(m, 0.0);
                for (Idx k = 0; k < c1.count; ++k) y1[c1.idx[k]] = c1.val[c1.idx[k]];
                for (Idx j = 0; j < n; ++j) {
                    if (status[j] == VarStatus::Basic) continue;
                    double ay = 0.0;
                    for (Off p = A.csc.start[j]; p < A.csc.start[j + 1]; ++p) {
                        ay += A.csc.value[p] * y1[A.csc.index[p]];
                    }
                    d_price[j] = -ay;
                }
                for (Idx i = 0; i < m; ++i) {
                    Idx log_v = n + i;
                    if (status[log_v] == VarStatus::Basic) continue;
                    d_price[log_v] = y1[i];
                }
            } else {
                compute_dual();
                d_price = dual;
            }

            // Choose entering nonbasic variable q (Bland's smallest-index rule when zero_step_streak > 12)
            Idx q = kNone;
            double best_merit = 0.0;
            double sigma = +1.0;
            bool use_bland = (zero_step_streak > 12);
            for (Idx j = 0; j < nt; ++j) {
                if (status[j] == VarStatus::Basic || status[j] == VarStatus::Fixed) continue;
                double dj = d_price[j];
                double viol = 0.0;
                double dir = 0.0;
                if (status[j] == VarStatus::AtLower && dj < -eps_d) {
                    viol = -dj; dir = +1.0;
                } else if (status[j] == VarStatus::AtUpper && dj > eps_d) {
                    viol = dj; dir = -1.0;
                } else if (status[j] == VarStatus::FreeZero && std::abs(dj) > eps_d) {
                    viol = std::abs(dj); dir = (dj < 0.0) ? +1.0 : -1.0;
                }
                if (viol > eps_d) {
                    if (use_bland) {
                        q = j;
                        sigma = dir;
                        break;
                    }
                    double merit = (viol * viol) / std::max(1.0, devex_weight[j]);
                    if (merit > best_merit) {
                        best_merit = merit;
                        q = j;
                        sigma = dir;
                    }
                }
            }

            if (q == kNone) {
                if (in_phase1) {
                    if (factor.num_updates() > 0) {
                        refactor_and_recompute();
                        continue;
                    }
                    return Status::Infeasible;
                }
                refactor_and_recompute();
                return Status::Optimal;
            }

            // FTRAN: alpha_q = B^{-1} a_q; basic variables move as X_B(theta) = X_B - sigma * theta * alpha_q
            extract_col_of_M(q, alpha_col);
            factor.ftran(alpha_col, &spike);

            // Ratio test
            double own_range = (is_finite_bound(lower[q]) && is_finite_bound(upper[q]))
                               ? (upper[q] - lower[q]) : kInf;
            Idx r_leave = kNone;
            double theta_min = own_range;
            double best_piv = 0.0;
            int leave_to_upper = 0;

            for (Idx k = 0; k < alpha_col.count; ++k) {
                Idx pos = alpha_col.idx[k];
                Idx bvar = basic_index[pos];
                double gk = sigma * alpha_col.val[pos]; // X_pos decreases at rate gk
                double xk = value[bvar];

                if (in_phase1) {
                    if (xk >= lower[bvar] - eps_p && xk <= upper[bvar] + eps_p) {
                        // Currently feasible basic variable: must not violate [lower, upper]
                        if (gk > piv_tol && is_finite_bound(lower[bvar])) {
                            double ratio = std::max(0.0, xk - lower[bvar]) / gk;
                            if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                                theta_min = ratio;
                                r_leave = pos;
                                best_piv = std::abs(gk);
                                leave_to_upper = 0;
                            }
                        } else if (gk < -piv_tol && is_finite_bound(upper[bvar])) {
                            double ratio = std::max(0.0, upper[bvar] - xk) / (-gk);
                            if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                                theta_min = ratio;
                                r_leave = pos;
                                best_piv = std::abs(gk);
                                leave_to_upper = 1;
                            }
                        }
                    } else if (xk > upper[bvar] + eps_p && gk > piv_tol) {
                        // Currently above upper bound and decreasing toward upper[bvar]
                        double ratio = (xk - upper[bvar]) / gk;
                        if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                            theta_min = ratio;
                            r_leave = pos;
                            best_piv = std::abs(gk);
                            leave_to_upper = 1;
                        }
                    } else if (xk < lower[bvar] - eps_p && gk < -piv_tol) {
                        // Currently below lower bound and increasing toward lower[bvar]
                        double ratio = (lower[bvar] - xk) / (-gk);
                        if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                            theta_min = ratio;
                            r_leave = pos;
                            best_piv = std::abs(gk);
                            leave_to_upper = 0;
                        }
                    }
                } else {
                    if (gk > piv_tol && is_finite_bound(lower[bvar])) {
                        double ratio = std::max(0.0, xk - lower[bvar]) / gk;
                        if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                            theta_min = ratio;
                            r_leave = pos;
                            best_piv = std::abs(gk);
                            leave_to_upper = 0;
                        }
                    } else if (gk < -piv_tol && is_finite_bound(upper[bvar])) {
                        double ratio = std::max(0.0, upper[bvar] - xk) / (-gk);
                        if (ratio < theta_min - 1e-12 || (std::abs(ratio - theta_min) <= 1e-12 && std::abs(gk) > best_piv)) {
                            theta_min = ratio;
                            r_leave = pos;
                            best_piv = std::abs(gk);
                            leave_to_upper = 1;
                        }
                    }
                }
            }

            if (r_leave == kNone && !is_finite_bound(theta_min)) {
                if (!in_phase1) {
                    return Status::Unbounded;
                }
                return Status::Infeasible;
            }

            // Bound flip of entering variable if it hits its own opposite bound first
            if (r_leave == kNone && is_finite_bound(own_range)) {
                for (Idx k = 0; k < alpha_col.count; ++k) {
                    Idx pos = alpha_col.idx[k];
                    value[basic_index[pos]] -= sigma * own_range * alpha_col.val[pos];
                }
                status[q] = (status[q] == VarStatus::AtLower) ? VarStatus::AtUpper : VarStatus::AtLower;
                value[q] = (status[q] == VarStatus::AtLower) ? lower[q] : upper[q];
                update_primal_infeas_all();
                recompute_objective();
                ++stats.bound_flips;
                ++stats.iterations;
                continue;
            }

            if (theta_min <= 1e-10) ++zero_step_streak;
            else zero_step_streak = 0;

            // Pivot q into position r_leave
            Idx p = basic_index[r_leave];
            for (Idx k = 0; k < alpha_col.count; ++k) {
                Idx pos = alpha_col.idx[k];
                value[basic_index[pos]] -= sigma * theta_min * alpha_col.val[pos];
            }
            value[q] += sigma * theta_min;
            value[p] = leave_to_upper ? upper[p] : lower[p];

            basic_index[r_leave] = q;
            pos_in_basis[q] = r_leave;
            status[q] = VarStatus::Basic;
            pos_in_basis[p] = kNone;
            if (lower[p] == upper[p] && is_finite_bound(lower[p])) {
                status[p] = VarStatus::Fixed;
            } else {
                status[p] = leave_to_upper ? VarStatus::AtUpper : VarStatus::AtLower;
            }

            if (factor.update(r_leave, alpha_col, &spike) == UpdateStatus::Singular) {
                refactor_and_recompute();
            } else {
                update_primal_infeas_all();
                recompute_objective();
            }
            ++stats.iterations;
            if (in_phase1) ++stats.phase1_iterations;
            else ++stats.phase2_iterations;
        }
        return Status::IterationLimit;
    }

    // Main LP solve entry point (Part 6.5 + 6.6 + 6.9)
    Status solve(const Deadline& deadline) {
        cost = cost_orig;
        refactor_and_recompute();

        if (opt.method == "primal") {
            return run_primal_simplex(deadline, opt.simplex_iter_limit);
        }

        // Prepare dual feasibility on current basis (6.4.3 & 6.5)
        Idx num_dual_infeas = make_boxed_dual_feasible();

        if (num_dual_infeas > 0) {
            // Dual Phase 1 via artificial bounding (Method A of 6.5):
            // Keep original costs so we find a dual-feasible basis for the true cost vector,
            // while boxing all variables in artificial intervals so every basis is dual-feasible.
            std::vector<double> saved_lo = lower;
            std::vector<double> saved_up = upper;

            for (Idx j = 0; j < nt; ++j) {
                bool has_lo = is_finite_bound(saved_lo[j]);
                bool has_up = is_finite_bound(saved_up[j]);
                if (!has_lo && !has_up) {
                    lower[j] = -1000.0;
                    upper[j] = +1000.0;
                } else if (has_lo && !has_up) {
                    lower[j] = 0.0;
                    upper[j] = 1.0;
                } else if (!has_lo && has_up) {
                    lower[j] = -1.0;
                    upper[j] = 0.0;
                } else {
                    lower[j] = 0.0;
                    upper[j] = 0.0;
                }
            }
            refactor_and_recompute();
            make_boxed_dual_feasible();
            Status p1_stat = run_dual_phase2(deadline, opt.simplex_iter_limit / 2, false);
            lower = saved_lo;
            upper = saved_up;
            refactor_and_recompute();
            Idx rem_infeas = make_boxed_dual_feasible();

            if (p1_stat != Status::Optimal || rem_infeas > 0) {
                // Fallback to Primal Simplex (Method B of 6.5: handles unbounded/infeasible and general Phase 1)
                return run_primal_simplex(deadline, opt.simplex_iter_limit);
            }
        }

        // Run Dual Simplex Phase 2
        Status st = run_dual_phase2(deadline, opt.simplex_iter_limit, true);
        if (st == Status::Infeasible || st == Status::IterationLimit) {
            // Confirm with Primal Simplex (Part 6.5 & 6.8 safeguard)
            init_slack_basis();
            Status pst = run_primal_simplex(deadline, opt.simplex_iter_limit);
            if (pst == Status::Optimal || pst == Status::Unbounded) {
                st = pst;
            }
        }
        if (st == Status::Optimal) {
            // Clean-up verification (6.9): check dual feasibility with original costs
            cost = cost_orig;
            refactor_and_recompute();
            double max_d_viol = 0.0;
            for (Idx j = 0; j < nt; ++j) {
                if (status[j] == VarStatus::Basic || status[j] == VarStatus::Fixed) continue;
                if (status[j] == VarStatus::AtLower && dual[j] < 0.0) max_d_viol = std::max(max_d_viol, -dual[j]);
                else if (status[j] == VarStatus::AtUpper && dual[j] > 0.0) max_d_viol = std::max(max_d_viol, dual[j]);
                else if (status[j] == VarStatus::FreeZero) max_d_viol = std::max(max_d_viol, std::abs(dual[j]));
            }
            if (max_d_viol > opt.tol.dual_feas_tol * 10.0) {
                st = run_primal_simplex(deadline, opt.simplex_iter_limit / 4);
            }
        }
        return st;
    }

    // Warm-start and MIP modification API (Part 6.14)
    void get_basis(BasisSnapshot& snap) const {
        snap.status = status;
        snap.basic_index = basic_index;
    }

    void set_basis(const BasisSnapshot& snap) {
        if (static_cast<Idx>(snap.status.size()) == nt &&
            static_cast<Idx>(snap.basic_index.size()) == m) {
            status = snap.status;
            basic_index = snap.basic_index;
            std::fill(pos_in_basis.begin(), pos_in_basis.end(), kNone);
            for (Idx k = 0; k < m; ++k) {
                pos_in_basis[basic_index[k]] = k;
            }
            sync_nonbasic_values();
        }
    }

    void change_col_bounds(Idx j, double lo, double up) {
        lower[j] = lo;
        upper[j] = up;
        if (status[j] != VarStatus::Basic) {
            if (lo == up && is_finite_bound(lo)) {
                status[j] = VarStatus::Fixed;
                value[j] = lo;
            } else if (status[j] == VarStatus::AtLower || status[j] == VarStatus::Fixed) {
                status[j] = VarStatus::AtLower;
                value[j] = lo;
            } else if (status[j] == VarStatus::AtUpper) {
                value[j] = up;
            }
        }
    }

    void add_rows(const std::vector<SparseRow>& new_rows) {
        if (new_rows.empty()) return;
        std::vector<Triplet> trips;
        for (Idx j = 0; j < n; ++j) {
            for (Off k = A.csc.start[j]; k < A.csc.start[j + 1]; ++k) {
                trips.push_back({A.csc.index[k], j, A.csc.value[k]});
            }
        }
        Idx old_m = m;
        Idx add_m = static_cast<Idx>(new_rows.size());
        for (Idx r = 0; r < add_m; ++r) {
            Idx row_idx = old_m + r;
            for (size_t k = 0; k < new_rows[r].idx.size(); ++k) {
                trips.push_back({row_idx, new_rows[r].idx[k], new_rows[r].val[k]});
            }
        }
        m = old_m + add_m;
        nt = n + m;
        A.build(CscMatrix::from_triplets(m, n, trips));

        for (Idx r = 0; r < add_m; ++r) {
            cost.push_back(0.0);
            cost_orig.push_back(0.0);
            cost_shift.push_back(0.0);
            lower.push_back(new_rows[r].lo);
            upper.push_back(new_rows[r].up);
            value.push_back(0.0);
            dual.push_back(0.0);
            status.push_back(VarStatus::Basic);
            Idx log_v = n + old_m + r;
            basic_index.push_back(log_v);
            pos_in_basis.push_back(old_m + r);
            primal_infeas.push_back(0.0);
            dse_weight.push_back(1.0);
            devex_weight.push_back(1.0);
        }
        row_duals.resize(m, 0.0);
        rho.resize(m);
        alpha_row.resize(nt);
        alpha_col.resize(m);
        spike.resize(m);
        tau_vec.resize(m);
        flip_rhs.resize(m);
    }

    // Compute tableau row for basis position k: X_{basic_var} + sum_{j nonbasic} alpha_j X_j = 0
    void get_tableau_row(Idx basis_pos, HVector& out_row, Idx& basic_var) const {
        basic_var = basic_index[basis_pos];
        HVector rho_tmp(m);
        rho_tmp.add(basis_pos, 1.0);
        factor.btran(rho_tmp);
        price_pivot_row(rho_tmp, out_row);
    }
};

} // namespace sov
