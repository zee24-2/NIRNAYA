#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../model/model.hpp"

namespace sov {

struct PresolveResult {
    Status status = Status::Optimal; // Optimal means "continue to solve reduced model" or "solved in presolve"
    bool solved_in_presolve = false;
    Model reduced;
    std::vector<Idx> orig_col_map; // reduced col j -> original col idx
    std::vector<Idx> orig_row_map; // reduced row i -> original row idx
    std::vector<double> fixed_x;   // size orig.ncol
    std::vector<bool> is_col_fixed;
    Idx rows_removed = 0;
    Idx cols_removed = 0;
    Idx bound_tightenings = 0;

    void postsolve(const std::vector<double>& x_red, std::vector<double>& x_orig) const {
        x_orig = fixed_x;
        for (size_t j = 0; j < orig_col_map.size() && j < x_red.size(); ++j) {
            x_orig[orig_col_map[j]] = x_red[j];
        }
    }
};

// LP & MIP Presolve and Postsolve Engine (Part 7: R1-R7, R14, R20)
class Presolver {
public:
    static PresolveResult presolve(const Model& orig, int level, int max_rounds = 20) {
        PresolveResult res;
        res.fixed_x.assign(orig.ncol, 0.0);
        res.is_col_fixed.assign(orig.ncol, false);

        if (level <= 0 || orig.is_qp()) {
            res.reduced = orig;
            res.orig_col_map.resize(orig.ncol);
            for (Idx j = 0; j < orig.ncol; ++j) res.orig_col_map[j] = j;
            res.orig_row_map.resize(orig.nrow);
            for (Idx i = 0; i < orig.nrow; ++i) res.orig_row_map[i] = i;
            return res;
        }

        Model work = orig;
        Idx m = work.nrow;
        Idx n = work.ncol;
        CsrMatrix csr = CsrMatrix::from_csc(work.A);

        std::vector<bool> row_del(m, false);
        std::vector<bool> col_del(n, false);

        auto fix_col = [&](Idx j, double val) {
            if (col_del[j]) return;
            col_del[j] = true;
            res.is_col_fixed[j] = true;
            res.fixed_x[j] = val;
            work.obj_offset += work.cost[j] * val;
            for (Off k = work.A.start[j]; k < work.A.start[j + 1]; ++k) {
                Idx i = work.A.index[k];
                if (row_del[i]) continue;
                double aij = work.A.value[k];
                if (is_finite_bound(work.rowlo[i])) work.rowlo[i] -= aij * val;
                if (is_finite_bound(work.rowup[i])) work.rowup[i] -= aij * val;
            }
            ++res.cols_removed;
        };

        for (int round = 0; round < max_rounds; ++round) {
            bool changed = false;

            // R4: Fixed columns
            for (Idx j = 0; j < n; ++j) {
                if (col_del[j]) continue;
                if (is_finite_bound(work.collo[j]) && is_finite_bound(work.colup[j]) &&
                    std::abs(work.colup[j] - work.collo[j]) <= 1e-11) {
                    fix_col(j, work.collo[j]);
                    changed = true;
                }
            }

            // R1 & R3: Empty rows and row singletons
            for (Idx i = 0; i < m; ++i) {
                if (row_del[i]) continue;
                Idx cnt = 0;
                Idx single_j = kNone;
                double single_a = 0.0;
                for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                    Idx j = csr.index[k];
                    if (!col_del[j] && std::abs(csr.value[k]) > 1e-13) {
                        ++cnt;
                        single_j = j;
                        single_a = csr.value[k];
                    }
                }
                if (cnt == 0) {
                    if ((is_finite_bound(work.rowlo[i]) && work.rowlo[i] > 1e-7) ||
                        (is_finite_bound(work.rowup[i]) && work.rowup[i] < -1e-7)) {
                        res.status = Status::Infeasible;
                        res.solved_in_presolve = true;
                        return res;
                    }
                    row_del[i] = true;
                    ++res.rows_removed;
                    changed = true;
                } else if (cnt == 1) {
                    double lo_b = -kInf, up_b = kInf;
                    if (single_a > 0.0) {
                        if (is_finite_bound(work.rowlo[i])) lo_b = work.rowlo[i] / single_a;
                        if (is_finite_bound(work.rowup[i])) up_b = work.rowup[i] / single_a;
                    } else {
                        if (is_finite_bound(work.rowup[i])) lo_b = work.rowup[i] / single_a;
                        if (is_finite_bound(work.rowlo[i])) up_b = work.rowlo[i] / single_a;
                    }
                    if (work.vartype[single_j] != VarType::Continuous) {
                        if (is_finite_bound(lo_b)) lo_b = std::ceil(lo_b - 1e-7);
                        if (is_finite_bound(up_b)) up_b = std::floor(up_b + 1e-7);
                    }
                    if (lo_b > work.collo[single_j] + 1e-9) {
                        work.collo[single_j] = lo_b;
                        ++res.bound_tightenings;
                    }
                    if (up_b < work.colup[single_j] - 1e-9) {
                        work.colup[single_j] = up_b;
                        ++res.bound_tightenings;
                    }
                    if (work.collo[single_j] > work.colup[single_j] + 1e-7) {
                        res.status = Status::Infeasible;
                        res.solved_in_presolve = true;
                        return res;
                    }
                    if (work.collo[single_j] > work.colup[single_j]) {
                        work.collo[single_j] = work.colup[single_j];
                    }
                    row_del[i] = true;
                    ++res.rows_removed;
                    changed = true;
                }
            }

            // R5, R6, R7, R14, R20: Activity bounds, redundancy, forcing rows, bound & coeff tightening
            for (Idx i = 0; i < m; ++i) {
                if (row_del[i]) continue;
                double min_act = 0.0, max_act = 0.0;
                int min_inf = 0, max_inf = 0;

                for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                    Idx j = csr.index[k];
                    if (col_del[j]) continue;
                    double a = csr.value[k];
                    if (a > 0.0) {
                        if (is_finite_bound(work.collo[j])) min_act += a * work.collo[j];
                        else ++min_inf;
                        if (is_finite_bound(work.colup[j])) max_act += a * work.colup[j];
                        else ++max_inf;
                    } else if (a < 0.0) {
                        if (is_finite_bound(work.colup[j])) min_act += a * work.colup[j];
                        else ++min_inf;
                        if (is_finite_bound(work.collo[j])) max_act += a * work.collo[j];
                        else ++max_inf;
                    }
                }

                // R6: Activity infeasibility
                if ((min_inf == 0 && is_finite_bound(work.rowup[i]) && min_act > work.rowup[i] + 1e-7) ||
                    (max_inf == 0 && is_finite_bound(work.rowlo[i]) && max_act < work.rowlo[i] - 1e-7)) {
                    res.status = Status::Infeasible;
                    res.solved_in_presolve = true;
                    return res;
                }

                // R5: Redundant row
                bool lo_red = (!is_finite_bound(work.rowlo[i])) || (min_inf == 0 && min_act >= work.rowlo[i] - 1e-9);
                bool up_red = (!is_finite_bound(work.rowup[i])) || (max_inf == 0 && max_act <= work.rowup[i] + 1e-9);
                if (lo_red && up_red) {
                    row_del[i] = true;
                    ++res.rows_removed;
                    changed = true;
                    continue;
                }

                // R14: Forcing row
                if (min_inf == 0 && is_finite_bound(work.rowup[i]) && std::abs(min_act - work.rowup[i]) <= 1e-9) {
                    for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                        Idx j = csr.index[k];
                        if (col_del[j]) continue;
                        fix_col(j, csr.value[k] > 0.0 ? work.collo[j] : work.colup[j]);
                    }
                    row_del[i] = true;
                    ++res.rows_removed;
                    changed = true;
                    continue;
                }
                if (max_inf == 0 && is_finite_bound(work.rowlo[i]) && std::abs(max_act - work.rowlo[i]) <= 1e-9) {
                    for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                        Idx j = csr.index[k];
                        if (col_del[j]) continue;
                        fix_col(j, csr.value[k] > 0.0 ? work.colup[j] : work.collo[j]);
                    }
                    row_del[i] = true;
                    ++res.rows_removed;
                    changed = true;
                    continue;
                }

                // R7: Bound tightening on integer variables (and tight bounds on continuous in level 2)
                if (min_inf == 0 && is_finite_bound(work.rowup[i])) {
                    double slack = work.rowup[i] - min_act;
                    for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                        Idx j = csr.index[k];
                        if (col_del[j] || work.vartype[j] == VarType::Continuous) continue;
                        double a = csr.value[k];
                        if (a > 1e-9 && is_finite_bound(work.collo[j])) {
                            double new_u = std::floor(work.collo[j] + (slack / a) + 1e-6);
                            if (new_u < work.colup[j] - 0.5) {
                                work.colup[j] = new_u;
                                ++res.bound_tightenings;
                                changed = true;
                            }
                        } else if (a < -1e-9 && is_finite_bound(work.colup[j])) {
                            double new_l = std::ceil(work.colup[j] + (slack / a) - 1e-6);
                            if (new_l > work.collo[j] + 0.5) {
                                work.collo[j] = new_l;
                                ++res.bound_tightenings;
                                changed = true;
                            }
                        }
                    }
                }
                if (max_inf == 0 && is_finite_bound(work.rowlo[i])) {
                    double slack = max_act - work.rowlo[i];
                    for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                        Idx j = csr.index[k];
                        if (col_del[j] || work.vartype[j] == VarType::Continuous) continue;
                        double a = csr.value[k];
                        if (a > 1e-9 && is_finite_bound(work.colup[j])) {
                            double new_l = std::ceil(work.colup[j] - (slack / a) - 1e-6);
                            if (new_l > work.collo[j] + 0.5) {
                                work.collo[j] = new_l;
                                ++res.bound_tightenings;
                                changed = true;
                            }
                        } else if (a < -1e-9 && is_finite_bound(work.collo[j])) {
                            double new_u = std::floor(work.collo[j] - (slack / a) + 1e-6);
                            if (new_u < work.colup[j] - 0.5) {
                                work.colup[j] = new_u;
                                ++res.bound_tightenings;
                                changed = true;
                            }
                        }
                    }
                }
            }

            // R2: Empty columns
            for (Idx j = 0; j < n; ++j) {
                if (col_del[j]) continue;
                bool has_active_row = false;
                for (Off k = work.A.start[j]; k < work.A.start[j + 1]; ++k) {
                    if (!row_del[work.A.index[k]]) {
                        has_active_row = true;
                        break;
                    }
                }
                if (!has_active_row) {
                    double cj = work.cost[j];
                    if (cj > 1e-12) {
                        if (!is_finite_bound(work.collo[j])) {
                            res.status = Status::Unbounded;
                            res.solved_in_presolve = true;
                            return res;
                        }
                        fix_col(j, work.collo[j]);
                    } else if (cj < -1e-12) {
                        if (!is_finite_bound(work.colup[j])) {
                            res.status = Status::Unbounded;
                            res.solved_in_presolve = true;
                            return res;
                        }
                        fix_col(j, work.colup[j]);
                    } else {
                        double v = 0.0;
                        if (is_finite_bound(work.collo[j]) && v < work.collo[j]) v = work.collo[j];
                        if (is_finite_bound(work.colup[j]) && v > work.colup[j]) v = work.colup[j];
                        fix_col(j, v);
                    }
                    changed = true;
                }
            }

            if (!changed) break;
        }

        // Finalize reduced model
        Model red;
        red.name = orig.name + "_pre";
        red.sense = orig.sense;
        red.obj_offset = work.obj_offset;

        std::vector<Idx> new_col_idx(n, kNone);
        std::vector<Idx> new_row_idx(m, kNone);

        for (Idx j = 0; j < n; ++j) {
            if (!col_del[j]) {
                Idx nj = red.ncol++;
                new_col_idx[j] = nj;
                res.orig_col_map.push_back(j);
                red.cost.push_back(work.cost[j]);
                red.collo.push_back(work.collo[j]);
                red.colup.push_back(work.colup[j]);
                red.vartype.push_back(work.vartype[j]);
                red.col_names.push_back(work.col_names[j]);
            }
        }
        for (Idx i = 0; i < m; ++i) {
            if (!row_del[i]) {
                Idx ni = red.nrow++;
                new_row_idx[i] = ni;
                res.orig_row_map.push_back(i);
                red.rowlo.push_back(work.rowlo[i]);
                red.rowup.push_back(work.rowup[i]);
                red.row_names.push_back(work.row_names[i]);
            }
        }

        std::vector<Triplet> red_trips;
        for (Idx j = 0; j < n; ++j) {
            if (col_del[j]) continue;
            Idx nj = new_col_idx[j];
            for (Off k = work.A.start[j]; k < work.A.start[j + 1]; ++k) {
                Idx i = work.A.index[k];
                if (row_del[i]) continue;
                red_trips.push_back({new_row_idx[i], nj, work.A.value[k]});
            }
        }
        red.A = CscMatrix::from_triplets(red.nrow, red.ncol, red_trips);
        red.Q = CscMatrix(red.ncol, red.ncol);
        res.reduced = std::move(red);
        if (res.reduced.ncol == 0 && res.reduced.nrow == 0) {
            res.solved_in_presolve = true;
            res.status = Status::Optimal;
        }
        return res;
    }
};

} // namespace sov
