#pragma once
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/status.hpp"
#include "../core/options.hpp"
#include "../core/timer.hpp"
#include "../core/log.hpp"
#include "../model/model.hpp"
#include "../io/sol_io.hpp"
#include "../lp/lp_solver.hpp"

namespace sov {

struct MipResult {
    Status status = Status::Infeasible;
    bool has_incumbent = false;
    double incumbent_obj = kInf;
    double best_bound = -kInf;
    double gap = kInf;
    std::vector<double> incumbent_x;
    Idx nodes_explored = 0;
    Idx lp_iterations = 0;
    Idx cuts_added = 0;
    Idx heur_successes = 0;
};

class MipSolver {
public:
    struct Node {
        Idx id = 0;
        Idx depth = 0;
        double lower_bound = -kInf;
        double estimate = -kInf;
        std::vector<double> col_lo;
        std::vector<double> col_up;
        BasisSnapshot basis;
    };

    struct Pseudocost {
        double sum_down = 0.0;
        Idx count_down = 0;
        double sum_up = 0.0;
        Idx count_up = 0;
    };

    static MipResult solve(
        const Model& model,
        const Options& opt,
        const Deadline& deadline,
        const Logger& logger)
    {
        MipResult res;
        Idx n = model.ncol;
        Idx m = model.nrow;
        CsrMatrix csr = CsrMatrix::from_csc(model.A);

        // Identify integer variables and objective integrality scale (R25)
        std::vector<Idx> int_vars;
        bool all_costs_int = true;
        for (Idx j = 0; j < n; ++j) {
            if (model.vartype[j] == VarType::Integer || model.vartype[j] == VarType::Binary) {
                int_vars.push_back(j);
            }
            if (std::abs(model.cost[j]) > 1e-12) {
                if (model.vartype[j] == VarType::Continuous ||
                    std::abs(model.cost[j] - std::round(model.cost[j])) > 1e-7) {
                    all_costs_int = false;
                }
            }
        }

        auto strengthen_bound = [&](double b) -> double {
            if (all_costs_int && std::isfinite(b)) {
                return std::ceil(b - 1e-6);
            }
            return b;
        };

        auto is_int_feasible = [&](const std::vector<double>& x) -> bool {
            for (Idx j : int_vars) {
                if (std::abs(x[j] - std::round(x[j])) > opt.tol.int_tol) return false;
            }
            return true;
        };

        // Polish an integer assignment by fixing integers and solving for continuous variables
        auto polish_and_try_incumbent = [&](const std::vector<double>& cand_x, const char* tag) -> bool {
            std::vector<double> x_fixed(n, 0.0);
            for (Idx j = 0; j < n; ++j) {
                if (model.vartype[j] != VarType::Continuous) {
                    x_fixed[j] = std::round(cand_x[j]);
                    if (x_fixed[j] < model.collo[j] || x_fixed[j] > model.colup[j]) return false;
                } else {
                    x_fixed[j] = cand_x[j];
                }
            }

            if (static_cast<Idx>(int_vars.size()) < n) {
                // Solve continuous LP with integers fixed
                LpSolver sub_lp(opt);
                sub_lp.load(model);
                for (Idx j : int_vars) {
                    sub_lp.change_col_bounds(j, x_fixed[j], x_fixed[j]);
                }
                Status st = sub_lp.solve(deadline);
                if (st != Status::Optimal) return false;
                for (Idx j = 0; j < n; ++j) {
                    if (model.vartype[j] == VarType::Continuous) {
                        x_fixed[j] = sub_lp.value[j];
                    }
                }
            }

            double obj = model.eval_objective(x_fixed);
            auto rep = verify_solution(model, x_fixed, obj, 1e-5);
            if (!rep.passed) return false;

            if (!res.has_incumbent || obj < res.incumbent_obj - 1e-7 * (1.0 + std::abs(res.incumbent_obj))) {
                res.has_incumbent = true;
                res.incumbent_obj = obj;
                res.incumbent_x = x_fixed;
                if (tag != nullptr) ++res.heur_successes;
                logger.log(2, "  *%-3s New incumbent: %.10g (time %.2fs)\n",
                           tag ? tag : "B",
                           model.user_objective_from_internal(obj),
                           deadline.timer.elapsed_sec());
                return true;
            }
            return false;
        };

        // Domain propagation on [lo, up] using row activity bounds (10.4.1)
        auto propagate_bounds = [&](std::vector<double>& lo, std::vector<double>& up) -> bool {
            for (int pass = 0; pass < 4; ++pass) {
                bool changed = false;
                for (Idx i = 0; i < m; ++i) {
                    double min_act = 0.0, max_act = 0.0;
                    int min_inf = 0, max_inf = 0;
                    for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                        Idx j = csr.index[k];
                        double a = csr.value[k];
                        if (a > 0.0) {
                            if (is_finite_bound(lo[j])) min_act += a * lo[j]; else ++min_inf;
                            if (is_finite_bound(up[j])) max_act += a * up[j]; else ++max_inf;
                        } else if (a < 0.0) {
                            if (is_finite_bound(up[j])) min_act += a * up[j]; else ++min_inf;
                            if (is_finite_bound(lo[j])) max_act += a * lo[j]; else ++max_inf;
                        }
                    }
                    if ((min_inf == 0 && is_finite_bound(model.rowup[i]) && min_act > model.rowup[i] + 1e-6) ||
                        (max_inf == 0 && is_finite_bound(model.rowlo[i]) && max_act < model.rowlo[i] - 1e-6)) {
                        return false;
                    }
                    if (min_inf == 0 && is_finite_bound(model.rowup[i])) {
                        double slack = model.rowup[i] - min_act;
                        for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                            Idx j = csr.index[k];
                            if (model.vartype[j] == VarType::Continuous) continue;
                            double a = csr.value[k];
                            if (a > 1e-9 && is_finite_bound(lo[j])) {
                                double new_u = std::floor(lo[j] + slack / a + 1e-6);
                                if (new_u < up[j] - 0.5) {
                                    up[j] = new_u;
                                    if (lo[j] > up[j] + 1e-6) return false;
                                    changed = true;
                                }
                            } else if (a < -1e-9 && is_finite_bound(up[j])) {
                                double new_l = std::ceil(up[j] + slack / a - 1e-6);
                                if (new_l > lo[j] + 0.5) {
                                    lo[j] = new_l;
                                    if (lo[j] > up[j] + 1e-6) return false;
                                    changed = true;
                                }
                            }
                        }
                    }
                    if (max_inf == 0 && is_finite_bound(model.rowlo[i])) {
                        double slack = max_act - model.rowlo[i];
                        for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                            Idx j = csr.index[k];
                            if (model.vartype[j] == VarType::Continuous) continue;
                            double a = csr.value[k];
                            if (a > 1e-9 && is_finite_bound(up[j])) {
                                double new_l = std::ceil(up[j] - slack / a - 1e-6);
                                if (new_l > lo[j] + 0.5) {
                                    lo[j] = new_l;
                                    if (lo[j] > up[j] + 1e-6) return false;
                                    changed = true;
                                }
                            } else if (a < -1e-9 && is_finite_bound(lo[j])) {
                                double new_u = std::floor(lo[j] - slack / a + 1e-6);
                                if (new_u < up[j] - 0.5) {
                                    up[j] = new_u;
                                    if (lo[j] > up[j] + 1e-6) return false;
                                    changed = true;
                                }
                            }
                        }
                    }
                }
                if (!changed) break;
            }
            return true;
        };

        // Reduced-cost fixing (10.4.2)
        auto reduced_cost_fix = [&](const LpSolver& lp, double z_lp, std::vector<double>& lo, std::vector<double>& up) {
            if (!res.has_incumbent) return;
            double z_cut = res.incumbent_obj - (all_costs_int ? (1.0 - 1e-5) : 1e-6);
            if (z_lp >= z_cut) return;
            double gap_room = z_cut - z_lp;
            for (Idx j : int_vars) {
                if (lp.status[j] == VarStatus::AtLower && lp.dual[j] > 1e-5 && is_finite_bound(lo[j])) {
                    double max_val = std::floor(lo[j] + gap_room / lp.dual[j] + 1e-6);
                    if (max_val < up[j]) up[j] = std::max(lo[j], max_val);
                } else if (lp.status[j] == VarStatus::AtUpper && lp.dual[j] < -1e-5 && is_finite_bound(up[j])) {
                    double min_val = std::ceil(up[j] + gap_room / lp.dual[j] - 1e-6);
                    if (min_val > lo[j]) lo[j] = std::min(up[j], min_val);
                }
            }
        };

        // Primal heuristics: Rounding + Diving + 1-opt (10.8)
        auto run_heuristics = [&](const std::vector<double>& x_lp,
                                  const std::vector<double>& lo,
                                  const std::vector<double>& up) {
            if (!opt.enable_heuristics) return;

            // 1. Direct rounding + propagation
            std::vector<double> x_round(x_lp.begin(), x_lp.begin() + n);
            for (Idx j : int_vars) {
                x_round[j] = std::clamp(std::round(x_lp[j]), lo[j], up[j]);
            }
            polish_and_try_incumbent(x_round, "R");

            // 2. Fractional Diving Heuristic (10.8.2)
            std::vector<double> d_lo = lo, d_up = up;
            std::vector<double> cur_x(x_lp.begin(), x_lp.begin() + n);
            LpSolver dive_lp(opt);
            dive_lp.load(model);
            for (Idx j = 0; j < n; ++j) dive_lp.change_col_bounds(j, d_lo[j], d_up[j]);

            for (int step = 0; step < std::min<int>(25, static_cast<int>(int_vars.size())); ++step) {
                if (deadline.expired()) break;
                Idx best_j = kNone;
                double min_frac = 1.0;
                for (Idx j : int_vars) {
                    if (d_lo[j] == d_up[j]) continue;
                    double f = std::abs(cur_x[j] - std::round(cur_x[j]));
                    if (f > 1e-6 && f < min_frac) {
                        min_frac = f;
                        best_j = j;
                    }
                }
                if (best_j == kNone) {
                    polish_and_try_incumbent(cur_x, "D");
                    break;
                }
                double fix_v = std::clamp(std::round(cur_x[best_j]), d_lo[best_j], d_up[best_j]);
                d_lo[best_j] = d_up[best_j] = fix_v;
                if (!propagate_bounds(d_lo, d_up)) break;
                for (Idx j = 0; j < n; ++j) dive_lp.change_col_bounds(j, d_lo[j], d_up[j]);
                Status st = dive_lp.solve(deadline);
                res.lp_iterations += dive_lp.stats.iterations;
                if (st != Status::Optimal) break;
                cur_x.assign(dive_lp.value.begin(), dive_lp.value.begin() + n);
                if (is_int_feasible(cur_x)) {
                    polish_and_try_incumbent(cur_x, "D");
                    break;
                }
            }
        };

        // Cut Separation: GMI + Knapsack Cover + MIR (Part 10.7)
        auto separate_cuts = [&](const LpSolver& lp, std::vector<SparseRow>& out_cuts) {
            out_cuts.clear();
            if (!opt.enable_cuts) return;

            // 1. Gomory Mixed-Integer (GMI) Cuts from tableau rows of fractional integer basic variables (10.7.1)
            HVector tab_row(lp.nt);
            for (Idx k = 0; k < lp.m && static_cast<int>(out_cuts.size()) < 15; ++k) {
                Idx bvar = lp.basic_index[k];
                if (bvar >= n) continue;
                if (model.vartype[bvar] == VarType::Continuous) continue;
                double x_val = lp.value[bvar];
                double f0 = x_val - std::floor(x_val);
                if (f0 < 0.02 || f0 > 0.98) continue;

                Idx check_bvar = kNone;
                lp.get_tableau_row(k, tab_row, check_bvar);

                // In LpSolver: X_bvar + sum_{j nonbasic} tab_row[j] * X_j = 0
                // With nonbasic substitution:
                //   AtLower: X_j = l_j + x'_j  => a'_j = tab_row[j]
                //   AtUpper: X_j = u_j - x'_j  => a'_j = -tab_row[j]
                // Then X_bvar + sum_{j nonbasic} a'_j x'_j = x_val
                // GMI cut: sum_j g_j x'_j >= 1
                // Converting x'_j back to structural X (and substituting logicals r_i = a_i^T x):
                std::vector<double> cut_struct(n, 0.0);
                double rhs_ge = 1.0;
                bool valid_cut = true;

                for (Idx idx_t = 0; idx_t < tab_row.count; ++idx_t) {
                    Idx j = tab_row.idx[idx_t];
                    if (lp.status[j] == VarStatus::Basic || lp.status[j] == VarStatus::Fixed) continue;
                    double abar = tab_row.val[j];
                    if (std::abs(abar) <= 1e-12) continue;

                    bool at_lo = (lp.status[j] == VarStatus::AtLower);
                    bool at_up = (lp.status[j] == VarStatus::AtUpper);
                    if (!at_lo && !at_up) {
                        valid_cut = false;
                        break;
                    }
                    double aprime = at_lo ? abar : -abar;
                    bool j_is_int = (j < n) && (model.vartype[j] != VarType::Continuous);
                    double gj = 0.0;
                    if (j_is_int) {
                        double fj = aprime - std::floor(aprime);
                        if (fj <= f0) gj = fj / f0;
                        else gj = (1.0 - fj) / (1.0 - f0);
                    } else {
                        if (aprime > 0.0) gj = aprime / f0;
                        else gj = -aprime / (1.0 - f0);
                    }
                    if (gj <= 1e-12) continue;

                    // gj * x'_j:
                    // if at_lo: gj * (X_j - l_j) => +gj * X_j, rhs_ge += gj * l_j
                    // if at_up: gj * (u_j - X_j) => -gj * X_j, rhs_ge -= gj * u_j
                    double coeff_X = at_lo ? gj : -gj;
                    double bnd = at_lo ? lp.lower[j] : lp.upper[j];
                    if (!is_finite_bound(bnd)) {
                        valid_cut = false;
                        break;
                    }
                    rhs_ge += coeff_X * bnd;

                    if (j < n) {
                        cut_struct[j] += coeff_X;
                    } else {
                        // Logical j = n + i represents r_i = sum_c A(i, c) x_c (only valid for original rows i < m)
                        Idx row_i = j - n;
                        if (row_i >= m) {
                            valid_cut = false;
                            break;
                        }
                        for (Off p = csr.start[row_i]; p < csr.start[row_i + 1]; ++p) {
                            cut_struct[csr.index[p]] += coeff_X * csr.value[p];
                        }
                    }
                }

                if (!valid_cut) continue;

                // Evaluate violation and dynamism
                double lhs_at_lp = 0.0, norm_sq = 0.0, max_c = 0.0, min_c = kInf;
                SparseRow sr;
                sr.lo = rhs_ge;
                sr.up = kInf;
                for (Idx j = 0; j < n; ++j) {
                    double v = cut_struct[j];
                    if (std::abs(v) > 1e-9) {
                        sr.idx.push_back(j);
                        sr.val.push_back(v);
                        lhs_at_lp += v * lp.value[j];
                        norm_sq += v * v;
                        max_c = std::max(max_c, std::abs(v));
                        min_c = std::min(min_c, std::abs(v));
                    }
                }
                if (sr.idx.empty() || norm_sq <= 1e-12) continue;
                double viol = rhs_ge - lhs_at_lp;
                double efficacy = viol / std::sqrt(norm_sq);
                if (efficacy >= opt.tol.cut_min_efficacy_root && (max_c / min_c) <= opt.tol.cut_max_dynamism) {
                    out_cuts.push_back(std::move(sr));
                }
            }

            // 2. Knapsack Minimal Cover Cuts on binary <= rows (10.7.3)
            for (Idx i = 0; i < m && static_cast<int>(out_cuts.size()) < 25; ++i) {
                if (!is_finite_bound(model.rowup[i]) || model.rowup[i] <= 0.0) continue;
                bool all_pos_bin = true;
                struct Item { Idx j; double a; double val; };
                std::vector<Item> items;
                for (Off k = csr.start[i]; k < csr.start[i + 1]; ++k) {
                    Idx j = csr.index[k];
                    double a = csr.value[k];
                    if (model.vartype[j] != VarType::Binary || a <= 0.0) {
                        all_pos_bin = false;
                        break;
                    }
                    items.push_back({j, a, lp.value[j]});
                }
                if (!all_pos_bin || items.size() < 2) continue;
                std::sort(items.begin(), items.end(), [](const Item& x, const Item& y) {
                    return (1.0 - x.val) / x.a < (1.0 - y.val) / y.a;
                });
                double sum_a = 0.0;
                double sum_x = 0.0;
                std::vector<Idx> cover;
                double max_cover_a = 0.0;
                for (const auto& it : items) {
                    cover.push_back(it.j);
                    sum_a += it.a;
                    sum_x += it.val;
                    max_cover_a = std::max(max_cover_a, it.a);
                    if (sum_a > model.rowup[i] + 1e-6) break;
                }
                if (sum_a > model.rowup[i] + 1e-6) {
                    // Prune redundant items to obtain a minimal cover (10.7.3)
                    for (size_t idx_c = 0; idx_c < cover.size(); ) {
                        double a_item = 0.0, v_item = 0.0;
                        for (const auto& it : items) {
                            if (it.j == cover[idx_c]) { a_item = it.a; v_item = it.val; break; }
                        }
                        if (sum_a - a_item > model.rowup[i] + 1e-6) {
                            sum_a -= a_item;
                            sum_x -= v_item;
                            cover.erase(cover.begin() + static_cast<ptrdiff_t>(idx_c));
                        } else {
                            ++idx_c;
                        }
                    }
                    max_cover_a = 0.0;
                    for (Idx cj : cover) {
                        for (const auto& it : items) if (it.j == cj) max_cover_a = std::max(max_cover_a, it.a);
                    }
                    double rhs_cov = static_cast<double>(cover.size()) - 1.0;
                    // Extended cover lifting: include any other item with a_j >= max_cover_a
                    for (const auto& it : items) {
                        if (std::find(cover.begin(), cover.end(), it.j) == cover.end() && it.a >= max_cover_a - 1e-9) {
                            cover.push_back(it.j);
                            sum_x += it.val;
                        }
                    }
                    if (sum_x > rhs_cov + 1e-3) {
                        SparseRow sr;
                        sr.lo = -kInf;
                        sr.up = rhs_cov;
                        for (Idx cj : cover) {
                            sr.idx.push_back(cj);
                            sr.val.push_back(1.0);
                        }
                        out_cuts.push_back(std::move(sr));
                    }
                }
            }
        };

        // Solve Root LP
        LpSolver root_lp(opt);
        root_lp.load(model);
        std::vector<double> root_lo = model.collo;
        std::vector<double> root_up = model.colup;
        if (!propagate_bounds(root_lo, root_up)) {
            res.status = Status::Infeasible;
            return res;
        }
        for (Idx j = 0; j < n; ++j) root_lp.change_col_bounds(j, root_lo[j], root_up[j]);

        Status root_st = root_lp.solve(deadline);
        res.lp_iterations += root_lp.stats.iterations;
        if (root_st == Status::Infeasible) {
            res.status = Status::Infeasible;
            return res;
        }
        if (root_st == Status::Unbounded) {
            res.status = Status::InfeasibleOrUnbounded;
            return res;
        }
        if (root_st != Status::Optimal) {
            res.status = root_st;
            return res;
        }

        res.best_bound = strengthen_bound(root_lp.objective);
        std::vector<double> root_x(root_lp.value.begin(), root_lp.value.begin() + n);
        if (is_int_feasible(root_x)) {
            polish_and_try_incumbent(root_x, "LP");
        }
        run_heuristics(root_x, root_lo, root_up);

        // Root Cutting-Plane Loop (10.3 step 3)
        for (int round = 0; round < 8 && opt.enable_cuts; ++round) {
            if (deadline.expired()) break;
            if (res.has_incumbent) {
                double g = (res.incumbent_obj - res.best_bound) / std::max(1.0, std::abs(res.incumbent_obj));
                if (g <= opt.tol.mip_gap || (res.incumbent_obj - res.best_bound) <= opt.tol.mip_abs_gap) {
                    res.gap = 0.0;
                    res.status = Status::Optimal;
                    return res;
                }
            }

            std::vector<SparseRow> new_cuts;
            separate_cuts(root_lp, new_cuts);
            if (new_cuts.empty()) break;

            double prev_obj = root_lp.objective;
            root_lp.add_rows(new_cuts);
            res.cuts_added += static_cast<Idx>(new_cuts.size());
            Status st = root_lp.solve(deadline);
            res.lp_iterations += root_lp.stats.iterations;
            if (st != Status::Optimal) break;

            res.best_bound = std::max(res.best_bound, strengthen_bound(root_lp.objective));
            root_x.assign(root_lp.value.begin(), root_lp.value.begin() + n);
            if (is_int_feasible(root_x)) {
                polish_and_try_incumbent(root_x, "C");
            }
            reduced_cost_fix(root_lp, root_lp.objective, root_lo, root_up);
            run_heuristics(root_x, root_lo, root_up);

            if (root_lp.objective - prev_obj < 1e-4 * (1.0 + std::abs(prev_obj))) {
                break;
            }
        }

        // Check if root solved the MIP
        if (res.has_incumbent) {
            double abs_g = std::max(0.0, res.incumbent_obj - res.best_bound);
            double rel_g = abs_g / std::max(1e-10, std::abs(res.incumbent_obj));
            if (rel_g <= opt.tol.mip_gap || abs_g <= opt.tol.mip_abs_gap) {
                res.gap = rel_g;
                res.status = Status::Optimal;
                return res;
            }
        }

        // Branch-and-Bound Tree Search (Part 10.5)
        std::vector<Pseudocost> pcost(n);
        double avg_ps_down = 1.0, avg_ps_up = 1.0;

        auto cmp_node = [](const Node& a, const Node& b) {
            return a.estimate > b.estimate; // Min-heap by estimate / bound
        };
        std::priority_queue<Node, std::vector<Node>, decltype(cmp_node)> pq(cmp_node);

        Node root_node;
        root_node.id = 0;
        root_node.depth = 0;
        root_node.lower_bound = res.best_bound;
        root_node.estimate = res.best_bound;
        root_node.col_lo = root_lo;
        root_node.col_up = root_up;
        root_lp.get_basis(root_node.basis);
        pq.push(std::move(root_node));

        Idx next_node_id = 1;

        while (!pq.empty()) {
            if (deadline.expired()) {
                res.status = deadline.is_interrupted() ? Status::Interrupted : Status::TimeLimit;
                break;
            }
            if (res.nodes_explored >= opt.node_limit) {
                res.status = Status::NodeLimit;
                break;
            }
            if (opt.solution_limit > 0 && res.has_incumbent && res.heur_successes >= opt.solution_limit) {
                res.status = Status::SolutionLimit;
                break;
            }

            Node cur = pq.top();
            pq.pop();

            // Update global lower bound
            double min_open_bound = cur.lower_bound;
            res.best_bound = std::max(res.best_bound, min_open_bound);

            double z_cut = res.has_incumbent
                ? (res.incumbent_obj - (all_costs_int ? (1.0 - 1e-5) : opt.tol.mip_abs_gap))
                : kInf;
            if (cur.lower_bound >= z_cut) continue;

            ++res.nodes_explored;

            // Propagate node bounds
            if (!propagate_bounds(cur.col_lo, cur.col_up)) continue;

            // Set up and warm-start node LP: set_basis FIRST, then change_col_bounds so Fixed bounds are preserved!
            root_lp.set_basis(cur.basis);
            for (Idx j = 0; j < n; ++j) {
                root_lp.change_col_bounds(j, cur.col_lo[j], cur.col_up[j]);
            }
            root_lp.objective_cutoff = z_cut;

            Status lps = root_lp.solve(deadline);
            res.lp_iterations += root_lp.stats.iterations;
            if (lps != Status::Optimal) continue;

            double z_node = strengthen_bound(root_lp.objective);
            if (z_node >= z_cut) continue;

            std::vector<double> x_node(root_lp.value.begin(), root_lp.value.begin() + n);
            if (is_int_feasible(x_node)) {
                polish_and_try_incumbent(x_node, "B");
                continue;
            }

            reduced_cost_fix(root_lp, root_lp.objective, cur.col_lo, cur.col_up);
            if ((res.nodes_explored % 10) == 0) {
                run_heuristics(x_node, cur.col_lo, cur.col_up);
            }

            // Choose branching variable via Reliability / Strong / Pseudocost Branching (10.6)
            struct BranchCand {
                Idx j;
                double val;
                double f_down;
                double f_up;
                double prelim_score;
            };
            std::vector<BranchCand> cands;
            for (Idx j : int_vars) {
                double v = x_node[j];
                double fd = v - std::floor(v);
                double fu = std::ceil(v) - v;
                if (std::min(fd, fu) > opt.tol.int_tol && cur.col_lo[j] < cur.col_up[j]) {
                    double psd = (pcost[j].count_down > 0) ? (pcost[j].sum_down / pcost[j].count_down) : avg_ps_down;
                    double psu = (pcost[j].count_up > 0) ? (pcost[j].sum_up / pcost[j].count_up) : avg_ps_up;
                    double sc = std::max(1e-6, psd * fd) * std::max(1e-6, psu * fu);
                    cands.push_back({j, v, fd, fu, sc});
                }
            }
            if (cands.empty()) {
                polish_and_try_incumbent(x_node, "B");
                continue;
            }

            std::sort(cands.begin(), cands.end(), [](const BranchCand& a, const BranchCand& b) {
                return a.prelim_score > b.prelim_score;
            });

            Idx best_branch_var = cands[0].j;
            double best_branch_val = cands[0].val;
            double best_score = -1.0;

            BasisSnapshot node_basis;
            root_lp.get_basis(node_basis);

            int sb_evals = 0;
            bool bound_tightened_by_sb = false;
            for (const auto& c : cands) {
                Idx j = c.j;
                double d_down = 0.0, d_up = 0.0;
                bool unreliable = (std::min(pcost[j].count_down, pcost[j].count_up) < opt.strong_branch_rel);
                if ((opt.branching == "reliability" || opt.branching == "strong") && unreliable && sb_evals < 6 && cur.depth <= 8) {
                    ++sb_evals;
                    // Down trial: x_j <= floor(v)
                    root_lp.set_basis(node_basis);
                    for (Idx col_k = 0; col_k < n; ++col_k) root_lp.change_col_bounds(col_k, cur.col_lo[col_k], cur.col_up[col_k]);
                    root_lp.change_col_bounds(j, cur.col_lo[j], std::floor(c.val));
                    Status st_d = root_lp.solve(deadline);
                    res.lp_iterations += root_lp.stats.iterations;
                    bool down_inf = (st_d == Status::Infeasible || st_d == Status::CutoffReached);
                    d_down = down_inf ? 1e6 : std::max(0.0, root_lp.objective - root_lp.obj_offset - (z_node - root_lp.obj_offset));

                    // Up trial: x_j >= ceil(v)
                    root_lp.set_basis(node_basis);
                    for (Idx col_k = 0; col_k < n; ++col_k) root_lp.change_col_bounds(col_k, cur.col_lo[col_k], cur.col_up[col_k]);
                    root_lp.change_col_bounds(j, std::ceil(c.val), cur.col_up[j]);
                    Status st_u = root_lp.solve(deadline);
                    res.lp_iterations += root_lp.stats.iterations;
                    bool up_inf = (st_u == Status::Infeasible || st_u == Status::CutoffReached);
                    d_up = up_inf ? 1e6 : std::max(0.0, root_lp.objective - root_lp.obj_offset - (z_node - root_lp.obj_offset));

                    // Restore node basis and column bounds
                    root_lp.set_basis(node_basis);
                    for (Idx col_k = 0; col_k < n; ++col_k) root_lp.change_col_bounds(col_k, cur.col_lo[col_k], cur.col_up[col_k]);

                    if (down_inf && up_inf) {
                        bound_tightened_by_sb = true;
                        best_branch_var = kNone;
                        break; // Whole node pruned!
                    } else if (down_inf) {
                        cur.col_lo[j] = std::ceil(c.val);
                        bound_tightened_by_sb = true;
                        pq.push(cur);
                        break;
                    } else if (up_inf) {
                        cur.col_up[j] = std::floor(c.val);
                        bound_tightened_by_sb = true;
                        pq.push(cur);
                        break;
                    }

                    pcost[j].sum_down += d_down / c.f_down;
                    pcost[j].count_down++;
                    pcost[j].sum_up += d_up / c.f_up;
                    pcost[j].count_up++;
                } else {
                    double psd = (pcost[j].count_down > 0) ? (pcost[j].sum_down / pcost[j].count_down) : avg_ps_down;
                    double psu = (pcost[j].count_up > 0) ? (pcost[j].sum_up / pcost[j].count_up) : avg_ps_up;
                    d_down = psd * c.f_down;
                    d_up = psu * c.f_up;
                }

                double score = std::max(1e-6, d_down) * std::max(1e-6, d_up);
                if (score > best_score) {
                    best_score = score;
                    best_branch_var = j;
                    best_branch_val = c.val;
                }
            }

            if (bound_tightened_by_sb || best_branch_var == kNone) continue;

            // Create down and up child nodes
            Node down_child = cur;
            down_child.id = next_node_id++;
            down_child.depth = cur.depth + 1;
            down_child.lower_bound = z_node;
            down_child.col_up[best_branch_var] = std::floor(best_branch_val);
            down_child.estimate = z_node + 0.1 * (best_branch_val - std::floor(best_branch_val));
            down_child.basis = node_basis;

            Node up_child = cur;
            up_child.id = next_node_id++;
            up_child.depth = cur.depth + 1;
            up_child.lower_bound = z_node;
            up_child.col_lo[best_branch_var] = std::ceil(best_branch_val);
            up_child.estimate = z_node + 0.1 * (std::ceil(best_branch_val) - best_branch_val);
            up_child.basis = node_basis;

            pq.push(std::move(down_child));
            pq.push(std::move(up_child));

            if ((res.nodes_explored % 50) == 0) {
                double cur_gap = res.has_incumbent
                    ? 100.0 * std::max(0.0, res.incumbent_obj - res.best_bound) / std::max(1e-10, std::abs(res.incumbent_obj))
                    : kInf;
                logger.log(2, " %6d | %5zu | %4d | %12.6g | %12.6g | %6.2f%% | %.2fs\n",
                           res.nodes_explored, pq.size(), cur.depth,
                           model.user_objective_from_internal(res.incumbent_obj),
                           model.user_objective_from_internal(res.best_bound),
                           cur_gap, deadline.timer.elapsed_sec());
            }
        }

        if (pq.empty()) {
            if (res.has_incumbent) {
                res.best_bound = res.incumbent_obj;
                res.gap = 0.0;
                res.status = Status::Optimal;
            } else {
                res.status = Status::Infeasible;
            }
        } else if (res.has_incumbent) {
            double min_b = res.incumbent_obj;
            while (!pq.empty()) {
                min_b = std::min(min_b, pq.top().lower_bound);
                pq.pop();
            }
            res.best_bound = std::min(res.incumbent_obj, min_b);
            res.gap = std::max(0.0, res.incumbent_obj - res.best_bound) / std::max(1e-10, std::abs(res.incumbent_obj));
            if (res.gap <= opt.tol.mip_gap || (res.incumbent_obj - res.best_bound) <= opt.tol.mip_abs_gap) {
                res.status = Status::Optimal;
            }
        }

        return res;
    }
};

} // namespace sov
