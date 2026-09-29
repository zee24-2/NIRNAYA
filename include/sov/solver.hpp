#pragma once
#include <string>
#include <vector>
#include <queue>
#include <atomic>
#include <fstream>
#include "../../src/core/types.hpp"
#include "../../src/core/status.hpp"
#include "../../src/core/options.hpp"
#include "../../src/core/timer.hpp"
#include "../../src/core/log.hpp"
#include "../../src/model/model.hpp"
#include "../../src/model/scaling.hpp"
#include "../../src/io/mps_reader.hpp"
#include "../../src/io/sol_io.hpp"
#include "../../src/presolve/presolve.hpp"
#include "../../src/lp/lp_solver.hpp"
#include "../../src/ipm/ipm_lp_qp.hpp"
#include "../../src/mip/mip_solver.hpp"
#include "../../src/gpu/pdhg_solver.hpp"

namespace sov {

struct SolveResult {
    Status status = Status::InternalError;
    double objective = kInf;
    double best_bound = -kInf;
    double gap = 0.0;
    bool has_incumbent = false;
    std::vector<double> x;
    std::vector<double> row_activity;
    std::vector<double> row_duals;
    std::vector<double> reduced_costs;
    Idx iterations = 0;
    Idx nodes = 0;
    Idx cuts_added = 0;
    double solve_time_sec = 0.0;
    VerificationReport verification;
    std::string message;
};

class Solver {
public:
    Options options;
    std::atomic<bool> stop_flag{false};

    explicit Solver(const Options& opt = Options()) : options(opt) {}

    void interrupt() {
        stop_flag.store(true, std::memory_order_relaxed);
    }

    // Convex MIQP Branch-and-Bound Engine (Slide 2 & Slide 3 Extension-Ready MIQP)
    static MipResult solve_miqp(const Model& model, const Options& opt, const Deadline& deadline) {
        MipResult res;
        Idx n = model.ncol;
        std::vector<Idx> int_vars;
        for (Idx j = 0; j < n; ++j) {
            if (model.vartype[j] != VarType::Continuous) int_vars.push_back(j);
        }

        struct QpNode {
            double lb = -kInf;
            std::vector<double> lo;
            std::vector<double> up;
        };
        auto cmp = [](const QpNode& a, const QpNode& b) { return a.lb > b.lb; };
        std::priority_queue<QpNode, std::vector<QpNode>, decltype(cmp)> pq(cmp);
        pq.push({-kInf, model.collo, model.colup});

        while (!pq.empty() && !deadline.expired() && res.nodes_explored < opt.node_limit) {
            QpNode cur = pq.top();
            pq.pop();
            if (res.has_incumbent && cur.lb >= res.incumbent_obj - opt.tol.mip_abs_gap) continue;
            ++res.nodes_explored;

            Model node_m = model;
            node_m.collo = cur.lo;
            node_m.colup = cur.up;
            auto ipm = IpmSolver::solve(node_m, opt, deadline);
            res.lp_iterations += ipm.iterations;
            if (ipm.status != Status::Optimal) continue;
            auto rep_node = verify_solution(node_m, ipm.x, ipm.objective, 1e-4);
            if (rep_node.max_row_viol > 1e-4 || rep_node.max_bound_viol > 1e-4) continue;
            if (!res.has_incumbent || res.nodes_explored == 1) {
                res.best_bound = std::max(res.best_bound, ipm.objective);
            }
            if (res.has_incumbent && ipm.objective >= res.incumbent_obj - opt.tol.mip_abs_gap) continue;

            Idx branch_j = kNone;
            double max_frac = 0.0;
            for (Idx j : int_vars) {
                double f = std::min(ipm.x[j] - std::floor(ipm.x[j]), std::ceil(ipm.x[j]) - ipm.x[j]);
                if (f > opt.tol.int_tol && cur.lo[j] < cur.up[j] && f > max_frac) {
                    max_frac = f;
                    branch_j = j;
                }
            }
            if (branch_j == kNone) {
                std::vector<double> x_int = ipm.x;
                for (Idx j : int_vars) x_int[j] = std::round(x_int[j]);
                double obj = model.eval_objective(x_int);
                if (verify_solution(model, x_int, obj, 1e-4).passed) {
                    if (!res.has_incumbent || obj < res.incumbent_obj) {
                        res.has_incumbent = true;
                        res.incumbent_obj = obj;
                        res.incumbent_x = x_int;
                    }
                }
                continue;
            }
            QpNode left = cur, right = cur;
            left.lb = ipm.objective;
            left.up[branch_j] = std::floor(ipm.x[branch_j]);
            right.lb = ipm.objective;
            right.lo[branch_j] = std::ceil(ipm.x[branch_j]);
            pq.push(std::move(left));
            pq.push(std::move(right));
        }
        if (res.has_incumbent) {
            res.status = Status::Optimal;
            res.best_bound = std::min(res.incumbent_obj, res.best_bound);
            res.gap = std::max(0.0, res.incumbent_obj - res.best_bound) / std::max(1.0, std::abs(res.incumbent_obj));
        }
        return res;
    }

    SolveResult solve(const Model& user_model) {
        SolveResult out;
        stop_flag.store(false, std::memory_order_relaxed);
        Logger logger(options.log_level);

        Deadline deadline;
        deadline.wall_limit_sec = options.time_limit;
        deadline.work_limit = options.work_limit;
        deadline.stop_flag = &stop_flag;

        Model valid_copy = user_model;
        Status vstat = valid_copy.validate(out.message);
        if (vstat != Status::Optimal) {
            out.status = vstat;
            out.solve_time_sec = deadline.timer.elapsed_sec();
            return out;
        }

        logger.log(1, "========================================================================\n");
        logger.log(1, "  NIRNAYA -- SOVEREIGN OPTIMIZATION SOLVER CORE (SIH26119 | Team 151198)\n");
        logger.log(1, "  Model: %-18s | Rows: %d | Cols: %d | NNZ: %lld | Ints: %d | Q-NNZ: %lld\n",
                   valid_copy.name.c_str(), valid_copy.nrow, valid_copy.ncol,
                   static_cast<long long>(valid_copy.A.nnz()),
                   valid_copy.num_integers(),
                   static_cast<long long>(valid_copy.Q.nnz()));
        logger.log(1, "========================================================================\n");

        // Convert to internal minimization form (C1)
        Model min_model = valid_copy.to_internal_min();

        // Step 3: Presolve (Part 7)
        auto pre = Presolver::presolve(min_model, options.presolve, options.presolve_rounds);
        if (options.presolve > 0 && (pre.rows_removed > 0 || pre.cols_removed > 0 || pre.bound_tightenings > 0)) {
            logger.log(1, "Presolve: %d rows, %d cols remaining (removed %d rows, %d cols, %d tightenings)\n",
                       pre.reduced.nrow, pre.reduced.ncol,
                       pre.rows_removed, pre.cols_removed, pre.bound_tightenings);
        }

        if (pre.solved_in_presolve) {
            out.status = pre.status;
            if (pre.status == Status::Optimal) {
                pre.postsolve({}, out.x);
                out.has_incumbent = true;
                out.objective = valid_copy.eval_objective(out.x);
                out.best_bound = out.objective;
                out.verification = verify_solution(valid_copy, out.x, out.objective, options.tol.checker_tol);
            }
            out.solve_time_sec = deadline.timer.elapsed_sec();
            return out;
        }

        // Step 4: Scale (Part 4)
        ScalingInfo sinfo;
        Model scaled_model = scale_model(pre.reduced, options.scale, sinfo);

        // Step 5 & 6: Dispatch to MIQP, QP, MILP, or LP (Dual Simplex / Primal Simplex / IPM / PDHG)
        std::vector<double> x_red(pre.reduced.ncol, 0.0);
        std::vector<double> y_red(pre.reduced.nrow, 0.0);
        std::vector<double> d_red(pre.reduced.ncol, 0.0);

        if (scaled_model.is_qp() && scaled_model.is_mip()) {
            logger.log(1, "Solving Convex MIQP via Branch-and-Bound + Mehrotra IPM...\n");
            auto miqp_res = solve_miqp(pre.reduced, options, deadline);
            out.status = miqp_res.status;
            out.nodes = miqp_res.nodes_explored;
            out.iterations = miqp_res.lp_iterations;
            out.has_incumbent = miqp_res.has_incumbent;
            if (miqp_res.has_incumbent) {
                x_red = miqp_res.incumbent_x;
                out.gap = miqp_res.gap;
            }
            out.best_bound = valid_copy.user_objective_from_internal(miqp_res.best_bound);
        } else if (scaled_model.is_qp()) {
            logger.log(1, "Solving Convex QP via Mehrotra Predictor-Corrector IPM + AMD Cholesky...\n");
            auto ipm_res = IpmSolver::solve(scaled_model, options, deadline);
            out.status = ipm_res.status;
            out.iterations = ipm_res.iterations;
            if (out.status == Status::Optimal) {
                sinfo.unscale_primal(ipm_res.x, x_red);
                sinfo.unscale_row_duals(ipm_res.y, y_red);
                out.has_incumbent = true;
            }
        } else if (scaled_model.is_mip()) {
            logger.log(1, "Solving MILP via Branch-and-Cut...\n");
            auto mip_res = MipSolver::solve(pre.reduced, options, deadline, logger);
            out.status = mip_res.status;
            out.nodes = mip_res.nodes_explored;
            out.iterations = mip_res.lp_iterations;
            out.cuts_added = mip_res.cuts_added;
            out.has_incumbent = mip_res.has_incumbent;
            if (mip_res.has_incumbent) {
                x_red = mip_res.incumbent_x;
                out.gap = mip_res.gap;
            }
            out.best_bound = valid_copy.user_objective_from_internal(mip_res.best_bound);
        } else {
            // Continuous LP
            if (options.method == "ipm") {
                logger.log(1, "Solving LP via Mehrotra Predictor-Corrector IPM + AMD Cholesky...\n");
                auto ipm_res = IpmSolver::solve(scaled_model, options, deadline);
                out.status = ipm_res.status;
                out.iterations = ipm_res.iterations;
                if (out.status == Status::Optimal) {
                    sinfo.unscale_primal(ipm_res.x, x_red);
                    sinfo.unscale_row_duals(ipm_res.y, y_red);
                    out.has_incumbent = true;
                }
            } else if (options.method == "pdhg" || options.method == "gpu") {
                logger.log(1, "Solving LP via Parallel First-Order PDHG (Chambolle-Pock PDLP)...\n");
                auto pdhg_res = PdhgSolver::solve(scaled_model, options, deadline);
                out.status = pdhg_res.status;
                out.iterations = pdhg_res.iterations;
                if (out.status == Status::Optimal) {
                    sinfo.unscale_primal(pdhg_res.x, x_red);
                    sinfo.unscale_row_duals(pdhg_res.y, y_red);
                    out.has_incumbent = true;
                }
            } else {
                logger.log(1, "Solving LP via Dual/Primal Simplex with Sparse Markowitz LU...\n");
                LpSolver lp(options);
                lp.load(scaled_model);
                out.status = lp.solve(deadline);
                out.iterations = lp.stats.iterations;
                if (out.status == Status::Optimal) {
                    std::vector<double> x_sc(lp.value.begin(), lp.value.begin() + scaled_model.ncol);
                    std::vector<double> d_sc(lp.dual.begin(), lp.dual.begin() + scaled_model.ncol);
                    sinfo.unscale_primal(x_sc, x_red);
                    sinfo.unscale_row_duals(lp.row_duals, y_red);
                    sinfo.unscale_reduced_costs(d_sc, d_red);
                    out.has_incumbent = true;
                }
            }
        }

        // Step 7: Unscale + Postsolve (Part 7.5)
        if (out.has_incumbent) {
            pre.postsolve(x_red, out.x);
            valid_copy.A.spmv(out.x, out.row_activity);
            out.objective = valid_copy.eval_objective(out.x);
            if (!scaled_model.is_mip()) {
                out.best_bound = out.objective;
                out.gap = 0.0;
                out.row_duals.assign(valid_copy.nrow, 0.0);
                for (size_t i = 0; i < pre.orig_row_map.size() && i < y_red.size(); ++i) {
                    double sign = (valid_copy.sense == Sense::Maximize) ? -1.0 : 1.0;
                    out.row_duals[pre.orig_row_map[i]] = sign * y_red[i];
                }
                out.reduced_costs.assign(valid_copy.ncol, 0.0);
                for (size_t j = 0; j < pre.orig_col_map.size() && j < d_red.size(); ++j) {
                    double sign = (valid_copy.sense == Sense::Maximize) ? -1.0 : 1.0;
                    out.reduced_costs[pre.orig_col_map[j]] = sign * d_red[j];
                }
            }

            // Step 8: Independent dd_real Verification (Part 12.6)
            out.verification = verify_solution(valid_copy, out.x, out.objective, options.tol.checker_tol);
            if (!out.verification.passed && out.status == Status::Optimal) {
                if (options.presolve > 0 || options.scale) {
                    Options fallback_opt = options;
                    fallback_opt.presolve = 0;
                    fallback_opt.scale = false;
                    fallback_opt.log_level = 0;
                    Solver fallback_solver(fallback_opt);
                    auto fb = fallback_solver.solve(user_model);
                    if (fb.verification.passed && fb.status == Status::Optimal) {
                        fb.solve_time_sec = deadline.timer.elapsed_sec();
                        return fb;
                    }
                }
                out.status = Status::NumericalTrouble;
                out.message = "Verification residual exceeded tolerance on " + out.verification.worst_item;
            }
        }

        out.solve_time_sec = deadline.timer.elapsed_sec();

        logger.log(1, "------------------------------------------------------------------------\n");
        logger.log(1, "  Status:     %s\n", status_to_string(out.status));
        if (out.has_incumbent) {
            logger.log(1, "  Objective:  %.12g\n", out.objective);
            logger.log(1, "  Best Bound: %.12g | Gap: %.4f%%\n", out.best_bound, 100.0 * out.gap);
            logger.log(1, "  Check (dd): PASS=%s | MaxRowViol=%.2e | MaxBndViol=%.2e | MaxIntViol=%.2e\n",
                       out.verification.passed ? "YES" : "NO",
                       out.verification.max_row_viol,
                       out.verification.max_bound_viol,
                       out.verification.max_int_viol);
        }
        logger.log(1, "  Iterations: %d | Nodes: %d | Cuts: %d | Time: %.4f s\n",
                   out.iterations, out.nodes, out.cuts_added, out.solve_time_sec);
        logger.log(1, "========================================================================\n");

        if (!options.sol_path.empty() && out.has_incumbent) {
            SolIo::write_sol(options.sol_path, valid_copy, out.status, out.objective, out.gap, out.x);
        }
        if (!options.stats_json.empty()) {
            write_stats_json(options.stats_json, valid_copy, out);
        }
        return out;
    }

    static bool write_stats_json(const std::string& path, const Model& m, const SolveResult& r) {
        std::ofstream out(path);
        if (!out.is_open()) return false;
        out << "{\n";
        out << "  \"solver\": \"NIRNAYA\",\n";
        out << "  \"problem_statement\": \"SIH26119\",\n";
        out << "  \"model\": \"" << m.name << "\",\n";
        out << "  \"status\": \"" << status_to_string(r.status) << "\",\n";
        out << "  \"objective\": " << (std::isfinite(r.objective) ? r.objective : 0.0) << ",\n";
        out << "  \"best_bound\": " << (std::isfinite(r.best_bound) ? r.best_bound : 0.0) << ",\n";
        out << "  \"gap\": " << (std::isfinite(r.gap) ? r.gap : 0.0) << ",\n";
        out << "  \"iterations\": " << r.iterations << ",\n";
        out << "  \"nodes\": " << r.nodes << ",\n";
        out << "  \"cuts_added\": " << r.cuts_added << ",\n";
        out << "  \"solve_time_sec\": " << r.solve_time_sec << ",\n";
        out << "  \"verified\": " << (r.verification.passed ? "true" : "false") << "\n";
        out << "}\n";
        return true;
    }
};

} // namespace sov
