#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/status.hpp"
#include "../core/options.hpp"
#include "../core/timer.hpp"
#include "../core/thread_pool.hpp"
#include "../model/model.hpp"

namespace sov {

struct PdhgResult {
    Status status = Status::Optimal;
    double objective = 0.0;
    std::vector<double> x;
    std::vector<double> y;
    Idx iterations = 0;
};

// GPU / Multi-Core Parallel First-Order PDHG (Primal-Dual Hybrid Gradient / PDLP) Solver (Part 11.4 & SIH26119)
class PdhgSolver {
public:
    static PdhgResult solve(const Model& model, const Options& opt, const Deadline& deadline) {
        PdhgResult res;
        Idx n = model.ncol;
        Idx m = model.nrow;
        CsrMatrix csr = CsrMatrix::from_csc(model.A);
        ThreadPool pool(opt.threads);

        // Estimate ||A||_2 via power iteration
        std::vector<double> u(n, 1.0), v(m, 0.0);
        double norm_est = 1.0;
        for (int p = 0; p < 12; ++p) {
            model.A.spmv(u, v);
            model.A.spmv_transpose(v, u);
            double n2 = 0.0;
            for (double val : u) n2 += val * val;
            norm_est = std::sqrt(std::max(1e-12, n2));
            for (double& val : u) val /= norm_est;
        }
        double sigma_max = std::sqrt(norm_est);
        double tau   = 0.85 / std::max(1e-6, sigma_max);
        double sigma = 0.85 / std::max(1e-6, sigma_max);

        std::vector<double> x(n, 0.0), x_next(n, 0.0), x_bar(n, 0.0);
        std::vector<double> y(m, 0.0), r_proj(m, 0.0);
        for (Idx j = 0; j < n; ++j) {
            if (is_finite_bound(model.collo[j]) && is_finite_bound(model.colup[j])) {
                x[j] = 0.5 * (model.collo[j] + model.colup[j]);
            } else if (is_finite_bound(model.collo[j])) {
                x[j] = std::max(0.0, model.collo[j]);
            } else if (is_finite_bound(model.colup[j])) {
                x[j] = std::min(0.0, model.colup[j]);
            }
        }

        std::vector<double> ATy(n, 0.0), Ax_bar(m, 0.0);
        Idx max_iters = std::max<Idx>(5000, opt.simplex_iter_limit / 10);

        for (Idx iter = 0; iter < max_iters; ++iter) {
            if ((iter & 255) == 0 && deadline.expired()) {
                res.status = Status::TimeLimit;
                break;
            }
            res.iterations = iter + 1;

            // 1. Primal step: x^{k+1} = proj_[l, u]( x^k - tau * (c - A^T y^k) )
            model.A.spmv_transpose(y, ATy);
            pool.parallel_for(0, n, [&](Idx j) {
                double cand = x[j] - tau * (model.cost[j] - ATy[j]);
                if (is_finite_bound(model.collo[j])) cand = std::max(model.collo[j], cand);
                if (is_finite_bound(model.colup[j])) cand = std::min(model.colup[j], cand);
                x_next[j] = cand;
                x_bar[j] = 2.0 * cand - x[j];
                x[j] = cand;
            });

            // 2. Dual step on row activities r in [rowlo, rowup]:
            //    y^{k+1} = y^k + sigma * ( proj_[rowlo, rowup](A x_bar - y^k / sigma) - A x_bar )
            model.A.spmv(x_bar, Ax_bar);
            double max_p_viol = 0.0;
            for (Idx i = 0; i < m; ++i) {
                double target = Ax_bar[i] - y[i] / sigma;
                double r_clamped = target;
                if (is_finite_bound(model.rowlo[i])) r_clamped = std::max(model.rowlo[i], r_clamped);
                if (is_finite_bound(model.rowup[i])) r_clamped = std::min(model.rowup[i], r_clamped);
                y[i] += sigma * (r_clamped - Ax_bar[i]);

                double viol = 0.0;
                if (is_finite_bound(model.rowlo[i])) viol = std::max(viol, model.rowlo[i] - Ax_bar[i]);
                if (is_finite_bound(model.rowup[i])) viol = std::max(viol, Ax_bar[i] - model.rowup[i]);
                max_p_viol = std::max(max_p_viol, viol);
            }

            if (iter > 20 && (iter & 31) == 0 && max_p_viol <= 1e-6) {
                res.status = Status::Optimal;
                break;
            }
        }

        res.x = x;
        res.y = y;
        res.objective = model.eval_objective(res.x);
        return res;
    }
};

} // namespace sov
