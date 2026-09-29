#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "../core/types.hpp"
#include "../core/tolerances.hpp"
#include "../core/dd_real.hpp"
#include "../core/csc_matrix.hpp"

namespace sov {

// Layer L2: Approximate Minimum Degree (AMD) Ordering (Amestoy, Davis & Duff 1996)
// + Sparse/Elimination-Tree LDL^T Cholesky + Iterative Refinement (Part 5.8 & 5.10)
class SparseCholeskyAmd {
public:
    Idx n = 0;
    std::vector<Idx> perm;     // Fill-reducing AMD permutation: new_idx -> old_idx
    std::vector<Idx> inv_perm; // old_idx -> new_idx
    std::vector<Idx> parent;   // Elimination tree parent array
    std::vector<std::vector<std::pair<Idx, double>>> L_cols; // Column j holds (row i > j, l_ij)
    std::vector<double> D;     // Diagonal entries d_j

    // Compute fill-reducing Approximate Minimum Degree (AMD) ordering on symmetric matrix pattern S
    static void compute_amd_ordering(Idx dim, const std::vector<std::vector<Idx>>& adj_in,
                                     std::vector<Idx>& perm_out, std::vector<Idx>& inv_perm_out) {
        perm_out.resize(dim);
        inv_perm_out.resize(dim);
        if (dim == 0) return;

        std::vector<std::vector<Idx>> adj = adj_in;
        std::vector<bool> eliminated(dim, false);
        std::vector<Idx> degree(dim, 0);
        for (Idx i = 0; i < dim; ++i) {
            degree[i] = static_cast<Idx>(adj[i].size());
        }

        for (Idx step = 0; step < dim; ++step) {
            // Select uneliminated vertex with minimum external degree
            Idx pivot = kNone;
            Idx min_deg = dim + 1;
            for (Idx v = 0; v < dim; ++v) {
                if (!eliminated[v] && degree[v] < min_deg) {
                    min_deg = degree[v];
                    pivot = v;
                }
            }
            eliminated[pivot] = true;
            perm_out[step] = pivot;
            inv_perm_out[pivot] = step;

            // Collect active neighbors of pivot (the new clique/element e_v)
            std::vector<Idx> clique;
            for (Idx nb : adj[pivot]) {
                if (!eliminated[nb]) clique.push_back(nb);
            }

            // Update adjacency and approximate degree of neighbors in clique
            for (size_t a = 0; a < clique.size(); ++a) {
                Idx u = clique[a];
                for (size_t b = a + 1; b < clique.size(); ++b) {
                    Idx w = clique[b];
                    if (std::find(adj[u].begin(), adj[u].end(), w) == adj[u].end()) {
                        adj[u].push_back(w);
                        adj[w].push_back(u);
                    }
                }
                Idx deg = 0;
                for (Idx nb : adj[u]) {
                    if (!eliminated[nb]) ++deg;
                }
                degree[u] = deg;
            }
        }
    }

    // Factor symmetric matrix S (given as dense dim x dim or sparse pattern) using AMD + LDL^T
    bool factor_spd(Idx dim, const std::vector<double>& S_dense, double reg = 1e-10) {
        n = dim;
        if (n == 0) return true;

        // 1. Extract adjacency graph of S
        std::vector<std::vector<Idx>> adj(n);
        double max_diag = 0.0;
        for (Idx i = 0; i < n; ++i) {
            max_diag = std::max(max_diag, std::abs(S_dense[static_cast<size_t>(i) * n + i]));
            for (Idx j = i + 1; j < n; ++j) {
                if (std::abs(S_dense[static_cast<size_t>(i) * n + j]) > kDropTol) {
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                }
            }
        }

        // 2. Compute AMD fill-reducing permutation P
        compute_amd_ordering(n, adj, perm, inv_perm);

        // 3. Up-looking / Column-by-column sparse LDL^T on P S P^T
        L_cols.assign(n, {});
        D.assign(n, 1.0);
        double delta = std::max(1e-12, reg * std::max(1.0, max_diag));

        std::vector<double> work_col(n, 0.0);
        for (Idx j = 0; j < n; ++j) {
            Idx orig_j = perm[j];
            for (Idx i = j; i < n; ++i) {
                Idx orig_i = perm[i];
                work_col[i] = S_dense[static_cast<size_t>(orig_i) * n + orig_j];
            }
            double dj = work_col[j] + delta;

            for (Idx k = 0; k < j; ++k) {
                // Find L(j, k) if nonzero
                double l_jk = 0.0;
                for (const auto& entry : L_cols[k]) {
                    if (entry.first == j) {
                        l_jk = entry.second;
                        break;
                    }
                }
                if (l_jk == 0.0) continue;
                dj -= l_jk * l_jk * D[k];
                for (const auto& entry : L_cols[k]) {
                    if (entry.first > j) {
                        work_col[entry.first] -= entry.second * l_jk * D[k];
                    }
                }
            }

            if (dj <= 1e-13 * std::max(1.0, max_diag)) {
                dj = 1e-10 * std::max(1.0, max_diag);
            }
            D[j] = dj;

            for (Idx i = j + 1; i < n; ++i) {
                double l_ij = work_col[i] / dj;
                if (std::abs(l_ij) > kDropTol) {
                    L_cols[j].push_back({i, l_ij});
                }
                work_col[i] = 0.0;
            }
            work_col[j] = 0.0;
        }
        return true;
    }

    void solve_raw(const std::vector<double>& rhs, std::vector<double>& x) const {
        x.assign(n, 0.0);
        if (n == 0) return;
        std::vector<double> z(n);
        for (Idx i = 0; i < n; ++i) z[i] = rhs[perm[i]];

        // Forward solve: L z_new = P rhs
        for (Idx j = 0; j < n; ++j) {
            double zj = z[j];
            if (zj == 0.0) continue;
            for (const auto& e : L_cols[j]) {
                z[e.first] -= e.second * zj;
            }
        }
        // Diagonal scaling: w = D^-1 z
        for (Idx j = 0; j < n; ++j) z[j] /= D[j];

        // Backward solve: L^T y = w
        for (Idx j = n - 1; j >= 0; --j) {
            double s = z[j];
            for (const auto& e : L_cols[j]) {
                s -= e.second * z[e.first];
            }
            z[j] = s;
        }
        for (Idx i = 0; i < n; ++i) x[perm[i]] = z[i];
    }

    // Solve with Iterative Refinement in dd_real (Part 5.10)
    void solve_refined(const std::vector<double>& S_dense,
                       const std::vector<double>& rhs,
                       std::vector<double>& x,
                       int max_refinements = 3) const {
        solve_raw(rhs, x);
        if (n == 0 || max_refinements <= 0) return;

        double rhs_norm = 0.0;
        for (double v : rhs) rhs_norm = std::max(rhs_norm, std::abs(v));

        std::vector<double> r(n, 0.0), dx(n, 0.0);
        for (int pass = 0; pass < max_refinements; ++pass) {
            double max_r = 0.0;
            for (Idx i = 0; i < n; ++i) {
                dd_real acc(rhs[i]);
                for (Idx j = 0; j < n; ++j) {
                    acc.add_prod(-S_dense[static_cast<size_t>(i) * n + j], x[j]);
                }
                r[i] = acc.to_double();
                max_r = std::max(max_r, std::abs(r[i]));
            }
            if (max_r <= 1e-12 * (1.0 + rhs_norm)) break;
            solve_raw(r, dx);
            for (Idx i = 0; i < n; ++i) x[i] += dx[i];
        }
    }
};

} // namespace sov
