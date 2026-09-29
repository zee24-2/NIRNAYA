#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../core/types.hpp"
#include "../core/tolerances.hpp"
#include "../core/hvector.hpp"
#include "../core/csc_matrix.hpp"

namespace sov {

struct BasisRepair {
    Idx basis_pos;
    Idx out_var;
    Idx in_logical_row;
};

enum class UpdateStatus : int8_t {
    Ok       = 0,
    Unstable = 1,
    Singular = 2
};

struct FactorResult {
    bool ok = true;
    Idx rank = 0;
    std::vector<BasisRepair> repairs;
};

// Sparse Markowitz LU Factorization + Forrest-Tomlin / Product-Form Eta Updates (Part 5)
class BasisFactor {
public:
    Idx m = 0;
    Idx n_struct = 0;
    double threshold_u = 0.1;
    int max_updates = 100;

    // L-etas: E_k eliminates row i using pivot row p_k with multiplier l_ik
    struct LEtaEntry {
        Idx row;
        double mult;
    };
    struct LEta {
        Idx pivot_row;
        std::vector<LEtaEntry> entries;
    };

    // U storage in pivot order k = 0..m-1:
    // Pivot k is at (row = pivot_row[k], basis_col = pivot_col[k]), diagonal u_diag[k]
    std::vector<Idx> pivot_row;      // k -> row p_k
    std::vector<Idx> pivot_col;      // k -> basis position q_k
    std::vector<Idx> row_to_step;    // row i -> step k
    std::vector<Idx> col_to_step;    // basis pos q -> step k
    std::vector<double> u_diag;      // k -> U(p_k, q_k)

    // U off-diagonal entries:
    // u_col_entries[k] holds entries (step j < k, val) in column pivot_col[k] at row pivot_row[j]
    struct UEntry {
        Idx step;
        double val;
    };
    std::vector<std::vector<UEntry>> u_col_entries; // indexed by step k (for FTRAN backsolve)
    std::vector<std::vector<UEntry>> u_row_entries; // indexed by step j (entries at step k > j for BTRAN forwardsolve)

    std::vector<LEta> l_etas;

    // Update etas (Product-Form / Forrest-Tomlin spike etas on basis coordinates):
    // B_new^{-1} = E_upd * B_old^{-1}
    struct UpdateEta {
        Idx r;                  // leaving basis position
        double inv_pivot;       // 1.0 / alpha_rq
        std::vector<Idx> idx;   // i != r with nonzeros
        std::vector<double> val;// -alpha_iq / alpha_rq
    };
    std::vector<UpdateEta> upd_etas;

    bool unstable_flag = false;
    Off nnz_L = 0;
    Off nnz_U = 0;
    Off initial_nnz_LU = 0;

    FactorResult factor(const CscMatrix& A, const std::vector<Idx>& basic_index) {
        m = A.nrow;
        n_struct = A.ncol;
        FactorResult res;
        res.rank = m;

        pivot_row.assign(m, kNone);
        pivot_col.assign(m, kNone);
        row_to_step.assign(m, kNone);
        col_to_step.assign(m, kNone);
        u_diag.assign(m, 1.0);
        u_col_entries.assign(m, {});
        u_row_entries.assign(m, {});
        l_etas.clear();
        upd_etas.clear();
        unstable_flag = false;

        if (m == 0) return res;

        // Active submatrix stored both column-wise (values) and row-wise (indices + values)
        std::vector<std::vector<std::pair<Idx, double>>> col_mat(m);
        std::vector<std::vector<std::pair<Idx, double>>> row_mat(m);
        std::vector<bool> row_active(m, true);
        std::vector<bool> col_active(m, true);

        for (Idx q = 0; q < m; ++q) {
            Idx var = basic_index[q];
            if (var < n_struct) {
                for (Off k = A.start[var]; k < A.start[var + 1]; ++k) {
                    Idx r = A.index[k];
                    double v = A.value[k];
                    if (std::abs(v) > kDropTol) {
                        col_mat[q].push_back({r, v});
                        row_mat[r].push_back({q, v});
                    }
                }
            } else {
                // Logical column n_struct + row_i is -e_{row_i}
                Idx r = var - n_struct;
                col_mat[q].push_back({r, -1.0});
                row_mat[r].push_back({q, -1.0});
            }
        }

        auto get_active_row_count = [&](Idx r) -> Idx {
            Idx c = 0;
            for (const auto& p : row_mat[r]) if (col_active[p.first]) ++c;
            return c;
        };

        for (Idx step = 0; step < m; ++step) {
            Idx best_r = kNone;
            Idx best_q = kNone;
            double best_piv = 0.0;

            // Phase A1: Check for column singletons first
            for (Idx q = 0; q < m; ++q) {
                if (!col_active[q]) continue;
                Idx r_only = kNone;
                double v_only = 0.0;
                Idx cnt = 0;
                for (const auto& p : col_mat[q]) {
                    if (row_active[p.first]) {
                        r_only = p.first;
                        v_only = p.second;
                        if (++cnt > 1) break;
                    }
                }
                if (cnt == 1 && std::abs(v_only) >= 1e-9) {
                    best_r = r_only;
                    best_q = q;
                    best_piv = v_only;
                    break;
                }
            }

            // Phase A2: Check for row singletons
            if (best_r == kNone) {
                for (Idx r = 0; r < m; ++r) {
                    if (!row_active[r]) continue;
                    Idx q_only = kNone;
                    double v_only = 0.0;
                    Idx cnt = 0;
                    for (const auto& p : row_mat[r]) {
                        if (col_active[p.first]) {
                            q_only = p.first;
                            v_only = p.second;
                            if (++cnt > 1) break;
                        }
                    }
                    if (cnt == 1 && std::abs(v_only) >= 1e-9) {
                        best_r = r;
                        best_q = q_only;
                        best_piv = v_only;
                        break;
                    }
                }
            }

            // Phase B: Markowitz threshold pivoting on nucleus
            if (best_r == kNone) {
                int64_t best_merit = static_cast<int64_t>(m) * static_cast<int64_t>(m) + 1;
                for (Idx q = 0; q < m; ++q) {
                    if (!col_active[q]) continue;
                    double cmax = 0.0;
                    Idx c_cnt = 0;
                    for (const auto& p : col_mat[q]) {
                        if (row_active[p.first]) {
                            cmax = std::max(cmax, std::abs(p.second));
                            ++c_cnt;
                        }
                    }
                    if (c_cnt == 0 || cmax < 1e-9) continue;
                    double thresh = std::max(1e-9, threshold_u * cmax);
                    for (const auto& p : col_mat[q]) {
                        Idx r = p.first;
                        double v = p.second;
                        if (!row_active[r] || std::abs(v) < thresh) continue;
                        Idx r_cnt = get_active_row_count(r);
                        int64_t merit = static_cast<int64_t>(r_cnt - 1) * static_cast<int64_t>(c_cnt - 1);
                        if (merit < best_merit || (merit == best_merit && std::abs(v) > std::abs(best_piv))) {
                            best_merit = merit;
                            best_r = r;
                            best_q = q;
                            best_piv = v;
                        }
                    }
                }
            }

            // Rank deficiency / singularity repair (5.3.4)
            if (best_r == kNone) {
                for (Idx r = 0; r < m; ++r) {
                    if (row_active[r]) { best_r = r; break; }
                }
                for (Idx q = 0; q < m; ++q) {
                    if (col_active[q]) { best_q = q; break; }
                }
                res.ok = false;
                res.rank--;
                res.repairs.push_back({best_q, basic_index[best_q], best_r});
                best_piv = -1.0; // Replaced by logical -e_{best_r}
                col_mat[best_q].clear();
                col_mat[best_q].push_back({best_r, -1.0});
                row_mat[best_r].push_back({best_q, -1.0});
            }

            // Record pivot step
            pivot_row[step] = best_r;
            pivot_col[step] = best_q;
            row_to_step[best_r] = step;
            col_to_step[best_q] = step;
            u_diag[step] = best_piv;

            row_active[best_r] = false;
            col_active[best_q] = false;

            // Extract remaining active entries of pivot row best_r -> U row
            std::vector<std::pair<Idx, double>> piv_row_other;
            for (const auto& p : row_mat[best_r]) {
                if (col_active[p.first] && std::abs(p.second) > kDropTol) {
                    piv_row_other.push_back(p);
                }
            }

            // Eliminate remaining active entries in pivot column best_q
            LEta eta;
            eta.pivot_row = best_r;
            for (const auto& p : col_mat[best_q]) {
                Idx r_elim = p.first;
                if (!row_active[r_elim]) continue;
                double a_iq = p.second;
                if (std::abs(a_iq) <= kDropTol) continue;
                double mult = a_iq / best_piv;
                eta.entries.push_back({r_elim, mult});

                // Update row r_elim and corresponding columns
                for (const auto& pr : piv_row_other) {
                    Idx c_other = pr.first;
                    double delta = -mult * pr.second;
                    // Update in row_mat[r_elim]
                    bool found_r = false;
                    for (auto& entry : row_mat[r_elim]) {
                        if (entry.first == c_other) {
                            entry.second += delta;
                            found_r = true;
                            break;
                        }
                    }
                    if (!found_r && std::abs(delta) > kDropTol) {
                        row_mat[r_elim].push_back({c_other, delta});
                    }
                    // Update in col_mat[c_other]
                    bool found_c = false;
                    for (auto& entry : col_mat[c_other]) {
                        if (entry.first == r_elim) {
                            entry.second += delta;
                            found_c = true;
                            break;
                        }
                    }
                    if (!found_c && std::abs(delta) > kDropTol) {
                        col_mat[c_other].push_back({r_elim, delta});
                    }
                }
            }
            if (!eta.entries.empty()) {
                nnz_L += static_cast<Off>(eta.entries.size());
                l_etas.push_back(std::move(eta));
            }

            // Store temporarily in row_mat[best_r] so we can assemble U after all steps are numbered
            row_mat[best_r] = std::move(piv_row_other);
        }

        // Assemble U in step indices now that col_to_step is complete
        nnz_U = m;
        for (Idx step = 0; step < m; ++step) {
            Idx r = pivot_row[step];
            for (const auto& p : row_mat[r]) {
                Idx c_step = col_to_step[p.first];
                if (c_step > step && std::abs(p.second) > kDropTol) {
                    u_row_entries[step].push_back({c_step, p.second});
                    u_col_entries[c_step].push_back({step, p.second});
                    ++nnz_U;
                }
            }
        }
        initial_nnz_LU = nnz_L + nnz_U;
        return res;
    }

    // FTRAN: Solve B x = rhs. rhs enters in row indices (0..m-1) and exits in basis positions (0..m-1).
    void ftran(HVector& rhs, HVector* spike = nullptr) const {
        if (m == 0) return;
        std::vector<double> b_work = rhs.val;
        // 1. Apply L-etas: b := E b
        for (const auto& eta : l_etas) {
            double bp = b_work[eta.pivot_row];
            if (bp == 0.0) continue;
            for (const auto& e : eta.entries) {
                b_work[e.row] -= e.mult * bp;
            }
        }
        if (spike != nullptr) {
            std::fill(spike->val.begin(), spike->val.end(), 0.0);
            spike->count = 0;
            spike->packed_valid = true;
            for (Idx i = 0; i < m; ++i) {
                if (std::abs(b_work[i]) > kDropTol) {
                    spike->add(i, b_work[i]);
                }
            }
        }

        // 2. Solve U x_base = b_work in reverse pivot order k = m-1 down to 0
        std::vector<double> x_base(m, 0.0);
        for (Idx k = m - 1; k >= 0; --k) {
            Idx r = pivot_row[k];
            double br = b_work[r];
            if (std::abs(br) > kDropTol) {
                double xk = br / u_diag[k];
                x_base[pivot_col[k]] = xk;
                for (const auto& ue : u_col_entries[k]) {
                    Idx r_above = pivot_row[ue.step];
                    b_work[r_above] -= ue.val * xk;
                }
            }
        }

        // 3. Apply update etas in forward creation order
        for (const auto& ueta : upd_etas) {
            double xr = x_base[ueta.r];
            if (xr == 0.0) continue;
            x_base[ueta.r] = xr * ueta.inv_pivot;
            for (size_t i = 0; i < ueta.idx.size(); ++i) {
                x_base[ueta.idx[i]] += xr * ueta.val[i];
            }
        }

        std::fill(rhs.val.begin(), rhs.val.end(), 0.0);
        rhs.count = 0;
        rhs.packed_valid = true;
        for (Idx q = 0; q < m; ++q) {
            if (std::abs(x_base[q]) > kDropTol) {
                rhs.add(q, x_base[q]);
            }
        }
        rhs.tighten();
    }

    // BTRAN: Solve B^T y = rhs. rhs enters in basis positions (0..m-1) and exits in row indices (0..m-1).
    void btran(HVector& rhs) const {
        if (m == 0) return;
        std::vector<double> c_work(m, 0.0);
        if (rhs.packed_valid) {
            for (Idx k = 0; k < rhs.count; ++k) {
                Idx q = rhs.idx[k];
                c_work[q] = rhs.val[q];
            }
        } else {
            c_work = rhs.val;
        }

        // 1. Apply update etas transposed in reverse creation order
        for (int idx_u = static_cast<int>(upd_etas.size()) - 1; idx_u >= 0; --idx_u) {
            const auto& ueta = upd_etas[idx_u];
            double sum = c_work[ueta.r] * ueta.inv_pivot;
            for (size_t i = 0; i < ueta.idx.size(); ++i) {
                sum += ueta.val[i] * c_work[ueta.idx[i]];
            }
            c_work[ueta.r] = sum;
        }

        // 2. Solve U^T w = c_work in forward pivot order k = 0..m-1
        std::vector<double> y_row(m, 0.0);
        for (Idx k = 0; k < m; ++k) {
            Idx q = pivot_col[k];
            double cq = c_work[q];
            if (std::abs(cq) > kDropTol) {
                double wk = cq / u_diag[k];
                y_row[pivot_row[k]] = wk;
                for (const auto& ue : u_row_entries[k]) {
                    Idx q_later = pivot_col[ue.step];
                    c_work[q_later] -= ue.val * wk;
                }
            }
        }

        // 3. Apply L-etas transposed in reverse order
        for (int idx_l = static_cast<int>(l_etas.size()) - 1; idx_l >= 0; --idx_l) {
            const auto& eta = l_etas[idx_l];
            double sum = 0.0;
            for (const auto& e : eta.entries) {
                sum += e.mult * y_row[e.row];
            }
            y_row[eta.pivot_row] -= sum;
        }

        rhs.clear();
        for (Idx i = 0; i < m; ++i) {
            if (std::abs(y_row[i]) > kDropTol) {
                rhs.add(i, y_row[i]);
            }
        }
        rhs.tighten();
    }

    // Update factorization when basis position r is replaced by column with FTRAN result alpha_q
    UpdateStatus update(Idx r, const HVector& alpha_q, const HVector* /*spike*/ = nullptr) {
        double alpha_rq = alpha_q.val[r];
        if (std::abs(alpha_rq) < 1e-9) {
            return UpdateStatus::Singular;
        }
        UpdateEta ueta;
        ueta.r = r;
        ueta.inv_pivot = 1.0 / alpha_rq;
        double max_mult = std::abs(ueta.inv_pivot);
        for (Idx k = 0; k < alpha_q.count; ++k) {
            Idx i = alpha_q.idx[k];
            if (i == r) continue;
            double v = alpha_q.val[i];
            if (std::abs(v) > kDropTol) {
                double mult = -v / alpha_rq;
                ueta.idx.push_back(i);
                ueta.val.push_back(mult);
                max_mult = std::max(max_mult, std::abs(mult));
            }
        }
        upd_etas.push_back(std::move(ueta));
        if (max_mult > 1e6) {
            unstable_flag = true;
            return UpdateStatus::Unstable;
        }
        return UpdateStatus::Ok;
    }

    bool need_refactor() const {
        if (unstable_flag) return true;
        if (static_cast<int>(upd_etas.size()) >= max_updates) return true;
        return false;
    }

    int num_updates() const {
        return static_cast<int>(upd_etas.size());
    }
};

} // namespace sov
