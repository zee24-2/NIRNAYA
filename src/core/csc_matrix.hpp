#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "types.hpp"
#include "tolerances.hpp"
#include "hvector.hpp"

namespace sov {

struct Triplet {
    Idx row;
    Idx col;
    double val;
};

struct CscMatrix {
    Idx nrow = 0;
    Idx ncol = 0;
    std::vector<Off> start;
    std::vector<Idx> index;
    std::vector<double> value;

    CscMatrix() : start(1, 0) {}
    CscMatrix(Idx r, Idx c) : nrow(r), ncol(c), start(c + 1, 0) {}

    Off nnz() const { return start.empty() ? 0 : start[ncol]; }

    static CscMatrix from_triplets(Idx nrow, Idx ncol, const std::vector<Triplet>& triplets) {
        CscMatrix M(nrow, ncol);
        std::vector<Off> col_counts(ncol, 0);
        for (const auto& t : triplets) {
            if (t.row >= 0 && t.row < nrow && t.col >= 0 && t.col < ncol && t.val != 0.0) {
                col_counts[t.col]++;
            }
        }
        M.start[0] = 0;
        for (Idx j = 0; j < ncol; ++j) {
            M.start[j + 1] = M.start[j] + col_counts[j];
        }
        Off total = M.start[ncol];
        M.index.resize(total);
        M.value.resize(total);
        std::vector<Off> write_pos = M.start;
        for (const auto& t : triplets) {
            if (t.row >= 0 && t.row < nrow && t.col >= 0 && t.col < ncol && t.val != 0.0) {
                Off p = write_pos[t.col]++;
                M.index[p] = t.row;
                M.value[p] = t.val;
            }
        }
        // Sort each column by row index and sum duplicates
        std::vector<Off> new_start(ncol + 1, 0);
        std::vector<Idx> new_index;
        std::vector<double> new_value;
        new_index.reserve(total);
        new_value.reserve(total);
        std::vector<std::pair<Idx, double>> buf;

        for (Idx j = 0; j < ncol; ++j) {
            new_start[j] = static_cast<Off>(new_index.size());
            Off s = M.start[j], e = M.start[j + 1];
            if (s == e) continue;
            buf.clear();
            for (Off k = s; k < e; ++k) {
                buf.push_back({M.index[k], M.value[k]});
            }
            std::sort(buf.begin(), buf.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
            });
            Idx cur_r = buf[0].first;
            double cur_v = buf[0].second;
            for (size_t k = 1; k < buf.size(); ++k) {
                if (buf[k].first == cur_r) {
                    cur_v += buf[k].second;
                } else {
                    if (cur_v != 0.0) {
                        new_index.push_back(cur_r);
                        new_value.push_back(cur_v);
                    }
                    cur_r = buf[k].first;
                    cur_v = buf[k].second;
                }
            }
            if (cur_v != 0.0) {
                new_index.push_back(cur_r);
                new_value.push_back(cur_v);
            }
        }
        new_start[ncol] = static_cast<Off>(new_index.size());
        M.start = std::move(new_start);
        M.index = std::move(new_index);
        M.value = std::move(new_value);
        return M;
    }

    // y = A * x
    void spmv(const std::vector<double>& x, std::vector<double>& y) const {
        y.assign(nrow, 0.0);
        for (Idx j = 0; j < ncol; ++j) {
            double xj = x[j];
            if (xj == 0.0) continue;
            for (Off k = start[j]; k < start[j + 1]; ++k) {
                y[index[k]] += value[k] * xj;
            }
        }
    }

    // y = A^T * x
    void spmv_transpose(const std::vector<double>& x, std::vector<double>& y) const {
        y.assign(ncol, 0.0);
        for (Idx j = 0; j < ncol; ++j) {
            double s = 0.0;
            for (Off k = start[j]; k < start[j + 1]; ++k) {
                s += value[k] * x[index[k]];
            }
            y[j] = s;
        }
    }

    void extract_col(Idx j, HVector& out) const {
        out.clear();
        for (Off k = start[j]; k < start[j + 1]; ++k) {
            out.add(index[k], value[k]);
        }
        out.tighten(0.0);
    }

    std::string validate() const {
        if (static_cast<Idx>(start.size()) != ncol + 1) return "start.size() != ncol + 1";
        if (start[0] != 0) return "start[0] != 0";
        for (Idx j = 0; j < ncol; ++j) {
            if (start[j] > start[j + 1]) return "start not non-decreasing";
            Idx prev = -1;
            for (Off k = start[j]; k < start[j + 1]; ++k) {
                Idx r = index[k];
                if (r < 0 || r >= nrow) return "row index out of bounds";
                if (r <= prev) return "row indices not strictly ascending";
                if (value[k] == 0.0) return "explicit zero stored";
                prev = r;
            }
        }
        return "";
    }
};

struct CsrMatrix {
    Idx nrow = 0;
    Idx ncol = 0;
    std::vector<Off> start;
    std::vector<Idx> index;
    std::vector<double> value;

    CsrMatrix() : start(1, 0) {}
    CsrMatrix(Idx r, Idx c) : nrow(r), ncol(c), start(r + 1, 0) {}

    Off nnz() const { return start.empty() ? 0 : start[nrow]; }

    static CsrMatrix from_csc(const CscMatrix& A) {
        CsrMatrix R(A.nrow, A.ncol);
        std::vector<Off> row_counts(A.nrow, 0);
        for (Off k = 0; k < A.nnz(); ++k) {
            row_counts[A.index[k]]++;
        }
        R.start[0] = 0;
        for (Idx i = 0; i < A.nrow; ++i) {
            R.start[i + 1] = R.start[i] + row_counts[i];
        }
        R.index.resize(A.nnz());
        R.value.resize(A.nnz());
        std::vector<Off> write_pos = R.start;
        for (Idx j = 0; j < A.ncol; ++j) {
            for (Off k = A.start[j]; k < A.start[j + 1]; ++k) {
                Idx i = A.index[k];
                Off p = write_pos[i]++;
                R.index[p] = j;
                R.value[p] = A.value[k];
            }
        }
        return R;
    }

    CscMatrix to_csc() const {
        CscMatrix C(nrow, ncol);
        std::vector<Off> col_counts(ncol, 0);
        for (Off k = 0; k < nnz(); ++k) {
            col_counts[index[k]]++;
        }
        C.start[0] = 0;
        for (Idx j = 0; j < ncol; ++j) {
            C.start[j + 1] = C.start[j] + col_counts[j];
        }
        C.index.resize(nnz());
        C.value.resize(nnz());
        std::vector<Off> write_pos = C.start;
        for (Idx i = 0; i < nrow; ++i) {
            for (Off k = start[i]; k < start[i + 1]; ++k) {
                Idx j = index[k];
                Off p = write_pos[j]++;
                C.index[p] = i;
                C.value[p] = value[k];
            }
        }
        return C;
    }
};

struct MatrixPair {
    CscMatrix csc;
    CsrMatrix csr;

    void build(const CscMatrix& A) {
        csc = A;
        csr = CsrMatrix::from_csc(A);
    }
};

} // namespace sov
