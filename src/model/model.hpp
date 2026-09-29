#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <unordered_map>
#include "../core/types.hpp"
#include "../core/status.hpp"
#include "../core/tolerances.hpp"
#include "../core/csc_matrix.hpp"
#include "../core/dd_real.hpp"

namespace sov {

struct SosSet {
    int type = 1; // 1 = SOS1, 2 = SOS2
    std::string name;
    std::vector<Idx> indices;
    std::vector<double> weights;
};

struct Model {
    std::string name = "unnamed";
    Sense sense = Sense::Minimize;
    double obj_offset = 0.0;

    Idx ncol = 0;
    Idx nrow = 0;

    std::vector<double> cost;
    std::vector<double> collo;
    std::vector<double> colup;
    std::vector<VarType> vartype;

    CscMatrix A;
    std::vector<double> rowlo;
    std::vector<double> rowup;

    std::vector<std::string> col_names;
    std::vector<std::string> row_names;

    // Symmetric quadratic objective 0.5 x^T Q_full x stored as upper triangle (i <= j)
    CscMatrix Q;
    std::vector<SosSet> sos;

    bool is_mip() const {
        for (auto t : vartype) {
            if (t == VarType::Integer || t == VarType::Binary) return true;
        }
        return false;
    }

    bool is_qp() const {
        return Q.nnz() > 0;
    }

    Idx num_integers() const {
        Idx c = 0;
        for (auto t : vartype) {
            if (t == VarType::Integer || t == VarType::Binary) ++c;
        }
        return c;
    }

    Idx num_binaries() const {
        Idx c = 0;
        for (auto t : vartype) {
            if (t == VarType::Binary) ++c;
        }
        return c;
    }

    // Evaluate objective in original sense using double-double precision
    double eval_objective(const std::vector<double>& x) const {
        dd_real sum(obj_offset);
        for (Idx j = 0; j < ncol; ++j) {
            if (cost[j] != 0.0 && x[j] != 0.0) {
                sum.add_prod(cost[j], x[j]);
            }
        }
        if (Q.nnz() > 0) {
            dd_real qsum(0.0);
            for (Idx j = 0; j < Q.ncol; ++j) {
                double xj = x[j];
                if (xj == 0.0) continue;
                for (Off k = Q.start[j]; k < Q.start[j + 1]; ++k) {
                    Idx i = Q.index[k];
                    double qij = Q.value[k];
                    if (i == j) {
                        qsum.add_prod(0.5 * qij, xj * xj);
                    } else if (i < j) {
                        qsum.add_prod(qij, x[i] * xj);
                    }
                }
            }
            sum += qsum;
        }
        return sum.to_double();
    }

    // Convert to internal minimization model (C1)
    Model to_internal_min() const {
        Model m = *this;
        if (sense == Sense::Maximize) {
            m.sense = Sense::Minimize;
            m.obj_offset = -obj_offset;
            for (double& c : m.cost) c = -c;
            for (double& v : m.Q.value) v = -v;
        }
        return m;
    }

    double user_objective_from_internal(double internal_obj) const {
        return (sense == Sense::Maximize) ? -internal_obj : internal_obj;
    }

    // Validate model and repair minor bound issues (V1-V7)
    Status validate(std::string& message) {
        if (static_cast<Idx>(cost.size()) != ncol ||
            static_cast<Idx>(collo.size()) != ncol ||
            static_cast<Idx>(colup.size()) != ncol ||
            static_cast<Idx>(vartype.size()) != ncol) {
            message = "Column array dimension mismatch";
            return Status::InputError;
        }
        if (static_cast<Idx>(rowlo.size()) != nrow ||
            static_cast<Idx>(rowup.size()) != nrow) {
            message = "Row array dimension mismatch";
            return Status::InputError;
        }
        for (Idx j = 0; j < ncol; ++j) {
            if (!std::isfinite(cost[j])) {
                message = "Non-finite objective coefficient on column " + std::to_string(j);
                return Status::InputError;
            }
            collo[j] = sanitize_bound_lo(collo[j]);
            colup[j] = sanitize_bound_up(colup[j]);
            if (vartype[j] == VarType::Binary) {
                collo[j] = std::max(0.0, collo[j]);
                colup[j] = std::min(1.0, colup[j]);
            }
            if (vartype[j] == VarType::Integer || vartype[j] == VarType::Binary) {
                if (is_finite_bound(collo[j])) collo[j] = std::ceil(collo[j] - 1e-9);
                if (is_finite_bound(colup[j])) colup[j] = std::floor(colup[j] + 1e-9);
            }
            if (collo[j] > colup[j]) {
                if (collo[j] <= colup[j] + 1e-9 * (1.0 + std::abs(colup[j]))) {
                    collo[j] = colup[j];
                } else {
                    message = "Infeasible column bounds on variable " +
                              (j < static_cast<Idx>(col_names.size()) ? col_names[j] : std::to_string(j));
                    return Status::Infeasible;
                }
            }
            if (vartype[j] == VarType::Integer && collo[j] >= 0.0 && colup[j] <= 1.0) {
                vartype[j] = VarType::Binary;
            }
        }
        for (Idx i = 0; i < nrow; ++i) {
            rowlo[i] = sanitize_bound_lo(rowlo[i]);
            rowup[i] = sanitize_bound_up(rowup[i]);
            if (rowlo[i] > rowup[i]) {
                if (rowlo[i] <= rowup[i] + 1e-9 * (1.0 + std::abs(rowup[i]))) {
                    rowlo[i] = rowup[i];
                } else {
                    message = "Infeasible row bounds on row " +
                              (i < static_cast<Idx>(row_names.size()) ? row_names[i] : std::to_string(i));
                    return Status::Infeasible;
                }
            }
        }
        for (Off k = 0; k < A.nnz(); ++k) {
            if (!std::isfinite(A.value[k])) {
                message = "Non-finite matrix entry";
                return Status::InputError;
            }
        }
        return Status::Optimal; // Valid
    }
};

class ModelBuilder {
public:
    explicit ModelBuilder(const std::string& name = "model") {
        model_.name = name;
    }

    Idx add_col(double cost, double lo, double up, VarType type = VarType::Continuous, const std::string& name = "") {
        Idx j = model_.ncol++;
        model_.cost.push_back(cost);
        model_.collo.push_back(lo);
        model_.colup.push_back(up);
        model_.vartype.push_back(type);
        std::string nm = name.empty() ? ("C" + std::to_string(j)) : name;
        model_.col_names.push_back(nm);
        return j;
    }

    Idx add_row(double lo, double up, const std::vector<Idx>& indices, const std::vector<double>& values, const std::string& name = "") {
        Idx i = model_.nrow++;
        model_.rowlo.push_back(lo);
        model_.rowup.push_back(up);
        std::string nm = name.empty() ? ("R" + std::to_string(i)) : name;
        model_.row_names.push_back(nm);
        for (size_t k = 0; k < indices.size(); ++k) {
            if (values[k] != 0.0) {
                triplets_.push_back({i, indices[k], values[k]});
            }
        }
        return i;
    }

    void add_coeff(Idx row, Idx col, double val) {
        if (val != 0.0) {
            triplets_.push_back({row, col, val});
        }
    }

    void add_quad(Idx i, Idx j, double val) {
        if (val == 0.0) return;
        if (i > j) std::swap(i, j);
        q_triplets_.push_back({i, j, val});
    }

    void set_sense(Sense s) { model_.sense = s; }
    void set_offset(double c0) { model_.obj_offset = c0; }

    Model finalize() {
        model_.A = CscMatrix::from_triplets(model_.nrow, model_.ncol, triplets_);
        if (!q_triplets_.empty()) {
            model_.Q = CscMatrix::from_triplets(model_.ncol, model_.ncol, q_triplets_);
        } else {
            model_.Q = CscMatrix(model_.ncol, model_.ncol);
        }
        return model_;
    }

private:
    Model model_;
    std::vector<Triplet> triplets_;
    std::vector<Triplet> q_triplets_;
};

} // namespace sov
