#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cstdio>
#include "../model/model.hpp"
#include "../core/dd_real.hpp"

namespace sov {

struct VerificationReport {
    bool passed = true;
    double max_bound_viol = 0.0;
    double max_row_viol   = 0.0;
    double max_int_viol   = 0.0;
    double computed_obj   = 0.0;
    double obj_rel_err    = 0.0;
    std::string worst_item;
};

// Independent double-double checker routine (Part 12.6 & Part 13.1)
inline VerificationReport verify_solution(
    const Model& model,
    const std::vector<double>& x,
    double claimed_obj,
    double tol = 1e-6)
{
    VerificationReport rep;
    if (static_cast<Idx>(x.size()) != model.ncol) {
        rep.passed = false;
        rep.worst_item = "Solution size mismatch";
        return rep;
    }

    // 1. Check column bounds and integrality
    for (Idx j = 0; j < model.ncol; ++j) {
        double v = x[j];
        if (!std::isfinite(v)) {
            rep.passed = false;
            rep.worst_item = "NaN/Inf in variable " + model.col_names[j];
            return rep;
        }
        if (is_finite_bound(model.collo[j])) {
            double viol = std::max(0.0, model.collo[j] - v) / (1.0 + std::abs(model.collo[j]));
            if (viol > rep.max_bound_viol) {
                rep.max_bound_viol = viol;
                if (viol > tol) rep.worst_item = "Lower bound on " + model.col_names[j];
            }
        }
        if (is_finite_bound(model.colup[j])) {
            double viol = std::max(0.0, v - model.colup[j]) / (1.0 + std::abs(model.colup[j]));
            if (viol > rep.max_bound_viol) {
                rep.max_bound_viol = viol;
                if (viol > tol) rep.worst_item = "Upper bound on " + model.col_names[j];
            }
        }
        if (model.vartype[j] == VarType::Integer || model.vartype[j] == VarType::Binary) {
            double iv = std::abs(v - std::round(v));
            if (iv > rep.max_int_viol) {
                rep.max_int_viol = iv;
                if (iv > tol) rep.worst_item = "Integrality on " + model.col_names[j];
            }
        }
    }

    // 2. Check rows in dd_real (~32 decimal digits)
    std::vector<dd_real> row_act(model.nrow, dd_real(0.0));
    for (Idx j = 0; j < model.ncol; ++j) {
        double xj = x[j];
        if (xj == 0.0) continue;
        for (Off k = model.A.start[j]; k < model.A.start[j + 1]; ++k) {
            row_act[model.A.index[k]].add_prod(model.A.value[k], xj);
        }
    }
    for (Idx i = 0; i < model.nrow; ++i) {
        double r = row_act[i].to_double();
        if (is_finite_bound(model.rowlo[i])) {
            double viol = std::max(0.0, model.rowlo[i] - r) / (1.0 + std::abs(model.rowlo[i]));
            if (viol > rep.max_row_viol) {
                rep.max_row_viol = viol;
                if (viol > tol) rep.worst_item = "Row lower bound on " + model.row_names[i];
            }
        }
        if (is_finite_bound(model.rowup[i])) {
            double viol = std::max(0.0, r - model.rowup[i]) / (1.0 + std::abs(model.rowup[i]));
            if (viol > rep.max_row_viol) {
                rep.max_row_viol = viol;
                if (viol > tol) rep.worst_item = "Row upper bound on " + model.row_names[i];
            }
        }
    }

    // 3. Recompute objective in dd_real
    rep.computed_obj = model.eval_objective(x);
    if (std::isfinite(claimed_obj)) {
        rep.obj_rel_err = std::abs(rep.computed_obj - claimed_obj) / (1.0 + std::abs(rep.computed_obj));
    }

    if (rep.max_bound_viol > tol || rep.max_row_viol > tol || rep.max_int_viol > tol || rep.obj_rel_err > 1e-5) {
        rep.passed = false;
    }
    return rep;
}

class SolIo {
public:
    static bool write_sol(
        const std::string& path,
        const Model& model,
        Status status,
        double obj,
        double gap,
        const std::vector<double>& x)
    {
        std::ofstream out(path);
        if (!out.is_open()) return false;
        char buf[256];
        out << "# status " << status_to_string(status) << "\n";
        std::snprintf(buf, sizeof(buf), "# objective %.17g\n", obj);
        out << buf;
        std::snprintf(buf, sizeof(buf), "# gap %.6e\n", gap);
        out << buf;
        for (Idx j = 0; j < model.ncol && j < static_cast<Idx>(x.size()); ++j) {
            std::string nm = (j < static_cast<Idx>(model.col_names.size())) ? model.col_names[j] : ("C" + std::to_string(j));
            std::snprintf(buf, sizeof(buf), "%s %.17g\n", nm.c_str(), x[j]);
            out << buf;
        }
        return true;
    }

    static bool read_sol(
        const std::string& path,
        const Model& model,
        Status& out_status,
        double& out_obj,
        std::vector<double>& out_x)
    {
        std::ifstream in(path);
        if (!in.is_open()) return false;
        out_x.assign(model.ncol, 0.0);
        out_status = Status::Optimal;
        out_obj = 0.0;
        std::unordered_map<std::string, Idx> cmap;
        for (Idx j = 0; j < model.ncol; ++j) {
            cmap[model.col_names[j]] = j;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            if (line[0] == '#') {
                std::string hash, key, val;
                iss >> hash >> key >> val;
                if (key == "objective") out_obj = std::stod(val);
                continue;
            }
            std::string nm;
            double v = 0.0;
            if (iss >> nm >> v) {
                auto it = cmap.find(nm);
                if (it != cmap.end()) {
                    out_x[it->second] = v;
                }
            }
        }
        return true;
    }
};

} // namespace sov
