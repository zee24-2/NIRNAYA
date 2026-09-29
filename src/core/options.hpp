#pragma once
#include <string>
#include <cmath>
#include "types.hpp"
#include "tolerances.hpp"

namespace sov {

struct Options {
    double time_limit         = 600.0;
    double work_limit         = kInf;
    int    threads            = 1;
    uint64_t seed             = 1;
    int    log_level          = 1;
    bool   deterministic      = true;
    int    memory_limit_mb    = 0;
    std::string method        = "auto"; // auto | dual | primal | ipm | ipm-crossover
    std::string preset        = "default"; // fast | default | aggressive | cautious
    std::string stats_json    = "";
    std::string sol_path      = "";
    std::string read_basis    = "";
    std::string write_basis   = "";

    Tolerances tol;

    int    presolve           = 2; // 0 off, 1 light, 2 full
    int    presolve_rounds    = 20;
    bool   scale              = true;
    int    refactor_interval  = 100;
    bool   bound_flipping     = true;
    bool   perturb            = true;
    Idx    simplex_iter_limit = 500000;
    int    ipm_max_iter       = 200;
    bool   qp_convexify       = false;

    std::string branching     = "reliability"; // reliability | pseudocost | mostfrac | strong
    Idx    node_limit         = 100000;
    Idx    solution_limit     = 0; // 0 = no limit
    bool   enable_cuts        = true;
    bool   enable_heuristics  = true;
    int    plunge_depth       = 10;
    int    strong_branch_iter = 200;
    int    strong_branch_rel  = 4;

    void apply_preset(const std::string& p) {
        preset = p;
        if (p == "cautious" || p == "strict") {
            tol.primal_feas_tol = 1e-8;
            tol.dual_feas_tol   = 1e-8;
            tol.lu_threshold_u  = 0.5;
            refactor_interval   = 30;
        } else if (p == "relaxed") {
            tol.primal_feas_tol = 1e-6;
            tol.dual_feas_tol   = 1e-6;
        } else if (p == "fast") {
            presolve_rounds     = 5;
            strong_branch_rel   = 2;
        }
    }

    bool set_by_name(const std::string& key, const std::string& val) {
        if (key == "time_limit")         { time_limit = std::stod(val); return true; }
        if (key == "work_limit")         { work_limit = std::stod(val); return true; }
        if (key == "threads")            { threads = std::stoi(val); return true; }
        if (key == "seed")               { seed = std::stoull(val); return true; }
        if (key == "log_level" || key == "log") { log_level = std::stoi(val); return true; }
        if (key == "deterministic")      { deterministic = (val == "1" || val == "true"); return true; }
        if (key == "method")             { method = val; return true; }
        if (key == "preset")             { apply_preset(val); return true; }
        if (key == "presolve")           { presolve = std::stoi(val); return true; }
        if (key == "scale")              { scale = (val != "0" && val != "false"); return true; }
        if (key == "refactor_interval")  { refactor_interval = std::stoi(val); return true; }
        if (key == "bound_flipping")     { bound_flipping = (val != "0" && val != "false"); return true; }
        if (key == "perturb")            { perturb = (val != "0" && val != "false"); return true; }
        if (key == "simplex_iter_limit") { simplex_iter_limit = std::stoi(val); return true; }
        if (key == "ipm_max_iter")       { ipm_max_iter = std::stoi(val); return true; }
        if (key == "mip_gap")            { tol.mip_gap = std::stod(val); return true; }
        if (key == "mip_abs_gap")        { tol.mip_abs_gap = std::stod(val); return true; }
        if (key == "primal_feas_tol")    { tol.primal_feas_tol = std::stod(val); return true; }
        if (key == "dual_feas_tol")      { tol.dual_feas_tol = std::stod(val); return true; }
        if (key == "int_tol")            { tol.int_tol = std::stod(val); return true; }
        if (key == "lu_threshold")       { tol.lu_threshold_u = std::stod(val); return true; }
        if (key == "node_limit")         { node_limit = std::stoi(val); return true; }
        if (key == "solution_limit")     { solution_limit = std::stoi(val); return true; }
        if (key == "branching")          { branching = val; return true; }
        if (key == "cuts")               { enable_cuts = (val != "0" && val != "false"); return true; }
        if (key == "heuristics")         { enable_heuristics = (val != "0" && val != "false"); return true; }
        return false;
    }
};

} // namespace sov
