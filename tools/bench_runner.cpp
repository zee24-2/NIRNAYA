#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "../include/sov/solver.hpp"

// Benchmark & Robustness Showcase Runner (Slide 3 Verification bar & Slide 4 Yardstick)
int main() {
    std::printf("================================================================================\n");
    std::printf("  NIRNAYA (SIH26119) -- BENCHMARK & ROBUSTNESS SHOWCASE HARNESS\n");
    std::printf("================================================================================\n");

    // 1. Ill-conditioned LP Showcase (8-decade coefficient range)
    sov::ModelBuilder mb_ill("IllConditioned_8Decade_LP");
    auto x1 = mb_ill.add_col(1e-4, 0.0, 1e5, sov::VarType::Continuous, "x1");
    auto x2 = mb_ill.add_col(2.5,  0.0, 1e3, sov::VarType::Continuous, "x2");
    auto x3 = mb_ill.add_col(1e3,  0.0, 10.0, sov::VarType::Continuous, "x3");
    mb_ill.add_row(15.0, sov::kInf, {x1, x2, x3}, {1e-4, 1.0, 1e2}, "R1");
    mb_ill.add_row(25.0, sov::kInf, {x1, x2, x3}, {3e-4, 2.0, 50.0}, "R2");
    auto m_ill = mb_ill.finalize();

    // 2. Degenerate Transportation / Assignment LP
    sov::ModelBuilder mb_deg("Degenerate_Transport_LP");
    const int N = 5;
    std::vector<std::vector<sov::Idx>> flow(N, std::vector<sov::Idx>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double c = static_cast<double>((i + 1) * (j + 2) % 7 + 1);
            flow[i][j] = mb_deg.add_col(c, 0.0, 1.0, sov::VarType::Continuous);
        }
    }
    for (int i = 0; i < N; ++i) {
        mb_deg.add_row(1.0, 1.0, flow[i], {1.0, 1.0, 1.0, 1.0, 1.0});
        std::vector<sov::Idx> col_vars;
        for (int r = 0; r < N; ++r) col_vars.push_back(flow[r][i]);
        mb_deg.add_row(1.0, 1.0, col_vars, {1.0, 1.0, 1.0, 1.0, 1.0});
    }
    auto m_deg = mb_deg.finalize();

    // 3. Root-Gap Closure Showcase (Multi-Constraint Fractional MILP: No Cuts vs With Cuts)
    sov::ModelBuilder mb_cut("WeakRelaxation_Cover_MILP");
    mb_cut.set_sense(sov::Sense::Minimize);
    std::vector<sov::Idx> bvars;
    for (int j = 0; j < 6; ++j) {
        bvars.push_back(mb_cut.add_col(10.0 + 3.0 * j, 0.0, 1.0, sov::VarType::Binary));
    }
    auto c0 = mb_cut.add_col(1.5, 0.0, 20.0, sov::VarType::Continuous);
    mb_cut.add_row(18.0, sov::kInf, {bvars[0], bvars[1], bvars[2], bvars[3], bvars[4], bvars[5], c0},
                   {7.0, 8.0, 6.0, 5.0, 7.0, 4.0, 0.2});
    mb_cut.add_row(-sov::kInf, 11.0, {bvars[0], bvars[1], bvars[2], bvars[4]}, {4.0, 5.0, 4.0, 5.0});
    mb_cut.add_row(1.0, sov::kInf, {bvars[2], bvars[3], bvars[5]}, {1.0, 1.0, 1.0});
    auto m_cut = mb_cut.finalize();

    sov::Options opt_full;
    opt_full.log_level = 0;

    sov::Options opt_nocuts = opt_full;
    opt_nocuts.enable_cuts = false;

    sov::Options opt_pdhg = opt_full;
    opt_pdhg.method = "pdhg";

    auto r_ill  = sov::Solver(opt_full).solve(m_ill);
    auto r_deg  = sov::Solver(opt_full).solve(m_deg);
    auto r_pdhg = sov::Solver(opt_pdhg).solve(m_deg);
    auto r_cut_off = sov::Solver(opt_nocuts).solve(m_cut);
    auto r_cut_on  = sov::Solver(opt_full).solve(m_cut);

    std::printf("  %-30s | %-10s | %-12s | %-8s | %-8s | %-8s\n",
                "Benchmark / Showcase", "Method", "Objective", "Nodes", "Cuts", "Verified");
    std::printf("  ------------------------------------------------------------------------------\n");
    std::printf("  %-30s | %-10s | %12.6g | %8d | %8d | %s\n",
                "IllConditioned_8Decade_LP", "Scale+Dual", r_ill.objective, r_ill.nodes, r_ill.cuts_added,
                r_ill.verification.passed ? "PASS" : "FAIL");
    std::printf("  %-30s | %-10s | %12.6g | %8d | %8d | %s\n",
                "Degenerate_Transport_LP", "Dual+DSE", r_deg.objective, r_deg.nodes, r_deg.cuts_added,
                r_deg.verification.passed ? "PASS" : "FAIL");
    std::printf("  %-30s | %-10s | %12.6g | %8d | %8d | %s\n",
                "Degenerate_Transport_LP", "PDHG (GPU)", r_pdhg.objective, r_pdhg.nodes, r_pdhg.cuts_added,
                r_pdhg.verification.passed ? "PASS" : "FAIL");
    std::printf("  %-30s | %-10s | %12.6g | %8d | %8d | %s\n",
                "WeakRelaxation_MILP (No Cuts)", "B&B Only", r_cut_off.objective, r_cut_off.nodes, r_cut_off.cuts_added,
                r_cut_off.verification.passed ? "PASS" : "FAIL");
    std::printf("  %-30s | %-10s | %12.6g | %8d | %8d | %s\n",
                "WeakRelaxation_MILP (With Cuts)", "Branch&Cut", r_cut_on.objective, r_cut_on.nodes, r_cut_on.cuts_added,
                r_cut_on.verification.passed ? "PASS" : "FAIL");
    std::printf("================================================================================\n");
    return 0;
}
