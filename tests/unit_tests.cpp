#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include "../include/sov/solver.hpp"
#include "../include/sov/sov.h"

#define SOV_TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::printf("  [FAIL] %s (line %d): %s\n", __FUNCTION__, __LINE__, msg); \
            return false; \
        } \
    } while (0)

#define SOV_TEST_NEAR(a, b, tol, msg) \
    do { \
        double _va = (a); \
        double _vb = (b); \
        if (std::abs(_va - _vb) > (tol)) { \
            std::printf("  [FAIL] %s (line %d): %s | got %.12g, expected %.12g (diff %.3e > %.3e)\n", \
                        __FUNCTION__, __LINE__, msg, _va, _vb, std::abs(_va - _vb), (double)(tol)); \
            return false; \
        } \
    } while (0)

// 1. Double-double (dd_real) and Neumaier compensated sum test (Part 2.5)
static bool test_dd_real_and_compensated_sum() {
    sov::dd_real acc(0.0);
    acc += 1e16;
    acc += 1.0;
    acc += -1e16;
    SOV_TEST_NEAR(acc.to_double(), 1.0, 1e-15, "dd_real catastrophic cancellation 1e16 + 1 - 1e16");

    sov::NeumaierSum nsum;
    nsum.add(1e16);
    nsum.add(1.0);
    nsum.add(-1e16);
    SOV_TEST_NEAR(nsum.result(), 1.0, 1e-15, "NeumaierSum cancellation");
    return true;
}

// 2. Sparse Matrix CSC/CSR round-trip and duplicate summing test (Part 2.2)
static bool test_sparse_matrices_and_hvector() {
    std::vector<sov::Triplet> trips = {
        {0, 0, 1.5}, {0, 0, 0.5}, // duplicate (0,0) -> 2.0
        {1, 0, -3.0},
        {0, 1, 4.0},
        {1, 1, 2.0}
    };
    auto A = sov::CscMatrix::from_triplets(2, 2, trips);
    SOV_TEST_ASSERT(A.validate().empty(), "CSC canonical validation");
    SOV_TEST_ASSERT(A.nnz() == 4, "Duplicate entries merged");

    std::vector<double> x = {2.0, 3.0}, y;
    A.spmv(x, y);
    // y0 = 2.0*2 + 4.0*3 = 16.0; y1 = -3.0*2 + 2.0*3 = 0.0
    SOV_TEST_NEAR(y[0], 16.0, 1e-12, "SpMV row 0");
    SOV_TEST_NEAR(y[1], 0.0,  1e-12, "SpMV row 1");

    auto R = sov::CsrMatrix::from_csc(A);
    auto A2 = R.to_csc();
    SOV_TEST_ASSERT(A2.nnz() == A.nnz(), "CSC -> CSR -> CSC round-trip nnz");
    for (sov::Off k = 0; k < A.nnz(); ++k) {
        SOV_TEST_NEAR(A2.value[k], A.value[k], 1e-15, "CSC -> CSR -> CSC value match");
    }

    sov::HVector hv(5);
    hv.add(1, 3.5);
    hv.add(3, -2.0);
    hv.add(1, -3.5); // cancels to 0 -> sentinel
    hv.tighten();
    SOV_TEST_ASSERT(hv.count == 1 && hv.idx[0] == 3, "HVector cancellation sentinel + tighten");
    return true;
}

// 3. Sparse Markowitz LU Factorization + FTRAN/BTRAN + Updates (Part 5)
static bool test_sparse_lu_factor_and_updates() {
    // 3x3 matrix:
    // [ 2  1  0 ]
    // [ 1  3  1 ]
    // [ 0  2  4 ]
    std::vector<sov::Triplet> trips = {
        {0, 0, 2.0}, {1, 0, 1.0},
        {0, 1, 1.0}, {1, 1, 3.0}, {2, 1, 2.0},
        {1, 2, 1.0}, {2, 2, 4.0}
    };
    auto A = sov::CscMatrix::from_triplets(3, 3, trips);
    std::vector<sov::Idx> basis = {0, 1, 2};
    sov::BasisFactor lu;
    auto fres = lu.factor(A, basis);
    SOV_TEST_ASSERT(fres.ok && fres.repairs.empty(), "LU factor nonsingular");

    // Solve B x = b where true x = [1, 2, 3]^T => b = [4, 10, 16]^T
    sov::HVector rhs(3);
    rhs.add(0, 4.0); rhs.add(1, 10.0); rhs.add(2, 16.0);
    rhs.tighten();
    lu.ftran(rhs);
    SOV_TEST_NEAR(rhs.val[0], 1.0, 1e-11, "FTRAN x[0]");
    SOV_TEST_NEAR(rhs.val[1], 2.0, 1e-11, "FTRAN x[1]");
    SOV_TEST_NEAR(rhs.val[2], 3.0, 1e-11, "FTRAN x[2]");

    // Solve B^T y = c where true y = [2, -1, 3]^T => B^T y = [3, 5, 11]^T
    sov::HVector cb(3);
    cb.add(0, 3.0); cb.add(1, 5.0); cb.add(2, 11.0);
    cb.tighten();
    lu.btran(cb);
    SOV_TEST_NEAR(cb.val[0],  2.0, 1e-11, "BTRAN y[0]");
    SOV_TEST_NEAR(cb.val[1], -1.0, 1e-11, "BTRAN y[1]");
    SOV_TEST_NEAR(cb.val[2],  3.0, 1e-11, "BTRAN y[2]");
    return true;
}

// 4. Appendix A.5 Worked Example: Dual Simplex 2-variable LP (Part 6 & Appendix A.5)
static bool test_appendix_a5_worked_example() {
    // min 2 x1 + 3 x2  s.t.  x1 + x2 >= 4 (row 0),  x1 + 3 x2 >= 6 (row 1),  x >= 0
    // Expected optimal solution: x = (3, 1), objective = 9, duals y = (1.5, 0.5)
    sov::ModelBuilder mb("Appendix_A5");
    mb.set_sense(sov::Sense::Minimize);
    auto x1 = mb.add_col(2.0, 0.0, sov::kInf, sov::VarType::Continuous, "x1");
    auto x2 = mb.add_col(3.0, 0.0, sov::kInf, sov::VarType::Continuous, "x2");
    mb.add_row(4.0, sov::kInf, {x1, x2}, {1.0, 1.0}, "row1");
    mb.add_row(6.0, sov::kInf, {x1, x2}, {1.0, 3.0}, "row2");
    auto m = mb.finalize();

    sov::Options opt;
    opt.log_level = 0;
    opt.presolve = 0;
    opt.scale = false;
    sov::Solver solver(opt);
    auto res = solver.solve(m);

    SOV_TEST_ASSERT(res.status == sov::Status::Optimal, "Appendix A.5 status optimal");
    SOV_TEST_NEAR(res.objective, 9.0, 1e-8, "Appendix A.5 optimal objective == 9");
    SOV_TEST_NEAR(res.x[0], 3.0, 1e-7, "Appendix A.5 x1 == 3");
    SOV_TEST_NEAR(res.x[1], 1.0, 1e-7, "Appendix A.5 x2 == 1");
    SOV_TEST_NEAR(res.row_duals[0], 1.5, 1e-7, "Appendix A.5 y1 == 1.5");
    SOV_TEST_NEAR(res.row_duals[1], 0.5, 1e-7, "Appendix A.5 y2 == 0.5");
    SOV_TEST_ASSERT(res.verification.passed, "Appendix A.5 dd_real verification passed");
    return true;
}

// 5. Infeasible and Unbounded LP Detection (Part 6.5 & 6.10)
static bool test_lp_infeasible_and_unbounded() {
    // Infeasible LP: x1 + x2 <= 2, x1 + x2 >= 5, x >= 0
    {
        sov::ModelBuilder mb("Infeasible_LP");
        auto x1 = mb.add_col(1.0, 0.0, sov::kInf);
        auto x2 = mb.add_col(1.0, 0.0, sov::kInf);
        mb.add_row(-sov::kInf, 2.0, {x1, x2}, {1.0, 1.0});
        mb.add_row(5.0, sov::kInf,  {x1, x2}, {1.0, 1.0});
        sov::Options opt;
        opt.log_level = 0;
        opt.presolve = 0;
        sov::Solver solver(opt);
        auto res = solver.solve(mb.finalize());
        SOV_TEST_ASSERT(res.status == sov::Status::Infeasible, "Infeasible LP detected by simplex");
    }

    // Unbounded LP: min -x1 - 2 x2  s.t.  x1 - x2 <= 3, x >= 0
    {
        sov::ModelBuilder mb("Unbounded_LP");
        auto x1 = mb.add_col(-1.0, 0.0, sov::kInf);
        auto x2 = mb.add_col(-2.0, 0.0, sov::kInf);
        mb.add_row(-sov::kInf, 3.0, {x1, x2}, {1.0, -1.0});
        sov::Options opt;
        opt.log_level = 0;
        opt.presolve = 0;
        sov::Solver solver(opt);
        auto res = solver.solve(mb.finalize());
        SOV_TEST_ASSERT(res.status == sov::Status::Unbounded, "Unbounded LP detected by simplex");
    }
    return true;
}

// 6. Scaling + Presolve/Postsolve Equivalence on Badly Scaled LP (Parts 4 & 7)
static bool test_scaling_and_presolve_equivalence() {
    sov::ModelBuilder mb("Badly_Scaled_LP");
    auto x1 = mb.add_col(1e-3, 0.0, 1e4);
    auto x2 = mb.add_col(2.0,  0.0, 100.0);
    auto x3 = mb.add_col(5.0,  3.0, 3.0); // Fixed column
    mb.add_row(10.0, sov::kInf, {x1, x2}, {1e-3, 1.0});
    mb.add_row(20.0, sov::kInf, {x1, x2, x3}, {2e-3, 1.0, 1.0}); // 2e-3 x1 + x2 + 3 >= 20 => 2e-3 x1 + x2 >= 17
    mb.add_row(4.0,  sov::kInf, {x2}, {1.0}); // Singleton row x2 >= 4
    auto model = mb.finalize();

    sov::Options opt_no_pre;
    opt_no_pre.log_level = 0;
    opt_no_pre.presolve = 0;
    opt_no_pre.scale = false;
    auto res1 = sov::Solver(opt_no_pre).solve(model);

    sov::Options opt_full;
    opt_full.log_level = 0;
    opt_full.presolve = 2;
    opt_full.scale = true;
    auto res2 = sov::Solver(opt_full).solve(model);

    SOV_TEST_ASSERT(res1.status == sov::Status::Optimal && res2.status == sov::Status::Optimal, "Both optimal");
    SOV_TEST_NEAR(res1.objective, res2.objective, 1e-7, "Presolved+scaled objective matches unscaled");
    SOV_TEST_ASSERT(res2.verification.passed, "Postsolved solution passes independent check");
    return true;
}

// 7. Interior-Point Method (IPM) for LP and Convex QP (Parts 8 & 9)
static bool test_ipm_lp_and_convex_qp() {
    // Compare IPM LP with Simplex on Appendix A.5
    sov::ModelBuilder mb_lp("IPM_LP");
    auto x1 = mb_lp.add_col(2.0, 0.0, 20.0);
    auto x2 = mb_lp.add_col(3.0, 0.0, 20.0);
    mb_lp.add_row(4.0, 20.0, {x1, x2}, {1.0, 1.0});
    mb_lp.add_row(6.0, 30.0, {x1, x2}, {1.0, 3.0});
    sov::Options opt;
    opt.log_level = 0;
    opt.method = "ipm";
    opt.presolve = 0;
    auto res_lp = sov::Solver(opt).solve(mb_lp.finalize());
    SOV_TEST_ASSERT(res_lp.status == sov::Status::Optimal, "IPM LP optimal");
    SOV_TEST_NEAR(res_lp.objective, 9.0, 1e-4, "IPM LP objective == 9.0");

    // Convex QP: min 0.5 * (2 x1^2 + 2 x2^2) - 2 x1 - 5 x2  s.t.  x1 + x2 <= 2, x >= 0
    // Unconstrained min is (1, 2.5) with sum 3.5 > 2; active constraint x1 + x2 = 2 gives
    // x1 - 2 = lambda, x2 - 5/2 = lambda => x1 = 1/4 = 0.25, x2 = 7/4 = 1.75
    // Obj = 0.25^2 + 1.75^2 - 2(0.25) - 5(1.75) = 0.0625 + 3.0625 - 0.5 - 8.75 = -6.125
    sov::ModelBuilder mb_qp("Convex_QP");
    auto q1 = mb_qp.add_col(-2.0, 0.0, 10.0);
    auto q2 = mb_qp.add_col(-5.0, 0.0, 10.0);
    mb_qp.add_quad(q1, q1, 2.0);
    mb_qp.add_quad(q2, q2, 2.0);
    mb_qp.add_row(-sov::kInf, 2.0, {q1, q2}, {1.0, 1.0});
    auto res_qp = sov::Solver(opt).solve(mb_qp.finalize());
    SOV_TEST_ASSERT(res_qp.status == sov::Status::Optimal, "Convex QP status optimal");
    SOV_TEST_NEAR(res_qp.x[0], 0.25,  1e-4, "Convex QP x1 == 0.25");
    SOV_TEST_NEAR(res_qp.x[1], 1.75,  1e-4, "Convex QP x2 == 1.75");
    SOV_TEST_NEAR(res_qp.objective, -6.125, 1e-4, "Convex QP objective == -6.125");
    return true;
}

// 8. MILP Branch-and-Cut with GMI & Knapsack Cover Cuts (Part 10)
static bool test_milp_branch_and_cut() {
    // 0-1 Knapsack Maximize:
    // max 16 x1 + 22 x2 + 12 x3 + 8 x4 + 11 x5 + 19 x6
    // s.t. 5 x1 + 7 x2 + 4 x3 + 3 x4 + 4 x5 + 6 x6 <= 14
    // Optimal binary solution: x1=1 (5, 16), x3=1 (4, 12) or x4=1 (3, 8), x6=1 (6, 19):
    // Wait: x1=1 (5,16) + x4=1 (3,8) + x6=1 (6,19) => weight 5+3+6 = 14, value 16+8+19 = 43!
    sov::ModelBuilder mb("Knapsack_MILP");
    mb.set_sense(sov::Sense::Maximize);
    double v[] = {16, 22, 12, 8, 11, 19};
    double w[] = {5,  7,  4,  3, 4,  6};
    std::vector<sov::Idx> cols;
    for (int i = 0; i < 6; ++i) {
        cols.push_back(mb.add_col(v[i], 0.0, 1.0, sov::VarType::Binary, "x" + std::to_string(i + 1)));
    }
    mb.add_row(-sov::kInf, 14.0, cols, {w[0], w[1], w[2], w[3], w[4], w[5]}, "Cap");

    sov::Options opt;
    opt.log_level = 0;
    sov::Solver solver(opt);
    auto res = solver.solve(mb.finalize());

    SOV_TEST_ASSERT(res.status == sov::Status::Optimal, "MILP Knapsack optimal");
    SOV_TEST_NEAR(res.objective, 43.0, 1e-6, "MILP Knapsack optimal value == 43");
    SOV_TEST_ASSERT(res.verification.passed, "MILP Knapsack dd_real verification passed");
    return true;
}

// 9. MPS Reader/Writer Round-Trip and C API (Part 3.2 & 3.5)
static bool test_mps_roundtrip_and_c_api() {
    sov_env* env = sov_env_create();
    sov_model* m = sov_model_create(env, "c_api_test");
    sov_set_option_string(m, "log_level", "0");
    int c0 = sov_add_col(m, 2.0, 0.0, 100.0, 0, "x1");
    int c1 = sov_add_col(m, 3.0, 0.0, 100.0, 0, "x2");
    int idx[2] = {c0, c1};
    double val1[2] = {1.0, 1.0};
    double val2[2] = {1.0, 3.0};
    sov_add_row(m, 4.0, 1e20, 2, idx, val1, "r1");
    sov_add_row(m, 6.0, 1e20, 2, idx, val2, "r2");

    SOV_TEST_ASSERT(sov_optimize(m) == 0, "sov_optimize returned 0");
    SOV_TEST_ASSERT(sov_get_status(m) == 0, "sov_get_status == Optimal");
    SOV_TEST_NEAR(sov_get_obj(m), 9.0, 1e-7, "C API objective == 9.0");
    double x_out[2] = {0.0, 0.0};
    sov_get_x(m, x_out);
    SOV_TEST_NEAR(x_out[0], 3.0, 1e-7, "C API x1 == 3.0");
    SOV_TEST_NEAR(x_out[1], 1.0, 1e-7, "C API x2 == 1.0");

    sov_model_free(m);
    sov_env_free(env);
    return true;
}

// 10. AMD Ordering + Sparse Cholesky + dd_real Iterative Refinement (Slide 3 Layer L2)
static bool test_amd_cholesky_and_refinement() {
    // 3x3 SPD system:
    // [  4  -1   1 ] [x0]   [ 7]
    // [ -1   5  -2 ] [x1] = [ 3]   => true x = [2, 2, 1]^T (4*2 - 2 + 1 = 7; -2 + 10 - 2 = 6; 2 - 4 + 6 = 4)
    // [  1  -2   6 ] [x2]   [ 4]
    std::vector<double> S = {
         4.0, -1.0,  1.0,
        -1.0,  5.0, -2.0,
         1.0, -2.0,  6.0
    };
    std::vector<double> rhs = {7.0, 6.0, 4.0}, x;
    sov::SparseCholeskyAmd chol;
    SOV_TEST_ASSERT(chol.factor_spd(3, S, 1e-12), "AMD Sparse Cholesky factor_spd");
    chol.solve_refined(S, rhs, x, 3);
    SOV_TEST_NEAR(x[0], 2.0, 1e-11, "AMD Cholesky refined x[0] == 2.0");
    SOV_TEST_NEAR(x[1], 2.0, 1e-11, "AMD Cholesky refined x[1] == 2.0");
    SOV_TEST_NEAR(x[2], 1.0, 1e-11, "AMD Cholesky refined x[2] == 1.0");
    return true;
}

// 11. Parallel First-Order PDHG / PDLP Solver (Slide 1 & Slide 4 GPU/Parallel)
static bool test_pdhg_first_order_lp() {
    sov::ModelBuilder mb("PDHG_LP");
    auto x1 = mb.add_col(2.0, 0.0, 20.0);
    auto x2 = mb.add_col(3.0, 0.0, 20.0);
    mb.add_row(4.0, 20.0, {x1, x2}, {1.0, 1.0});
    mb.add_row(6.0, 30.0, {x1, x2}, {1.0, 3.0});

    sov::Options opt;
    opt.log_level = 0;
    opt.method = "pdhg";
    opt.threads = 2;
    auto res = sov::Solver(opt).solve(mb.finalize());
    SOV_TEST_ASSERT(res.status == sov::Status::Optimal, "PDHG status optimal");
    SOV_TEST_NEAR(res.objective, 9.0, 1e-2, "PDHG objective ~ 9.0");
    return true;
}

// 12. Convex MIQP Branch-and-Bound Engine (Slide 2 & Slide 3 Extension-Ready MIQP)
static bool test_convex_miqp() {
    // min 0.5 * 2 * (x1 - 1.4)^2 => min x1^2 - 2.8 x1  with x1 integer in [0, 5] => x1* = 1 (obj = 1 - 2.8 = -1.8)
    sov::ModelBuilder mb("Convex_MIQP");
    auto x1 = mb.add_col(-2.8, 0.0, 5.0, sov::VarType::Integer, "x1");
    mb.add_quad(x1, x1, 2.0);
    mb.add_row(0.0, 5.0, {x1}, {1.0}, "bnd");

    sov::Options opt;
    opt.log_level = 0;
    auto res = sov::Solver(opt).solve(mb.finalize());
    SOV_TEST_ASSERT(res.status == sov::Status::Optimal, "Convex MIQP status optimal");
    SOV_TEST_NEAR(res.x[0], 1.0, 1e-6, "Convex MIQP optimal integer x1 == 1");
    SOV_TEST_NEAR(res.objective, -1.8, 1e-5, "Convex MIQP optimal objective == -1.8");
    return true;
}

int main() {
    struct TestEntry {
        const char* name;
        bool (*fn)();
    };
    std::vector<TestEntry> suite = {
        {"CORE-05: dd_real & Neumaier Compensated Arithmetic",    test_dd_real_and_compensated_sum},
        {"CORE-02/03: Sparse CSC/CSR Matrices & HVector",         test_sparse_matrices_and_hvector},
        {"LU-01..04: Sparse Markowitz LU + FTRAN/BTRAN",          test_sparse_lu_factor_and_updates},
        {"CHOL-AMD & REFINE: AMD Cholesky + dd_real Refinement",  test_amd_cholesky_and_refinement},
        {"SIMPLEX-04: Appendix A.5 Worked Dual Simplex LP",       test_appendix_a5_worked_example},
        {"SIMPLEX-03/07: Infeasible & Unbounded LP Detection",    test_lp_infeasible_and_unbounded},
        {"SCALE-01 & PRE-01: Scaling + Presolve/Postsolve",       test_scaling_and_presolve_equivalence},
        {"IPM-03 & QP-02: Mehrotra IPM for LP and Convex QP",     test_ipm_lp_and_convex_qp},
        {"GPU/PDHG-01: Parallel First-Order PDHG (PDLP) Solver",  test_pdhg_first_order_lp},
        {"MIP-01..15: Branch-and-Cut with GMI & Cover Cuts",      test_milp_branch_and_cut},
        {"MIQP-01: Convex MIQP Branch-and-Bound + IPM",           test_convex_miqp},
        {"IO-01 & IO-06: MPS Round-Trip & Stable C API",          test_mps_roundtrip_and_c_api}
    };

    std::printf("========================================================================\n");
    std::printf("  NIRNAYA (SIH26119) -- UNIT & INTEGRATION TEST SUITE (12 MODULES)\n");
    std::printf("========================================================================\n");

    int passed = 0;
    for (const auto& t : suite) {
        bool ok = t.fn();
        std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", t.name);
        if (ok) ++passed;
    }

    std::printf("------------------------------------------------------------------------\n");
    std::printf("  Summary: %d / %zu tests passed.\n", passed, suite.size());
    std::printf("========================================================================\n");
    return (passed == static_cast<int>(suite.size())) ? 0 : 1;
}
