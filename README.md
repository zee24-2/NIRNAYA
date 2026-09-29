# NIRNAYA — Sovereign Optimization Solver Core (`SIH26119`)

**Smart India Hackathon 2026 | Problem Statement ID: `SIH26119` | Team ID: `151198` (*Visionaries for Change*)**  
**Problem Owner:** Mangalore Refinery and Petrochemicals Limited (**MRPL**)  
**Theme:** Smart Automation — *High-Performance Mathematical Optimization Solver for LP, MILP, Convex QP & MIQP*

---

## 1. Executive Summary & PDF Commitment Coherence

**NIRNAYA** (*Sovereign Optimization Solver Core*) is a clean-room, zero-external-dependency mathematical optimization engine engineered in **C++17** from first principles. Every commitment highlighted in the 6-page **SIH26119 Presentation PDF** is implemented, tested, and verified:

| PDF Presentation Slide & Commitment | NIRNAYA Implementation & Verification |
| :--- | :--- |
| **Slide 1: Title & Sovereign Vision**<br>Problem `SIH26119`, Team `151198` (*Visionaries for Change*), GPU-Accelerated Optimization Solver for LP & MILP | Full **NIRNAYA** CLI (`bin/nirnaya.exe` & `bin/sovsolve.exe`), shared library (`bin/nirnaya.dll`), Python ctypes SDK (`python/nirnaya.py`), and Parallel First-Order **PDHG (Chambolle-Pock PDLP)** solver (`src/gpu/pdhg_solver.hpp`). |
| **Slide 2: All 5 Industrial Domains + MIQP**<br>1. *Refinery scheduling*<br>2. *Crude blending*<br>3. *Production planning*<br>4. *Power dispatch*<br>5. *Logistics & supply chain*<br>+ *Convex QP/MIQP Extension-Ready* | All **6 industrial models** implemented in `tools/gen_refinery_models.cpp`, exported to `data/refinery/*.mps`, solved to proven optimality (`0.0000%` gap), and verified by `bin/checker.exe` in ~32-digit `dd_real` arithmetic. |
| **Slide 3: 6-Layer Sovereign Architecture**<br>- **L6 Interfaces**: CLI, C API, Python (`ctypes`), MPS/LP<br>- **L5 Driver & Strategy**: Problem classification, auto-dispatch<br>- **L4 MIP Engine**: Presolve, B&B tree, GMI/Cover cuts, Heuristics<br>- **L3 Continuous**: Dual/Primal Simplex, Mehrotra IPM, PDHG<br>- **L2 Linear Algebra**: Sparse Markowitz LU, AMD Cholesky, `dd_real` refinement<br>- **L1 Core**: `CSC`/`CSR`, `HVector`, `ThreadPool`, deterministic RNG | 100% faithful 6-layer header/source hierarchy under `src/` and `include/sov/`, with **Sparse AMD $LDL^T$ Cholesky + `dd_real` Iterative Refinement** (`src/linalg/sparse_cholesky_amd.hpp`) and **Convex MIQP Branch-and-Bound**. |
| **Slide 4: Key Design Highlights & Differentiators**<br>- Power-of-2 IEEE-754 exponent scaling<br>- ~32-digit `dd_real` extended precision & independent checker<br>- Deterministic Parallel `ThreadPool`<br>- Multi-Language Embeddability (`sov.h`, `nirnaya.dll`, `nirnaya.py`) | Verified across **12 unit & integration test suites** (`bin/unit_tests.exe`: `12 / 12 PASS`) and the **Benchmark & Robustness Showcase Harness** (`bin/bench_runner.exe`). |

---

## 2. Six-Layer Architecture Map

```
+-----------------------------------------------------------------------------------+
| L6: INTERFACES & VERIFICATION                                                     |
|   - CLI: bin/nirnaya.exe & bin/sovsolve.exe (--method=auto|simplex|ipm|pdhg)      |
|   - C ABI & Shared DLL: include/sov/sov.h, src/driver/c_api.cpp, bin/nirnaya.dll  |
|   - Python Bindings: python/nirnaya.py (Zero-dependency ctypes wrapper)           |
|   - Independent Verifier: tools/checker.cpp (~32-digit dd_real certificate check) |
+-----------------------------------------------------------------------------------+
| L5: DRIVER, STRATEGY & CONVEX MIQP ENGINE                                         |
|   - Unified Facade: include/sov/solver.hpp (LP / MILP / Convex QP / Convex MIQP)  |
+-----------------------------------------------------------------------------------+
| L4: MIXED-INTEGER PROGRAMMING (MIP) & PRESOLVE ENGINE                             |
|   - Presolve/Postsolve: src/presolve/presolve.hpp, src/model/scaling.hpp          |
|   - Branch-and-Cut: src/mip/mip_solver.hpp (Reliability branching, GMI & Cover)   |
+-----------------------------------------------------------------------------------+
| L3: CONTINUOUS OPTIMIZATION ENGINES                                               |
|   - Revised Dual & Primal Simplex: src/lp/lp_solver.hpp (DSE, Harris ratio test)  |
|   - Mehrotra Predictor-Corrector IPM: src/ipm/ipm_lp_qp.hpp (LP & Convex QP)      |
|   - Parallel First-Order PDHG (PDLP): src/gpu/pdhg_solver.hpp (SpMV/SpMTV)        |
+-----------------------------------------------------------------------------------+
| L2: SPARSE NUMERICAL LINEAR ALGEBRA                                               |
|   - Sparse Markowitz LU + Eta Updates: src/linalg/lu_factor.hpp                   |
|   - Sparse AMD Cholesky + dd_real Refinement: src/linalg/sparse_cholesky_amd.hpp  |
+-----------------------------------------------------------------------------------+
| L1: CORE NUMERICAL & PARALLEL PRIMITIVES                                          |
|   - dd_real (~32-digit arithmetic) & Neumaier Summation: src/core/dd_real.hpp     |
|   - Sparse CSC/CSR & HVector: src/core/csc_matrix.hpp, src/core/hvector.hpp       |
|   - Deterministic Work-Stealing ThreadPool: src/core/thread_pool.hpp              |
+-----------------------------------------------------------------------------------+
```

---

## 3. Build & Verification Commands

Compile and run all executables using `g++` (C++17, `-fno-fast-math` to preserve IEEE-754 extended precision error-free transformations):

```powershell
# Run 12-module unit & integration test suite
.\bin\unit_tests.exe

# Generate & solve all 6 MRPL Industrial Models (LP, MILP, Convex QP, Convex MIQP)
.\bin\gen_refinery_models.exe

# Run Benchmark & Robustness Showcase (8-decade scaling, degeneracy, PDHG, root-gap cuts)
.\bin\bench_runner.exe

# Solve any MPS model via the NIRNAYA CLI (methods: auto, simplex, ipm, pdhg)
.\bin\nirnaya.exe data/refinery/MRPL_Crude_Scheduling_MILP.mps --sol=sched.sol

# Independently verify any solution in ~32-digit dd_real arithmetic
.\bin\checker.exe data/refinery/MRPL_UnitCommitment_Dispatch_MIQP.mps data/refinery/MRPL_UnitCommitment_Dispatch_MIQP.sol
```

---

## 4. Industrial Showcase Results (`gen_refinery_models.exe`)

| Industrial Domain (Slide 2) | Model Name | Class | Rows × Cols | Optimal Objective | `dd_real` Max Row Viol | Status |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| **1. Crude Blending** | `MRPL_Crude_Blending_LP` | LP | `6 × 6` | `22171.5977742` | `2.55e-14` | **PASS** |
| **2. Production Planning** | `MRPL_MultiPeriod_LotSizing_MILP` | MILP | `8 × 12` (4 int) | `9720.0` | `0.00e+00` | **PASS** |
| **3. Refinery Scheduling** | `MRPL_Crude_Scheduling_MILP` | MILP | `35 × 24` (12 bin) | `5290.0` | `0.00e+00` | **PASS** |
| **4. Logistics & Supply Chain** | `MRPL_Logistics_SupplyChain_MILP` | MILP | `7 × 15` (3 bin) | `3190.0` | `0.00e+00` | **PASS** |
| **5. Power Dispatch** | `MRPL_Power_Economic_Dispatch_QP` | Convex QP | `1 × 3` (`3 Q-nnz`) | `4583.77777778` | `7.85e-16` | **PASS** |
| **6. Unit Commitment + Dispatch** | `MRPL_UnitCommitment_Dispatch_MIQP` | Convex MIQP | `3 × 4` (2 bin, `2 Q-nnz`) | `2238.63636364` | `7.20e-14` | **PASS** |
