#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include "../include/sov/solver.hpp"

namespace sov {

// 1. Crude Blending (LP) -- Slide 2 "Crude blending"
Model build_crude_blending_lp() {
    ModelBuilder mb("MRPL_Crude_Blending_LP");
    mb.set_sense(Sense::Maximize);

    Idx u0 = mb.add_col(-72.0, 0.0, 500.0, VarType::Continuous, "Crude_ArabLight");
    Idx u1 = mb.add_col(-78.0, 0.0, 500.0, VarType::Continuous, "Crude_BonnyLight");
    Idx u2 = mb.add_col(-68.0, 0.0, 500.0, VarType::Continuous, "Crude_IranHeavy");

    Idx s0 = mb.add_col(95.0,  0.0, 1000.0, VarType::Continuous, "Sell_Naphtha");
    Idx s1 = mb.add_col(105.0, 0.0, 1000.0, VarType::Continuous, "Sell_JetKero");
    Idx s2 = mb.add_col(102.0, 0.0, 1000.0, VarType::Continuous, "Sell_Gasoil");

    mb.add_row(-kInf, 1000.0, {u0, u1, u2}, {1.0, 1.0, 1.0}, "CDU_Capacity");
    mb.add_row(0.0, 0.0, {u0, u1, u2, s0}, {0.22, 0.28, 0.18, -1.0}, "Yield_Naphtha");
    mb.add_row(0.0, 0.0, {u0, u1, u2, s1}, {0.30, 0.34, 0.25, -1.0}, "Yield_JetKero");
    mb.add_row(0.0, 0.0, {u0, u1, u2, s2}, {0.42, 0.35, 0.48, -1.0}, "Yield_Gasoil");
    mb.add_row(-kInf, 0.0, {u0, u1, u2},
               {(0.60 - 0.50) * 0.42, (0.15 - 0.50) * 0.35, (0.90 - 0.50) * 0.48},
               "Gasoil_Sulfur_Spec");
    mb.add_row(300.0, kInf, {s2}, {1.0}, "Min_Gasoil_Demand");

    return mb.finalize();
}

// 2. Production Planning (MILP) -- Slide 2 "Production planning"
Model build_multiperiod_lotsizing_milp() {
    ModelBuilder mb("MRPL_MultiPeriod_LotSizing_MILP");
    mb.set_sense(Sense::Minimize);

    const int T = 4;
    const double demand[T] = {120.0, 180.0, 150.0, 210.0};
    const double cap = 250.0;
    std::vector<Idx> x(T), inv(T), z(T);

    for (int t = 0; t < T; ++t) {
        x[t]   = mb.add_col(12.0,  0.0, cap,   VarType::Continuous, "Prod_t" + std::to_string(t + 1));
        inv[t] = mb.add_col(2.5,   0.0, 300.0, VarType::Continuous, "Inv_t" + std::to_string(t + 1));
        z[t]   = mb.add_col(450.0, 0.0, 1.0,   VarType::Binary,     "Setup_t" + std::to_string(t + 1));
    }

    for (int t = 0; t < T; ++t) {
        if (t == 0) {
            mb.add_row(demand[t], demand[t], {x[t], inv[t]}, {1.0, -1.0}, "Bal_t1");
        } else {
            mb.add_row(demand[t], demand[t], {inv[t - 1], x[t], inv[t]}, {1.0, 1.0, -1.0}, "Bal_t" + std::to_string(t + 1));
        }
        mb.add_row(-kInf, 0.0, {x[t], z[t]}, {1.0, -cap}, "CapLink_t" + std::to_string(t + 1));
    }
    return mb.finalize();
}

// 3. Refinery Scheduling (MILP) -- Slide 2 "Refinery scheduling"
Model build_crude_scheduling_milp() {
    ModelBuilder mb("MRPL_Crude_Scheduling_MILP");
    mb.set_sense(Sense::Minimize);

    const int K = 3, T = 4;
    const double tank_cost[K] = {10.0, 14.0, 18.0};
    const double tank_init_inv[K] = {200.0, 180.0, 150.0};

    std::vector<std::vector<Idx>> b(K, std::vector<Idx>(T));
    std::vector<std::vector<Idx>> q(K, std::vector<Idx>(T));

    for (int k = 0; k < K; ++k) {
        for (int t = 0; t < T; ++t) {
            b[k][t] = mb.add_col(50.0 * (k + 1), 0.0, 1.0, VarType::Binary,
                                 "FeedBin_k" + std::to_string(k + 1) + "_t" + std::to_string(t + 1));
            q[k][t] = mb.add_col(tank_cost[k], 0.0, 120.0, VarType::Continuous,
                                 "Flow_k" + std::to_string(k + 1) + "_t" + std::to_string(t + 1));
        }
    }

    for (int t = 0; t < T; ++t) {
        mb.add_row(1.0, 1.0, {b[0][t], b[1][t], b[2][t]}, {1.0, 1.0, 1.0}, "OneTank_t" + std::to_string(t + 1));
        mb.add_row(95.0, 120.0, {q[0][t], q[1][t], q[2][t]}, {1.0, 1.0, 1.0}, "CDU_Demand_t" + std::to_string(t + 1));
        for (int k = 0; k < K; ++k) {
            mb.add_row(-kInf, 0.0, {q[k][t], b[k][t]}, {1.0, -120.0}, "MaxRate_k" + std::to_string(k) + "_t" + std::to_string(t));
            mb.add_row(0.0, kInf,  {q[k][t], b[k][t]}, {1.0, -80.0},  "MinRate_k" + std::to_string(k) + "_t" + std::to_string(t));
        }
    }

    for (int k = 0; k < K; ++k) {
        mb.add_row(-kInf, tank_init_inv[k], {q[k][0], q[k][1], q[k][2], q[k][3]}, {1.0, 1.0, 1.0, 1.0},
                   "TankInv_k" + std::to_string(k + 1));
    }
    return mb.finalize();
}

// 4. Logistics & Supply Chain Network (Fixed-Charge Facility Location MILP) -- Slide 2 "Logistics and supply chain"
Model build_logistics_supply_chain_milp() {
    ModelBuilder mb("MRPL_Logistics_SupplyChain_MILP");
    mb.set_sense(Sense::Minimize);

    // 3 Coastal/Inland Depots (Mangaluru, Hassan, Mysuru) -> 4 Demand Zones
    const int I = 3, J = 4;
    const double fixed_cost[I] = {1200.0, 900.0, 1050.0};
    const double cap[I] = {160.0, 140.0, 150.0};
    const double demand[J] = {70.0, 85.0, 65.0, 60.0}; // total demand 280 <= any two depots
    const double ship_cost[I][J] = {
        {4.0, 7.0, 8.5, 6.0},
        {6.5, 3.5, 5.0, 7.5},
        {8.0, 5.5, 3.0, 4.5}
    };

    std::vector<Idx> y(I);
    std::vector<std::vector<Idx>> f(I, std::vector<Idx>(J));
    for (int i = 0; i < I; ++i) {
        y[i] = mb.add_col(fixed_cost[i], 0.0, 1.0, VarType::Binary, "OpenDepot_" + std::to_string(i + 1));
        for (int j = 0; j < J; ++j) {
            f[i][j] = mb.add_col(ship_cost[i][j], 0.0, demand[j], VarType::Continuous,
                                 "Ship_d" + std::to_string(i + 1) + "_z" + std::to_string(j + 1));
        }
    }

    // Demand satisfaction per zone j
    for (int j = 0; j < J; ++j) {
        mb.add_row(demand[j], demand[j], {f[0][j], f[1][j], f[2][j]}, {1.0, 1.0, 1.0}, "ZoneDemand_" + std::to_string(j + 1));
    }
    // Depot capacity linked to open binary y_i: sum_j f_{i,j} - cap_i * y_i <= 0
    for (int i = 0; i < I; ++i) {
        mb.add_row(-kInf, 0.0, {f[i][0], f[i][1], f[i][2], f[i][3], y[i]},
                   {1.0, 1.0, 1.0, 1.0, -cap[i]}, "DepotCap_" + std::to_string(i + 1));
    }
    return mb.finalize();
}

// 5. Power Dispatch (Convex QP) -- Slide 2 "Power dispatch"
Model build_economic_dispatch_qp() {
    ModelBuilder mb("MRPL_Power_Economic_Dispatch_QP");
    mb.set_sense(Sense::Minimize);

    Idx p0 = mb.add_col(18.0, 20.0, 90.0, VarType::Continuous, "Gen1_GT");
    Idx p1 = mb.add_col(22.0, 15.0, 80.0, VarType::Continuous, "Gen2_STG");
    Idx p2 = mb.add_col(26.0, 10.0, 70.0, VarType::Continuous, "Gen3_Aux");

    mb.add_quad(p0, p0, 0.12);
    mb.add_quad(p1, p1, 0.16);
    mb.add_quad(p2, p2, 0.20);

    mb.add_row(180.0, 180.0, {p0, p1, p2}, {1.0, 1.0, 1.0}, "Power_Balance_180MW");
    return mb.finalize();
}

// 6. Unit Commitment + Quadratic Dispatch (Convex MIQP) -- Slide 2 & Slide 3 "Extension-ready: MIQP"
Model build_unit_commitment_miqp() {
    ModelBuilder mb("MRPL_UnitCommitment_Dispatch_MIQP");
    mb.set_sense(Sense::Minimize);

    // 2 Generators with binary commitment u_0, u_1 and quadratic dispatch p_0, p_1 meeting 100 MW
    Idx u0 = mb.add_col(200.0, 0.0, 1.0, VarType::Binary, "Commit_GT1");
    Idx u1 = mb.add_col(150.0, 0.0, 1.0, VarType::Binary, "Commit_GT2");
    Idx p0 = mb.add_col(15.0,  0.0, 80.0, VarType::Continuous, "Gen_GT1");
    Idx p1 = mb.add_col(18.0,  0.0, 80.0, VarType::Continuous, "Gen_GT2");

    mb.add_quad(p0, p0, 0.10);
    mb.add_quad(p1, p1, 0.12);

    mb.add_row(100.0, 100.0, {p0, p1}, {1.0, 1.0}, "Load_100MW");
    mb.add_row(-kInf, 0.0,   {p0, u0}, {1.0, -80.0}, "Cap_GT1");
    mb.add_row(-kInf, 0.0,   {p1, u1}, {1.0, -80.0}, "Cap_GT2");
    return mb.finalize();
}

} // namespace sov

int main() {
    std::vector<sov::Model> models = {
        sov::build_crude_blending_lp(),
        sov::build_multiperiod_lotsizing_milp(),
        sov::build_crude_scheduling_milp(),
        sov::build_logistics_supply_chain_milp(),
        sov::build_economic_dispatch_qp(),
        sov::build_unit_commitment_miqp()
    };

    std::printf("========================================================================\n");
    std::printf("  NIRNAYA (SIH26119) -- ALL 5 INDUSTRIAL DOMAINS + MIQP SHOWCASE\n");
    std::printf("========================================================================\n");

    int failures = 0;
    for (const auto& m : models) {
        std::string mps_file = "data/refinery/" + m.name + ".mps";
        std::string sol_file = "data/refinery/" + m.name + ".sol";
        sov::MpsReader::write(mps_file, m);

        sov::Options opt;
        opt.log_level = 1;
        opt.sol_path = sol_file;
        sov::Solver solver(opt);
        auto res = solver.solve(m);
        if (res.status != sov::Status::Optimal || !res.verification.passed) {
            ++failures;
        }
    }
    return failures;
}
