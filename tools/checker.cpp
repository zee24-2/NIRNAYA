#include <iostream>
#include <cstdio>
#include "../src/io/mps_reader.hpp"
#include "../src/io/sol_io.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: checker <model.mps> <solution.sol> [tol=1e-6]\n";
        return 1;
    }
    std::string mps_path = argv[1];
    std::string sol_path = argv[2];
    double tol = (argc >= 4) ? std::stod(argv[3]) : 1e-6;

    sov::Model model;
    std::string err;
    if (sov::MpsReader::read(mps_path, model, err) != sov::Status::Optimal) {
        std::cerr << "FAIL: Could not read model " << mps_path << ": " << err << "\n";
        return 1;
    }

    sov::Status claimed_status;
    double claimed_obj = 0.0;
    std::vector<double> x;
    if (!sov::SolIo::read_sol(sol_path, model, claimed_status, claimed_obj, x)) {
        std::cerr << "FAIL: Could not read solution " << sol_path << "\n";
        return 1;
    }

    auto rep = sov::verify_solution(model, x, claimed_obj, tol);
    std::printf("================================================================\n");
    std::printf("  SOVEREIGN INDEPENDENT SOLUTION CHECKER (dd_real ~32 digits)\n");
    std::printf("  Model:             %s (%d rows, %d cols)\n", model.name.c_str(), model.nrow, model.ncol);
    std::printf("  Claimed Objective: %.17g\n", claimed_obj);
    std::printf("  True Objective:    %.17g (rel err: %.3e)\n", rep.computed_obj, rep.obj_rel_err);
    std::printf("  Max Bound Viol:    %.3e\n", rep.max_bound_viol);
    std::printf("  Max Row Viol:      %.3e\n", rep.max_row_viol);
    std::printf("  Max Int Viol:      %.3e\n", rep.max_int_viol);
    std::printf("  Verdict:           %s\n", rep.passed ? "PASS" : ("FAIL (" + rep.worst_item + ")").c_str());
    std::printf("================================================================\n");
    return rep.passed ? 0 : 1;
}
