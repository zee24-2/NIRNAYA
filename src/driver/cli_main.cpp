#include <iostream>
#include <string>
#include <vector>
#include "../../include/sov/solver.hpp"

static void print_usage() {
    std::cout << "Usage: sovsolve <model.mps> [options]\n"
              << "Options:\n"
              << "  --sol <path>              Write solution file (.sol)\n"
              << "  --stats-json <path>       Write solve statistics (.json)\n"
              << "  --time-limit <sec>        Wall-clock time limit in seconds (default 600)\n"
              << "  --method <m>              auto | dual | primal | ipm\n"
              << "  --presolve <0|1|2>        Presolve level (default 2)\n"
              << "  --mip-gap <val>           Relative MIP optimality gap (default 1e-4)\n"
              << "  --mip-abs-gap <val>       Absolute MIP optimality gap (default 1e-6)\n"
              << "  --node-limit <N>          Maximum B&B nodes\n"
              << "  --seed <N>                Random seed (default 1)\n"
              << "  --log <0|1|2|3>           Logging level (default 1)\n"
              << "  --no-scale                Disable matrix equilibration scaling\n"
              << "  --deterministic           Enable deterministic execution\n"
              << "  --tolerance-preset <p>    strict | default | relaxed\n"
              << "  --set <name>=<value>      Set arbitrary solver option\n"
              << "  --version                 Display version info\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 4;
    }

    sov::Options opt;
    std::string mps_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage();
            return 0;
        } else if (arg == "--version") {
            std::cout << "Sovereign Optimization Solver (sovsolve) v1.0.0 (C++17)\n";
            return 0;
        } else if (arg == "--sol" && i + 1 < argc) {
            opt.sol_path = argv[++i];
        } else if (arg == "--stats-json" && i + 1 < argc) {
            opt.stats_json = argv[++i];
        } else if (arg == "--time-limit" && i + 1 < argc) {
            opt.time_limit = std::stod(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            opt.threads = std::stoi(argv[++i]);
        } else if (arg == "--method" && i + 1 < argc) {
            opt.method = argv[++i];
        } else if (arg == "--presolve" && i + 1 < argc) {
            opt.presolve = std::stoi(argv[++i]);
        } else if (arg == "--mip-gap" && i + 1 < argc) {
            opt.tol.mip_gap = std::stod(argv[++i]);
        } else if (arg == "--mip-abs-gap" && i + 1 < argc) {
            opt.tol.mip_abs_gap = std::stod(argv[++i]);
        } else if (arg == "--node-limit" && i + 1 < argc) {
            opt.node_limit = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            opt.seed = std::stoull(argv[++i]);
        } else if (arg == "--log" && i + 1 < argc) {
            opt.log_level = std::stoi(argv[++i]);
        } else if (arg == "--no-scale") {
            opt.scale = false;
        } else if (arg == "--deterministic") {
            opt.deterministic = true;
        } else if (arg == "--tolerance-preset" && i + 1 < argc) {
            opt.apply_preset(argv[++i]);
        } else if (arg == "--set" && i + 1 < argc) {
            std::string kv = argv[++i];
            auto eq = kv.find('=');
            if (eq != std::string::npos) {
                opt.set_by_name(kv.substr(0, eq), kv.substr(eq + 1));
            }
        } else if (arg[0] != '-') {
            mps_path = arg;
        }
    }

    if (mps_path.empty()) {
        std::cerr << "Error: No input MPS model specified.\n";
        return 4;
    }

    sov::Model model;
    std::string err_msg;
    sov::Status rstat = sov::MpsReader::read(mps_path, model, err_msg);
    if (rstat != sov::Status::Optimal) {
        std::cerr << "Input Error reading " << mps_path << ": " << err_msg << "\n";
        return sov::status_to_exit_code(rstat);
    }

    sov::Solver solver(opt);
    sov::SolveResult res = solver.solve(model);
    return sov::status_to_exit_code(res.status);
}
