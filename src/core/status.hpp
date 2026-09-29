#pragma once
#include <cstdint>
#include <string>

namespace sov {

enum class Status : int32_t {
    Optimal               = 0,
    Infeasible            = 1,
    Unbounded             = 2,
    InfeasibleOrUnbounded = 3,
    TimeLimit             = 4,
    WorkLimit             = 5,
    IterationLimit        = 6,
    NodeLimit             = 7,
    SolutionLimit         = 8,
    GapLimit              = 9,
    Interrupted           = 10,
    NumericalTrouble      = 11,
    InputError            = 12,
    NotSupported          = 13,
    InternalError         = 14,
    CutoffReached         = 15
};

inline const char* status_to_string(Status s) {
    switch (s) {
        case Status::Optimal:               return "optimal";
        case Status::Infeasible:            return "infeasible";
        case Status::Unbounded:             return "unbounded";
        case Status::InfeasibleOrUnbounded: return "infeasible_or_unbounded";
        case Status::TimeLimit:             return "time_limit";
        case Status::WorkLimit:             return "work_limit";
        case Status::IterationLimit:        return "iteration_limit";
        case Status::NodeLimit:             return "node_limit";
        case Status::SolutionLimit:         return "solution_limit";
        case Status::GapLimit:              return "gap_limit";
        case Status::Interrupted:           return "interrupted";
        case Status::NumericalTrouble:      return "numerical_trouble";
        case Status::InputError:            return "input_error";
        case Status::NotSupported:          return "not_supported";
        case Status::InternalError:         return "internal_error";
        case Status::CutoffReached:         return "cutoff_reached";
    }
    return "unknown";
}

inline int status_to_exit_code(Status s) {
    switch (s) {
        case Status::Optimal:
        case Status::GapLimit:
            return 0;
        case Status::TimeLimit:
        case Status::WorkLimit:
        case Status::IterationLimit:
        case Status::NodeLimit:
        case Status::SolutionLimit:
        case Status::Interrupted:
            return 1;
        case Status::Infeasible:
            return 2;
        case Status::Unbounded:
        case Status::InfeasibleOrUnbounded:
            return 3;
        case Status::NumericalTrouble:
            return 5;
        default:
            return 4;
    }
}

} // namespace sov
