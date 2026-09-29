#include "../../include/sov/sov.h"
#include "../../include/sov/solver.hpp"
#include <cstring>

struct sov_env {
    sov::Options default_opt;
};

struct sov_model {
    sov_env* env = nullptr;
    sov::ModelBuilder builder;
    sov::Model finalized;
    bool is_finalized = false;
    sov::Solver solver;
    sov::SolveResult result;
    std::string last_error;

    explicit sov_model(sov_env* e, const char* name)
        : env(e), builder(name ? name : "model"), solver(e ? e->default_opt : sov::Options()) {}
};

extern "C" {

sov_env* sov_env_create(void) {
    try {
        return new sov_env();
    } catch (...) {
        return nullptr;
    }
}

void sov_env_free(sov_env* env) {
    delete env;
}

sov_model* sov_model_create(sov_env* env, const char* name) {
    try {
        return new sov_model(env, name);
    } catch (...) {
        return nullptr;
    }
}

void sov_model_free(sov_model* model) {
    delete model;
}

int sov_add_col(sov_model* model, double cost, double lo, double up, int type, const char* name) {
    if (!model) return -1;
    try {
        model->is_finalized = false;
        sov::VarType vt = sov::VarType::Continuous;
        if (type == 1) vt = sov::VarType::Integer;
        else if (type == 2) vt = sov::VarType::Binary;
        return static_cast<int>(model->builder.add_col(cost, lo, up, vt, name ? name : ""));
    } catch (const std::exception& e) {
        model->last_error = e.what();
        return -1;
    }
}

int sov_add_row(sov_model* model, double lo, double up, int nnz, const int* idx, const double* val, const char* name) {
    if (!model) return -1;
    try {
        model->is_finalized = false;
        std::vector<sov::Idx> indices(idx, idx + nnz);
        std::vector<double> values(val, val + nnz);
        return static_cast<int>(model->builder.add_row(lo, up, indices, values, name ? name : ""));
    } catch (const std::exception& e) {
        model->last_error = e.what();
        return -1;
    }
}

int sov_set_sense(sov_model* model, int sense) {
    if (!model) return -1;
    model->builder.set_sense(sense < 0 ? sov::Sense::Maximize : sov::Sense::Minimize);
    return 0;
}

int sov_set_offset(sov_model* model, double offset) {
    if (!model) return -1;
    model->builder.set_offset(offset);
    return 0;
}

int sov_set_quad(sov_model* model, int nnz, const int* i, const int* j, const double* v) {
    if (!model) return -1;
    try {
        model->is_finalized = false;
        for (int k = 0; k < nnz; ++k) {
            model->builder.add_quad(i[k], j[k], v[k]);
        }
        return 0;
    } catch (const std::exception& e) {
        model->last_error = e.what();
        return -1;
    }
}

int sov_read_mps(sov_model* model, const char* path) {
    if (!model || !path) return -1;
    std::string err;
    auto st = sov::MpsReader::read(path, model->finalized, err);
    if (st != sov::Status::Optimal) {
        model->last_error = err;
        return -1;
    }
    model->is_finalized = true;
    return 0;
}

int sov_write_mps(sov_model* model, const char* path) {
    if (!model || !path) return -1;
    if (!model->is_finalized) {
        model->finalized = model->builder.finalize();
        model->is_finalized = true;
    }
    return sov::MpsReader::write(path, model->finalized) ? 0 : -1;
}

int sov_set_option_string(sov_model* model, const char* name, const char* value) {
    if (!model || !name || !value) return -1;
    return model->solver.options.set_by_name(name, value) ? 0 : -1;
}

int sov_optimize(sov_model* model) {
    if (!model) return -1;
    try {
        if (!model->is_finalized) {
            model->finalized = model->builder.finalize();
            model->is_finalized = true;
        }
        model->result = model->solver.solve(model->finalized);
        model->last_error = model->result.message;
        return 0;
    } catch (const std::exception& e) {
        model->last_error = e.what();
        return -1;
    }
}

int sov_interrupt(sov_model* model) {
    if (!model) return -1;
    model->solver.interrupt();
    return 0;
}

int sov_get_status(sov_model* model) {
    if (!model) return static_cast<int>(sov::Status::InternalError);
    return static_cast<int>(model->result.status);
}

double sov_get_obj(sov_model* model) {
    if (!model) return 0.0;
    return model->result.objective;
}

double sov_get_bound(sov_model* model) {
    if (!model) return 0.0;
    return model->result.best_bound;
}

int sov_get_x(sov_model* model, double* out_x) {
    if (!model || !out_x) return -1;
    for (size_t j = 0; j < model->result.x.size(); ++j) {
        out_x[j] = model->result.x[j];
    }
    return 0;
}

const char* sov_error_message(sov_model* model) {
    if (!model) return "Null model handle";
    return model->last_error.c_str();
}

} // extern "C"
