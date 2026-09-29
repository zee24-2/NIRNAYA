#pragma once
#ifdef __cplusplus
extern "C" {
#endif

typedef struct sov_env sov_env;
typedef struct sov_model sov_model;

sov_env*   sov_env_create(void);
void       sov_env_free(sov_env* env);

sov_model* sov_model_create(sov_env* env, const char* name);
void       sov_model_free(sov_model* model);

int sov_add_col(sov_model* model, double cost, double lo, double up, int type, const char* name);
int sov_add_row(sov_model* model, double lo, double up, int nnz, const int* idx, const double* val, const char* name);
int sov_set_sense(sov_model* model, int sense); // 1 = Minimize, -1 = Maximize
int sov_set_offset(sov_model* model, double offset);
int sov_set_quad(sov_model* model, int nnz, const int* i, const int* j, const double* v);

int sov_read_mps(sov_model* model, const char* path);
int sov_write_mps(sov_model* model, const char* path);

int sov_set_option_string(sov_model* model, const char* name, const char* value);
int sov_optimize(sov_model* model);
int sov_interrupt(sov_model* model);

int    sov_get_status(sov_model* model);
double sov_get_obj(sov_model* model);
double sov_get_bound(sov_model* model);
int    sov_get_x(sov_model* model, double* out_x);
const char* sov_error_message(sov_model* model);

#ifdef __cplusplus
}
#endif
