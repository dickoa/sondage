#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP C_up_brewer(SEXP pik, SEXP eps);
extern SEXP C_inclusion_prob(SEXP a, SEXP n);
extern SEXP C_maxent_design_create(SEXP pik, SEXP eps);
extern SEXP C_maxent_sample(SEXP design, SEXP max_attempts);
extern SEXP C_maxent_sample_batch(SEXP design, SEXP n_samples, SEXP max_attempts);
extern SEXP C_cube(SEXP prob, SEXP X, SEXP eps);
extern SEXP C_cube_batch(SEXP prob, SEXP X, SEXP eps, SEXP n_samples);
extern SEXP C_cube_stratified(SEXP prob, SEXP X, SEXP strata, SEXP eps);

static const R_CallMethodDef CallEntries[] = {
    {"C_up_brewer",            (DL_FUNC) &C_up_brewer,            2},
    {"C_inclusion_prob",       (DL_FUNC) &C_inclusion_prob,       2},
    {"C_maxent_design_create", (DL_FUNC) &C_maxent_design_create, 2},
    {"C_maxent_sample",        (DL_FUNC) &C_maxent_sample,        2},
    {"C_maxent_sample_batch",  (DL_FUNC) &C_maxent_sample_batch,  3},
    {"C_cube",                 (DL_FUNC) &C_cube,                 3},
    {"C_cube_batch",           (DL_FUNC) &C_cube_batch,           4},
    {"C_cube_stratified",      (DL_FUNC) &C_cube_stratified,      4},
    {NULL, NULL, 0}
};

void R_init_sondage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
