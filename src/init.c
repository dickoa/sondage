#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP C_up_brewer(SEXP pik, SEXP eps);
extern SEXP C_inclusion_prob(SEXP a, SEXP n);
extern SEXP C_cps_single(SEXP pik, SEXP eps);
extern SEXP C_cps_design(SEXP pik, SEXP eps);
extern SEXP C_cps_draw_batch(SEXP design, SEXP n_samples);
extern SEXP C_up_chromy(SEXP x, SEXP n);
extern SEXP C_chromy_joint_exp(SEXP x, SEXP n, SEXP nsim);
extern SEXP C_cps_jip(SEXP pik, SEXP eps);
extern SEXP C_up_systematic_jip(SEXP pik, SEXP eps);
extern SEXP C_high_entropy_jip(SEXP pik, SEXP eps);
extern SEXP C_cube(SEXP prob, SEXP X, SEXP eps);
extern SEXP C_cube_stratified(SEXP prob, SEXP X, SEXP strata, SEXP eps);

static const R_CallMethodDef CallEntries[] = {
    {"C_up_brewer",               (DL_FUNC) &C_up_brewer,               2},
    {"C_inclusion_prob",          (DL_FUNC) &C_inclusion_prob,          2},
    {"C_cps_single",              (DL_FUNC) &C_cps_single,              2},
    {"C_cps_design",              (DL_FUNC) &C_cps_design,              2},
    {"C_cps_draw_batch",          (DL_FUNC) &C_cps_draw_batch,          2},
    {"C_up_chromy",               (DL_FUNC) &C_up_chromy,               2},
    {"C_chromy_joint_exp",        (DL_FUNC) &C_chromy_joint_exp,        3},
    {"C_cps_jip",                 (DL_FUNC) &C_cps_jip,                 2},
    {"C_up_systematic_jip",     (DL_FUNC) &C_up_systematic_jip,     2},
    {"C_high_entropy_jip",      (DL_FUNC) &C_high_entropy_jip,      2},
    {"C_cube",                  (DL_FUNC) &C_cube,                  3},
    {"C_cube_stratified",       (DL_FUNC) &C_cube_stratified,       4},
    {NULL, NULL, 0}
};

void R_init_sondage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
