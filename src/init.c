#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP C_up_brewer(SEXP pik, SEXP eps);
extern SEXP C_inclusion_prob(SEXP a, SEXP n);
/* Maximum entropy (sequential algorithm) */
extern SEXP C_maxent_single(SEXP pik, SEXP eps);
extern SEXP C_maxent_design(SEXP pik, SEXP eps);
extern SEXP C_maxent_draw_batch(SEXP design, SEXP n_samples);
/* Chromy sequential sampling */
extern SEXP C_up_chromy(SEXP x, SEXP n);
extern SEXP C_up_chromy_joint(SEXP x, SEXP n, SEXP nsim);
/* Joint inclusion probabilities */
extern SEXP C_up_maxent_joint(SEXP pik, SEXP eps);
extern SEXP C_up_systematic_joint(SEXP pik, SEXP eps);
extern SEXP C_up_brewer_joint(SEXP pik);

static const R_CallMethodDef CallEntries[] = {
    {"C_up_brewer",               (DL_FUNC) &C_up_brewer,               2},
    {"C_inclusion_prob",          (DL_FUNC) &C_inclusion_prob,          2},
    /* Maximum entropy */
    {"C_maxent_single",           (DL_FUNC) &C_maxent_single,           2},
    {"C_maxent_design",           (DL_FUNC) &C_maxent_design,           2},
    {"C_maxent_draw_batch",       (DL_FUNC) &C_maxent_draw_batch,       2},
    /* Chromy sequential */
    {"C_up_chromy",               (DL_FUNC) &C_up_chromy,               2},
    {"C_up_chromy_joint",         (DL_FUNC) &C_up_chromy_joint,         3},
    /* Joint prob */
    {"C_up_maxent_joint",         (DL_FUNC) &C_up_maxent_joint,         2},
    {"C_up_systematic_joint",     (DL_FUNC) &C_up_systematic_joint,     2},
    {"C_up_brewer_joint",         (DL_FUNC) &C_up_brewer_joint,         1},
    {NULL, NULL, 0}
};

void R_init_sondage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
