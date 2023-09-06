#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP inclusion_probs(SEXP a, SEXP n);
extern SEXP up_brewer(SEXP pik, SEXP eps);

static const R_CallMethodDef CallEntries[] = {
    {"C_inclusion_probs", (DL_FUNC) &inclusion_probs, 2},
    {"C_up_brewer", (DL_FUNC) &up_brewer, 2},
    {NULL, NULL, 0}
};

void R_init_sondage(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
