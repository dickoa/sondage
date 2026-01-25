#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "maxent.h"

typedef struct {
    const double *pi_poisson;
    const int *idx;
    const int *certain_idx;
    int N;
    int n;
    int N_certain;
    int n_total;
} MaxentDesign;

static void extract_design(SEXP design, MaxentDesign *d) {
    d->pi_poisson = REAL(VECTOR_ELT(design, DESIGN_PI_POISSON));
    d->idx = INTEGER(VECTOR_ELT(design, DESIGN_IDX));
    d->N = INTEGER(VECTOR_ELT(design, DESIGN_N))[0];
    d->n = INTEGER(VECTOR_ELT(design, DESIGN_N_SAMPLE))[0];
    d->N_certain = INTEGER(VECTOR_ELT(design, DESIGN_N_CERTAIN))[0];
    d->certain_idx = (d->N_certain > 0) ? 
        INTEGER(VECTOR_ELT(design, DESIGN_CERTAIN_IDX)) : NULL;
    d->n_total = d->n + d->N_certain;
}

static int draw_poisson_sample(const MaxentDesign *d, int *out) {
    int count = 0;
    for (int k = 0; k < d->N; k++) {
        if (unif_rand() < d->pi_poisson[k]) {
            if (count < d->n) {
                out[count] = d->idx[k];
            }
            count++;
        }
    }
    return count;
}

SEXP C_maxent_sample(SEXP design, SEXP max_attempts_sexp) {
    MaxentDesign d;
    extract_design(design, &d);
    const int max_attempts = INTEGER(max_attempts_sexp)[0];

    SEXP result = PROTECT(allocVector(INTSXP, d.n_total));
    int *out = INTEGER(result);
    int out_idx = 0;

    for (int k = 0; k < d.N_certain; k++) {
        out[out_idx++] = d.certain_idx[k];
    }

    if (d.N == 0 || d.n == 0) {
        UNPROTECT(1);
        return result;
    }

    int *temp = (int *) R_alloc(d.n, sizeof(int));
    
    GetRNGstate();
    
    int attempts;
    for (attempts = 0; attempts < max_attempts; attempts++) {
        if (draw_poisson_sample(&d, temp) == d.n) {
            break;
        }
    }
    
    if (attempts == max_attempts) {
        warning("Maximum entropy sampling: no valid sample after %d attempts", max_attempts);
        draw_poisson_sample(&d, temp);
    }

    for (int i = 0; i < d.n; i++) {
        out[out_idx++] = temp[i];
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}

SEXP C_maxent_sample_batch(SEXP design, SEXP n_samples_sexp, SEXP max_attempts_sexp) {
    MaxentDesign d;
    extract_design(design, &d);
    const int n_samples = INTEGER(n_samples_sexp)[0];
    const int max_attempts = INTEGER(max_attempts_sexp)[0];

    SEXP result = PROTECT(allocMatrix(INTSXP, d.n_total, n_samples));
    int *result_ptr = INTEGER(result);

    for (int s = 0; s < n_samples; s++) {
        int *col = result_ptr + s * d.n_total;
        for (int k = 0; k < d.N_certain; k++) {
            col[k] = d.certain_idx[k];
        }
    }

    if (d.N == 0 || d.n == 0) {
        UNPROTECT(1);
        return result;
    }

    GetRNGstate();

    for (int s = 0; s < n_samples; s++) {
        int *col = result_ptr + s * d.n_total + d.N_certain;
        
        for (int attempts = 0; attempts < max_attempts; attempts++) {
            if (draw_poisson_sample(&d, col) == d.n) {
                break;
            }
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}
