#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>

/*
 * Maximum Entropy Sampling (Conditional Poisson Sampling)
 *
 * Sampling functions that use precomputed design objects.
 * Design objects are created by C_maxent_design_create in maxent_design.c
 */

/* ============================================================
 * Sample from precomputed design (rejective method)
 *
 * Uses rejection sampling with Poisson probabilities.
 * Expected attempts: O(sqrt(n))
 *
 * Handles certainty selections (pik >= 1-eps)
 * Returns integer vector of selected indices (1-indexed for R)
 * ============================================================ */
SEXP C_maxent_sample(SEXP design, SEXP max_attempts_sexp) {
    /* Extract design components */
    SEXP pi_poisson_sexp = VECTOR_ELT(design, 3);
    const double *pi_poisson = REAL(pi_poisson_sexp);
    const int N = INTEGER(VECTOR_ELT(design, 4))[0];
    const int n = INTEGER(VECTOR_ELT(design, 5))[0];
    const int *idx = INTEGER(VECTOR_ELT(design, 6));
    const int max_attempts = INTEGER(max_attempts_sexp)[0];

    /* Get certainty selection info */
    int N_certain = 0;
    const int *certain_idx = NULL;
    if (LENGTH(design) >= 11) {
        N_certain = INTEGER(VECTOR_ELT(design, 9))[0];
        if (N_certain > 0) {
            certain_idx = INTEGER(VECTOR_ELT(design, 10));
        }
    }

    const int total_selected = n + N_certain;

    /* Allocate output - indices */
    SEXP s = PROTECT(allocVector(INTSXP, total_selected));
    int *s_ptr = INTEGER(s);
    int s_idx = 0;

    /* Always include certainty selections first */
    for (int k = 0; k < N_certain; k++) {
        s_ptr[s_idx++] = certain_idx[k];  /* Already 1-indexed */
    }

    /* If no random units to sample, we're done */
    if (N == 0 || n == 0) {
        UNPROTECT(1);
        return s;
    }

    /* Temporary array for random unit selection */
    int *temp_selected = (int *) R_alloc(n, sizeof(int));

    /* Rejection sampling for random units */
    GetRNGstate();

    int attempts = 0;
    int success = 0;

    while (attempts < max_attempts && !success) {
        attempts++;

        int sample_size = 0;
        for (int k = 0; k < N; k++) {
            if (unif_rand() < pi_poisson[k]) {
                if (sample_size < n) {
                    temp_selected[sample_size] = idx[k];  /* idx is 1-indexed */
                }
                sample_size++;
            }
        }

        if (sample_size == n) {
            success = 1;
        }
    }

    /* If not successful, warn */
    if (!success) {
        /* Try one more time */
        int sample_size = 0;
        for (int k = 0; k < N; k++) {
            if (unif_rand() < pi_poisson[k]) {
                if (sample_size < n) {
                    temp_selected[sample_size] = idx[k];
                }
                sample_size++;
            }
        }
        warning("Maximum entropy sampling: failed after %d attempts", max_attempts);
    }

    /* Copy selected indices to output */
    for (int i = 0; i < n; i++) {
        s_ptr[s_idx++] = temp_selected[i];
    }

    PutRNGstate();

    UNPROTECT(1);
    return s;
}

/* ============================================================
 * Batch sampling - draw multiple samples efficiently
 *
 * Handles certainty selections (pik >= 1-eps)
 * Returns integer matrix of indices (n_total x n_samples)
 * where n_total = n + N_certain
 * ============================================================ */
SEXP C_maxent_sample_batch(SEXP design, SEXP n_samples_sexp, SEXP max_attempts_sexp) {
    const double *pi_poisson = REAL(VECTOR_ELT(design, 3));
    const int N = INTEGER(VECTOR_ELT(design, 4))[0];
    const int n = INTEGER(VECTOR_ELT(design, 5))[0];
    const int *idx = INTEGER(VECTOR_ELT(design, 6));
    const int n_samples = INTEGER(n_samples_sexp)[0];
    const int max_attempts = INTEGER(max_attempts_sexp)[0];

    /* Get certainty selection info */
    int N_certain = 0;
    const int *certain_idx = NULL;
    if (LENGTH(design) >= 11) {
        N_certain = INTEGER(VECTOR_ELT(design, 9))[0];
        if (N_certain > 0) {
            certain_idx = INTEGER(VECTOR_ELT(design, 10));
        }
    }

    const int n_total = n + N_certain;

    /* Allocate output matrix (n_total x n_samples) */
    SEXP result = PROTECT(allocMatrix(INTSXP, n_total, n_samples));
    int *result_ptr = INTEGER(result);

    /* Set certainty selections in all columns (first N_certain rows) */
    for (int s = 0; s < n_samples; s++) {
        int *col = result_ptr + s * n_total;
        for (int k = 0; k < N_certain; k++) {
            col[k] = certain_idx[k];  /* Already 1-indexed */
        }
    }

    /* If no random units to sample, we're done */
    if (N == 0 || n == 0) {
        UNPROTECT(1);
        return result;
    }

    GetRNGstate();

    for (int s = 0; s < n_samples; s++) {
        int *col = result_ptr + s * n_total;
        int col_idx = N_certain;  /* Start after certainty selections */

        int attempts = 0;
        int success = 0;

        while (attempts < max_attempts && !success) {
            attempts++;
            col_idx = N_certain;  /* Reset */

            int sample_size = 0;
            for (int k = 0; k < N; k++) {
                if (unif_rand() < pi_poisson[k]) {
                    if (sample_size < n) {
                        col[col_idx++] = idx[k];  /* idx is 1-indexed */
                    }
                    sample_size++;
                }
            }

            if (sample_size == n) {
                success = 1;
            }
        }
    }

    PutRNGstate();

    UNPROTECT(1);
    return result;
}
