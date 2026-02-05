/*
 * Maximum Entropy sampling / Conditional Poisson sampling
 *
 * References:
 * - Chen, S.X., Dempster, A.P., Liu, J.S. (1994)
 * - Tillé, Y. (2006). Sampling Algorithms. Springer.
 */

#include "cps_core.h"

SEXP C_maxent_single(SEXP pik_sexp, SEXP eps_sexp) {
    const double *pik_full = REAL(pik_sexp);
    const double epsilon = REAL(eps_sexp)[0];
    const int N_full = LENGTH(pik_sexp);

    int N = 0;
    int N_certain = 0;
    double n_sum = 0.0;

    for (int k = 0; k < N_full; k++) {
        double pk = pik_full[k];
        if (pk >= 1.0 - epsilon) {
            N_certain++;
        } else if (pk > epsilon) {
            N++;
            n_sum += pk;
        }
    }

    const int n = (int)(n_sum + 0.5);
    const int n_total = n + N_certain;

    SEXP result = PROTECT(allocVector(INTSXP, n_total));
    int *out = INTEGER(result);
    int out_idx = 0;

    double *pik_valid = NULL;
    int *idx_map = NULL;

    if (N > 0) {
        pik_valid = (double *) R_alloc(N, sizeof(double));
        idx_map = (int *) R_alloc(N, sizeof(int));
    }

    int j = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_full[k];
        if (pk >= 1.0 - epsilon) {
            out[out_idx++] = k + 1;
        } else if (pk > epsilon) {
            pik_valid[j] = pk;
            idx_map[j] = k;
            j++;
        }
    }

    if (N == 0 || n == 0) {
        UNPROTECT(1);
        return result;
    }

    double *w = (double *) R_alloc(N, sizeof(double));
    double *expa = (double *) R_alloc((size_t)N * n, sizeof(double));

    cps_calibrate(pik_valid, N, n, w, expa, 1e-9, 100);

    int *sample_idx = (int *) R_alloc(n, sizeof(int));

    GetRNGstate();
    cps_sample(w, expa, N, n, sample_idx);
    PutRNGstate();

    for (int i = 0; i < n; i++) {
        out[out_idx++] = idx_map[sample_idx[i]] + 1;
    }

    UNPROTECT(1);
    return result;
}

SEXP C_maxent_design(SEXP pik_sexp, SEXP eps_sexp) {
    const double *pik_full = REAL(pik_sexp);
    const double epsilon = REAL(eps_sexp)[0];
    const int N_full = LENGTH(pik_sexp);

    int N = 0;
    int N_certain = 0;
    double n_sum = 0.0;

    for (int k = 0; k < N_full; k++) {
        double pk = pik_full[k];
        if (pk >= 1.0 - epsilon) {
            N_certain++;
        } else if (pk > epsilon) {
            N++;
            n_sum += pk;
        }
    }

    const int n = (int)(n_sum + 0.5);

    double *pik_valid = NULL;
    int *idx_map = NULL;
    int *certain_idx = NULL;

    if (N > 0) {
        pik_valid = (double *) R_alloc(N, sizeof(double));
        idx_map = (int *) R_alloc(N, sizeof(int));
    }
    if (N_certain > 0) {
        certain_idx = (int *) R_alloc(N_certain, sizeof(int));
    }

    int j = 0, c = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_full[k];
        if (pk >= 1.0 - epsilon) {
            certain_idx[c++] = k;
        } else if (pk > epsilon) {
            pik_valid[j] = pk;
            idx_map[j] = k;
            j++;
        }
    }

    double *w_arr = NULL;
    double *expa_arr = NULL;

    if (N > 0 && n > 0) {
        w_arr = (double *) R_alloc(N, sizeof(double));
        expa_arr = (double *) R_alloc((size_t)N * n, sizeof(double));

        cps_calibrate(pik_valid, N, n, w_arr, expa_arr, 1e-9, 100);
    }

    SEXP result = PROTECT(allocVector(VECSXP, 8));
    SEXP names = PROTECT(allocVector(STRSXP, 8));

    SET_VECTOR_ELT(result, 0, ScalarInteger(N));
    SET_STRING_ELT(names, 0, mkChar("N"));

    SET_VECTOR_ELT(result, 1, ScalarInteger(n));
    SET_STRING_ELT(names, 1, mkChar("n"));

    SET_VECTOR_ELT(result, 2, ScalarInteger(N_full));
    SET_STRING_ELT(names, 2, mkChar("N_full"));

    SEXP idx_sexp = PROTECT(allocVector(INTSXP, N));
    for (int i = 0; i < N; i++) {
        INTEGER(idx_sexp)[i] = idx_map[i] + 1;
    }
    SET_VECTOR_ELT(result, 3, idx_sexp);
    SET_STRING_ELT(names, 3, mkChar("idx"));

    SET_VECTOR_ELT(result, 4, ScalarInteger(N_certain));
    SET_STRING_ELT(names, 4, mkChar("N_certain"));

    SEXP certain_sexp = PROTECT(allocVector(INTSXP, N_certain));
    for (int i = 0; i < N_certain; i++) {
        INTEGER(certain_sexp)[i] = certain_idx[i] + 1;
    }
    SET_VECTOR_ELT(result, 5, certain_sexp);
    SET_STRING_ELT(names, 5, mkChar("certain_idx"));

    SEXP w_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) {
        memcpy(REAL(w_sexp), w_arr, N * sizeof(double));
    }
    SET_VECTOR_ELT(result, 6, w_sexp);
    SET_STRING_ELT(names, 6, mkChar("w"));

    SEXP expa_sexp = PROTECT(allocMatrix(REALSXP, N, (n > 0) ? n : 1));
    if (N > 0 && n > 0) {
        memcpy(REAL(expa_sexp), expa_arr, (size_t)N * n * sizeof(double));
    }
    SET_VECTOR_ELT(result, 7, expa_sexp);
    SET_STRING_ELT(names, 7, mkChar("expa"));

    setAttrib(result, R_NamesSymbol, names);
    setAttrib(result, R_ClassSymbol, mkString("maxent_design_v2"));

    UNPROTECT(6);
    return result;
}

SEXP C_maxent_draw_batch(SEXP design, SEXP n_samples_sexp) {
    int N = INTEGER(VECTOR_ELT(design, 0))[0];
    int n = INTEGER(VECTOR_ELT(design, 1))[0];
    int N_certain = INTEGER(VECTOR_ELT(design, 4))[0];
    const int *idx_map = INTEGER(VECTOR_ELT(design, 3));
    const int *certain_idx = INTEGER(VECTOR_ELT(design, 5));
    const double *w = REAL(VECTOR_ELT(design, 6));
    const double *expa = REAL(VECTOR_ELT(design, 7));
    const int n_samples = INTEGER(n_samples_sexp)[0];

    int n_total = n + N_certain;

    SEXP result = PROTECT(allocMatrix(INTSXP, n_total, n_samples));
    int *result_ptr = INTEGER(result);

    for (int s = 0; s < n_samples; s++) {
        int *col = result_ptr + s * n_total;
        for (int k = 0; k < N_certain; k++) {
            col[k] = certain_idx[k];
        }
    }

    if (N == 0 || n == 0) {
        UNPROTECT(1);
        return result;
    }

    int *sample_idx = (int *) R_alloc(n, sizeof(int));

    GetRNGstate();

    for (int s = 0; s < n_samples; s++) {
        cps_sample(w, expa, N, n, sample_idx);

        int *col = result_ptr + s * n_total + N_certain;
        for (int i = 0; i < n; i++) {
            col[i] = idx_map[sample_idx[i]];
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}
