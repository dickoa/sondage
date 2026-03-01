/*
 * Conditional Poisson Sampling (CPS)
 *
 * References:
 * - Chen, S.X., Dempster, A.P., Liu, J.S. (1994)
 * - Tillé, Y. (2006). Sampling Algorithms. Springer.
 */

#include "cps_core.h"

/* Convert complement indices (size m) to selected indices (size N - m). */
static int cps_complement_to_selected(const int *comp_idx, int m, int N,
                                      int *sel_idx, int n_expected) {
    int cpos = 0;
    int spos = 0;

    for (int i = 0; i < N; i++) {
        if (cpos < m) {
            int ci = comp_idx[cpos];
            if (ci < 0 || ci >= N) return -1;
            if (ci == i) {
                cpos++;
                continue;
            }
            if (ci < i) return -2; /* unsorted/duplicate complement index */
        }
        sel_idx[spos++] = i;
    }

    if (cpos != m) return -3;
    if (spos != n_expected) return -4;
    return 0;
}

SEXP C_cps_single(SEXP pik_sexp, SEXP eps_sexp) {
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
    const int use_complement = (n > N - n);
    const int n_work = use_complement ? (N - n) : n;

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

    /* If complement size is 0, all uncertain units are selected. */
    if (n_work == 0) {
        for (int i = 0; i < N; i++) {
            out[out_idx++] = idx_map[i] + 1;
        }
        UNPROTECT(1);
        return result;
    }

    double *pik_work = (double *) R_alloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        pik_work[i] = use_complement ? (1.0 - pik_valid[i]) : pik_valid[i];
    }

    double *w = (double *) R_alloc(N, sizeof(double));
    double *expa = (double *) R_alloc((size_t)N * n_work, sizeof(double));

    int cal_iters = cps_calibrate(pik_work, N, n_work, w, expa, 1e-9, 100);
    if (cal_iters >= 100) {
        Rf_warning("CPS calibration did not converge after %d iterations", 100);
    }

    int *sample_idx = (int *) R_alloc(n_work, sizeof(int));
    int *selected_idx = sample_idx;
    if (use_complement) {
        selected_idx = (int *) R_alloc(n, sizeof(int));
    }

    int selected = 0;
    GetRNGstate();
    selected = cps_sample(w, expa, N, n_work, sample_idx);
    PutRNGstate();

    if (selected != n_work) {
        UNPROTECT(1);
        error("CPS sampling produced %d units, expected %d", selected, n_work);
    }

    if (use_complement) {
        int rc = cps_complement_to_selected(sample_idx, n_work, N,
                                            selected_idx, n);
        if (rc != 0) {
            UNPROTECT(1);
            error("CPS complement conversion failed (code %d)", rc);
        }
    }

    for (int i = 0; i < n; i++) {
        if (selected_idx[i] < 0 || selected_idx[i] >= N) {
            UNPROTECT(1);
            error("CPS internal index out of range: %d (N=%d)",
                  selected_idx[i], N);
        }
        out[out_idx++] = idx_map[selected_idx[i]] + 1;
    }

    UNPROTECT(1);
    return result;
}

SEXP C_cps_design(SEXP pik_sexp, SEXP eps_sexp) {
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
    const int use_complement = (n > N - n);
    const int n_work = use_complement ? (N - n) : n;

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

    if (N > 0 && n_work > 0) {
        double *pik_work = (double *) R_alloc(N, sizeof(double));
        for (int i = 0; i < N; i++) {
            pik_work[i] = use_complement ? (1.0 - pik_valid[i]) : pik_valid[i];
        }

        w_arr = (double *) R_alloc(N, sizeof(double));
        expa_arr = (double *) R_alloc((size_t)N * n_work, sizeof(double));

        int cal_iters = cps_calibrate(pik_work, N, n_work, w_arr, expa_arr,
                                      1e-9, 100);
        if (cal_iters >= 100) {
            Rf_warning("CPS calibration did not converge after %d iterations",
                       100);
        }
    }

    SEXP result = PROTECT(allocVector(VECSXP, 10));
    SEXP names = PROTECT(allocVector(STRSXP, 10));

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
    if (N > 0 && n_work > 0) {
        memcpy(REAL(w_sexp), w_arr, N * sizeof(double));
    }
    SET_VECTOR_ELT(result, 6, w_sexp);
    SET_STRING_ELT(names, 6, mkChar("w"));

    SEXP expa_sexp = PROTECT(allocMatrix(REALSXP, N,
                                         (n_work > 0) ? n_work : 1));
    if (N > 0 && n_work > 0) {
        memcpy(REAL(expa_sexp), expa_arr,
               (size_t)N * n_work * sizeof(double));
    }
    SET_VECTOR_ELT(result, 7, expa_sexp);
    SET_STRING_ELT(names, 7, mkChar("expa"));

    SET_VECTOR_ELT(result, 8, ScalarInteger(use_complement));
    SET_STRING_ELT(names, 8, mkChar("use_complement"));

    SET_VECTOR_ELT(result, 9, ScalarInteger(n_work));
    SET_STRING_ELT(names, 9, mkChar("n_work"));

    setAttrib(result, R_NamesSymbol, names);
    setAttrib(result, R_ClassSymbol, mkString("cps_design"));

    UNPROTECT(6);
    return result;
}

SEXP C_cps_draw_batch(SEXP design, SEXP n_samples_sexp) {
    int N = INTEGER(VECTOR_ELT(design, 0))[0];
    int n = INTEGER(VECTOR_ELT(design, 1))[0];
    int N_certain = INTEGER(VECTOR_ELT(design, 4))[0];
    int use_complement = INTEGER(VECTOR_ELT(design, 8))[0];
    int n_work = INTEGER(VECTOR_ELT(design, 9))[0];
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

    if (n_work == 0) {
        for (int s = 0; s < n_samples; s++) {
            int *col = result_ptr + s * n_total + N_certain;
            for (int i = 0; i < N; i++) {
                col[i] = idx_map[i];
            }
        }
        UNPROTECT(1);
        return result;
    }

    int *sample_idx = (int *) R_alloc(n_work, sizeof(int));
    int *selected_idx = sample_idx;
    if (use_complement) {
        selected_idx = (int *) R_alloc(n, sizeof(int));
    }

    GetRNGstate();

    for (int s = 0; s < n_samples; s++) {
        if (s % 100 == 0) R_CheckUserInterrupt();
        int selected = cps_sample(w, expa, N, n_work, sample_idx);
        if (selected != n_work) {
            PutRNGstate();
            UNPROTECT(1);
            error(
                "CPS sampling produced %d units, expected %d (batch draw %d)",
                selected, n_work, s + 1
            );
        }

        if (use_complement) {
            int rc = cps_complement_to_selected(sample_idx, n_work, N,
                                                selected_idx, n);
            if (rc != 0) {
                PutRNGstate();
                UNPROTECT(1);
                error(
                    "CPS complement conversion failed (code %d, batch draw %d)",
                    rc, s + 1
                );
            }
        }

        int *col = result_ptr + s * n_total + N_certain;
        for (int i = 0; i < n; i++) {
            if (selected_idx[i] < 0 || selected_idx[i] >= N) {
                PutRNGstate();
                UNPROTECT(1);
                error(
                    "CPS internal index out of range: %d (N=%d, batch draw %d)",
                    selected_idx[i], N, s + 1
                );
            }
            col[i] = idx_map[selected_idx[i]];
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}
