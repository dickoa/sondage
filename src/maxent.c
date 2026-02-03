/*
 * Maximum Entropy sampling / Conditional Poisson sampling
 *
 * References:
 * - Chen, S.X., Dempster, A.P., Liu, J.S. (1994)
 * - Tillé, Y. (2006). Sampling Algorithms. Springer.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>

static void compute_expa(const double *w, int N, int n, double *expa) {
    memset(expa, 0, (size_t)N * n * sizeof(double));

    double cumsum = 0.0;
    for (int i = N - 1; i >= 0; i--) {
        cumsum += w[i];
        expa[i * n + 0] = cumsum;
    }

    double logprod = 0.0;
    for (int i = N - 1; i >= N - n && i >= 0; i--) {
        logprod += log(w[i]);
        int z = N - i - 1;
        if (z < n) {
            expa[i * n + z] = exp(logprod);
        }
    }

    for (int i = N - 3; i >= 0; i--) {
        int max_z = N - i - 1;
        if (max_z > n - 1) max_z = n - 1;
        for (int z = 1; z <= max_z; z++) {
            expa[i * n + z] = w[i] * expa[(i + 1) * n + (z - 1)] +
                              expa[(i + 1) * n + z];
        }
    }
}

static inline double compute_q(const double *w, const double *expa,
                                int N, int n, int i, int z) {
    if (z == 0) {
        double denom = expa[i * n + 0];
        return (denom > 0) ? w[i] / denom : 0.0;
    }

    if (z >= N - i - 1) {
        return 1.0;
    }

    double denom = expa[i * n + z];
    if (denom <= 0) return 0.0;

    double numer = w[i] * expa[(i + 1) * n + (z - 1)];
    return numer / denom;
}

static void compute_pik_from_expa(const double *w, const double *expa,
                                   int N, int n, double *pik) {
    double *pro = (double *) R_alloc((size_t)N * n, sizeof(double));
    memset(pro, 0, (size_t)N * n * sizeof(double));

    pro[0 * n + (n - 1)] = 1.0;

    pik[0] = 0.0;
    for (int z = 0; z < n; z++) {
        double q_0z = compute_q(w, expa, N, n, 0, z);
        pik[0] += pro[0 * n + z] * q_0z;
    }

    for (int i = 1; i < N; i++) {
        for (int z = n - 1; z >= 0; z--) {
            double q_prev = compute_q(w, expa, N, n, i - 1, z);

            pro[i * n + z] += pro[(i - 1) * n + z] * (1.0 - q_prev);

            if (z + 1 < n) {
                double q_prev_zp1 = compute_q(w, expa, N, n, i - 1, z + 1);
                pro[i * n + z] += pro[(i - 1) * n + (z + 1)] * q_prev_zp1;
            }
        }

        pik[i] = 0.0;
        for (int z = 0; z < n; z++) {
            double q_iz = compute_q(w, expa, N, n, i, z);
            pik[i] += pro[i * n + z] * q_iz;
        }
    }
}

static int calibrate_piktilde(const double *pik_target, int N, int n,
                               double *piktilde, double *w, double *expa,
                               double eps, int max_iter) {
    double *pik_implied = (double *) R_alloc(N, sizeof(double));

    memcpy(piktilde, pik_target, N * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < N; i++) {
            double p = piktilde[i];
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            w[i] = p / (1.0 - p);
        }

        compute_expa(w, N, n, expa);

        compute_pik_from_expa(w, expa, N, n, pik_implied);

        double max_diff = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = pik_target[i] - pik_implied[i];
            piktilde[i] += diff;

            if (piktilde[i] < 1e-10) piktilde[i] = 1e-10;
            if (piktilde[i] > 1.0 - 1e-10) piktilde[i] = 1.0 - 1e-10;

            if (fabs(diff) > max_diff) max_diff = fabs(diff);
        }

        if (max_diff < eps) {
            for (int i = 0; i < N; i++) {
                double p = piktilde[i];
                if (p < 1e-10) p = 1e-10;
                if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
                w[i] = p / (1.0 - p);
            }
            compute_expa(w, N, n, expa);
            return iter + 1;
        }
    }

    return max_iter;
}

static void sample_sequential(const double *w, const double *expa,
                               int N, int n, int *sample_idx) {
    int remaining = n;
    int count = 0;

    for (int k = 0; k < N && remaining > 0; k++) {
        double q = compute_q(w, expa, N, n, k, remaining - 1);
        if (unif_rand() < q) {
            sample_idx[count++] = k;
            remaining--;
        }
    }
}

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

    double *piktilde = (double *) R_alloc(N, sizeof(double));
    double *w = (double *) R_alloc(N, sizeof(double));
    double *expa = (double *) R_alloc((size_t)N * n, sizeof(double));

    calibrate_piktilde(pik_valid, N, n, piktilde, w, expa, 1e-9, 100);

    int *sample_idx = (int *) R_alloc(n, sizeof(int));

    GetRNGstate();
    sample_sequential(w, expa, N, n, sample_idx);
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

    /* Allocate and fill arrays */
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

    double *piktilde = NULL;
    double *w_arr = NULL;
    double *expa_arr = NULL;

    if (N > 0 && n > 0) {
        piktilde = (double *) R_alloc(N, sizeof(double));
        w_arr = (double *) R_alloc(N, sizeof(double));
        expa_arr = (double *) R_alloc((size_t)N * n, sizeof(double));

        calibrate_piktilde(pik_valid, N, n, piktilde, w_arr, expa_arr, 1e-9, 100);
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

SEXP C_maxent_draw(SEXP design) {
    int N = INTEGER(VECTOR_ELT(design, 0))[0];
    int n = INTEGER(VECTOR_ELT(design, 1))[0];
    int N_certain = INTEGER(VECTOR_ELT(design, 4))[0];
    const int *idx_map = INTEGER(VECTOR_ELT(design, 3));
    const int *certain_idx = INTEGER(VECTOR_ELT(design, 5));
    const double *w = REAL(VECTOR_ELT(design, 6));
    const double *expa = REAL(VECTOR_ELT(design, 7));

    int n_total = n + N_certain;

    SEXP result = PROTECT(allocVector(INTSXP, n_total));
    int *out = INTEGER(result);
    int out_idx = 0;

    for (int k = 0; k < N_certain; k++) {
        out[out_idx++] = certain_idx[k];
    }

    if (N == 0 || n == 0) {
        UNPROTECT(1);
        return result;
    }

    int *sample_idx = (int *) R_alloc(n, sizeof(int));

    GetRNGstate();
    sample_sequential(w, expa, N, n, sample_idx);
    PutRNGstate();

    for (int i = 0; i < n; i++) {
        out[out_idx++] = idx_map[sample_idx[i]];
    }

    UNPROTECT(1);
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
        sample_sequential(w, expa, N, n, sample_idx);

        int *col = result_ptr + s * n_total + N_certain;
        for (int i = 0; i < n; i++) {
            col[i] = idx_map[sample_idx[i]];
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}
