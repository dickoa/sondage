/*
 * joint_probs.c - Joint Inclusion Probabilities
 *
 * Implements:
 * - C_up_maxent_jip: Exact CPS joint probabilities (Aires' formula)
 * - C_up_systematic_jip: Exact systematic sampling joint probabilities
 * - C_up_brewer_jip: Brewer approximation (equation 18, Brewer & Donadio 2003)
 *
 * References:
 * - Aires (1999). Algorithms to Find Exact Inclusion Probabilities for
 *   Conditional Poisson Sampling and Pareto pips Sampling Design
 * - Tillé (2006). Sampling Algorithms. Springer.
 * - Brewer & Donadio (2003). The High Entropy Variance of the HT Estimator.
 */

#include "cps_core.h"
#include <float.h>

static int compare_double(const void *a, const void *b) {
    double diff = *(const double *)a - *(const double *)b;
    return (diff > 0) - (diff < 0);
}

SEXP C_up_maxent_jip(SEXP pik_sexp, SEXP eps_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);
    const double eps = REAL(eps_sexp)[0];

    double n_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        n_sum += pik_full[k];
    }
    const int n = (int)(n_sum + 0.5);

    SEXP result = PROTECT(allocMatrix(REALSXP, N_full, N_full));
    double *pikl = REAL(result);
    memset(pikl, 0, (size_t)N_full * N_full * sizeof(double));

    for (int k = 0; k < N_full; k++) {
        pikl[k * N_full + k] = pik_full[k];
    }

    if (n < 2) {
        UNPROTECT(1);
        return result;
    }

    int N_valid = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            N_valid++;
        }
    }

    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] >= 1.0 - eps) {
            for (int l = 0; l < N_full; l++) {
                pikl[k * N_full + l] = pik_full[l];
                pikl[l * N_full + k] = pik_full[l];
            }
            pikl[k * N_full + k] = pik_full[k];
        }
    }

    if (N_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    int *valid_idx = (int *) R_alloc(N_valid, sizeof(int));
    double *pik_valid = (double *) R_alloc(N_valid, sizeof(double));

    int j = 0;
    double n_valid_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            valid_idx[j] = k;
            pik_valid[j] = pik_full[k];
            n_valid_sum += pik_full[k];
            j++;
        }
    }
    const int n_valid = (int)(n_valid_sum + 0.5);

    if (n_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    double *w = (double *) R_alloc(N_valid, sizeof(double));
    double *expa = (double *) R_alloc((size_t)N_valid * n_valid, sizeof(double));
    cps_calibrate(pik_valid, N_valid, n_valid, w, expa, eps, 100);

    for (int i = 0; i < N_valid; i++) {
        for (int jj = i + 1; jj < N_valid; jj++) {
            double r_i = w[i];
            double r_j = w[jj];
            double r_diff = r_j - r_i;
            double pi_ij;

            if (fabs(r_diff) > eps * (r_i + r_j + 1e-10)) {
                pi_ij = (r_j * pik_valid[i] - r_i * pik_valid[jj]) / r_diff;
            } else {
                double d = 0.0;
                for (int m = 0; m < N_valid; m++) {
                    d += pik_valid[m] * (1.0 - pik_valid[m]);
                }
                pi_ij = pik_valid[i] * pik_valid[jj] *
                        (1.0 - (1.0 - pik_valid[i]) * (1.0 - pik_valid[jj]) / d);
            }

            /* Clamp to valid range */
            if (pi_ij < 0.0) pi_ij = 0.0;
            if (pi_ij > pik_valid[i]) pi_ij = pik_valid[i];
            if (pi_ij > pik_valid[jj]) pi_ij = pik_valid[jj];

            int k = valid_idx[i];
            int l = valid_idx[jj];
            pikl[l * N_full + k] = pi_ij;  /* Column-major */
            pikl[k * N_full + l] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_up_systematic_jip(SEXP pik_sexp, SEXP eps_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);
    const double eps = REAL(eps_sexp)[0];

    SEXP result = PROTECT(allocMatrix(REALSXP, N_full, N_full));
    double *pikl = REAL(result);

    for (int i = 0; i < N_full; i++) {
        for (int j = 0; j < N_full; j++) {
            pikl[j * N_full + i] = pik_full[i] * pik_full[j];
        }
        pikl[i * N_full + i] = pik_full[i];
    }

    int N = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            N++;
        }
    }

    if (N < 2) {
        UNPROTECT(1);
        return result;
    }

    int *valid_idx = (int *) R_alloc(N, sizeof(int));
    double *pik1 = (double *) R_alloc(N, sizeof(double));

    int j = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            valid_idx[j] = k;
            pik1[j] = pik_full[k];
            j++;
        }
    }

    /* Cumulative sums */
    double *Vk = (double *) R_alloc(N, sizeof(double));
    double *Vk1 = (double *) R_alloc(N, sizeof(double));

    double cumsum = 0.0;
    for (int i = 0; i < N; i++) {
        cumsum += pik1[i];
        Vk[i] = cumsum;
        Vk1[i] = fmod(Vk[i], 1.0);
    }

    if (Vk1[N-1] < eps || Vk1[N-1] > 1.0 - eps) {
        Vk1[N-1] = 0.0;
    }

    double *sorted_Vk1 = (double *) R_alloc(N, sizeof(double));
    memcpy(sorted_Vk1, Vk1, N * sizeof(double));
    qsort(sorted_Vk1, N, sizeof(double), compare_double);

    double *r = (double *) R_alloc(N + 1, sizeof(double));
    memcpy(r, sorted_Vk1, N * sizeof(double));
    r[N] = 1.0;

    double *cent = (double *) R_alloc(N, sizeof(double));
    double *p = (double *) R_alloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        cent[i] = (r[i] + r[i + 1]) / 2.0;
        p[i] = r[i + 1] - r[i];
    }

    int *M = (int *) R_alloc((size_t)N * N, sizeof(int));

    for (int i = 0; i < N; i++) {
        for (int jj = 0; jj < N; jj++) {
            double A_curr = fmod(Vk[i] - cent[jj] + 10.0, 1.0);
            double A_prev = (i == 0) ? fmod(-cent[jj] + 10.0, 1.0)
                                     : fmod(Vk[i-1] - cent[jj] + 10.0, 1.0);
            M[i * N + jj] = (A_prev > A_curr) ? 1 : 0;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int jj = i + 1; jj < N; jj++) {
            double pi_ij = 0.0;
            for (int m = 0; m < N; m++) {
                pi_ij += p[m] * M[i * N + m] * M[jj * N + m];
            }

            int k = valid_idx[i];
            int l = valid_idx[jj];
            pikl[l * N_full + k] = pi_ij;
            pikl[k * N_full + l] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_up_brewer_jip(SEXP pik_sexp) {
    const int N = LENGTH(pik_sexp);
    const double *pik = REAL(pik_sexp);

    double n = 0.0;
    double sum_pik2 = 0.0;
    for (int k = 0; k < N; k++) {
        n += pik[k];
        sum_pik2 += pik[k] * pik[k];
    }

    SEXP result = PROTECT(allocMatrix(REALSXP, N, N));
    double *pikl = REAL(result);

    if (n <= 1.0 + 1e-10) {
        memset(pikl, 0, (size_t)N * N * sizeof(double));
        for (int k = 0; k < N; k++) {
            pikl[k * N + k] = pik[k];
        }
        UNPROTECT(1);
        return result;
    }

    double *c = (double *) R_alloc(N, sizeof(double));
    const double nm1 = n - 1.0;
    const double coef1 = (2.0 * n - 1.0) / nm1;
    const double coef2 = sum_pik2 / nm1;

    for (int k = 0; k < N; k++) {
        double denom = n - coef1 * pik[k] + coef2;
        c[k] = nm1 / denom;
    }

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (i == j) {
                pikl[i * N + i] = pik[i];
            } else {
                double pi_ij = pik[i] * pik[j] * (c[i] + c[j]) / 2.0;
                pikl[j * N + i] = pi_ij;
                pikl[i * N + j] = pi_ij;
            }
        }
    }

    UNPROTECT(1);
    return result;
}
