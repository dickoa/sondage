/*
 * joint_probs.c - Joint Inclusion Probabilities
 *
 * Implements:
 * - C_up_maxent_joint: Exact CPS joint probabilities (Aires' formula)
 * - C_up_systematic_joint: Exact systematic sampling joint probabilities
 * - C_up_brewer_joint: Brewer approximation (equation 18, Brewer & Donadio 2003)
 *
 * References:
 * - Aires (1999). Algorithms to Find Exact Inclusion Probabilities for Conditional Poisson Sampling and Pareto πps Sampling Design
 * - Tillé (2006). Sampling Algorithms. Springer.
 * - Brewer & Donadio (2003). The High Entropy Variance of the HT Estimator.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <float.h>

static int compare_double(const void *a, const void *b) {
    double diff = *(const double *)a - *(const double *)b;
    return (diff > 0) - (diff < 0);
}

static void compute_expa_joint(const double *w, int N, int n, double *expa) {
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

static inline double compute_q_joint(const double *w, const double *expa,
                                      int N, int n, int i, int z) {
    if (z == 0) {
        double denom = expa[i * n + 0];
        return (denom > 0) ? w[i] / denom : 0.0;
    }
    if (z >= N - i - 1) return 1.0;

    double denom = expa[i * n + z];
    if (denom <= 0) return 0.0;

    return w[i] * expa[(i + 1) * n + (z - 1)] / denom;
}

static void compute_pik_from_w(const double *w, const double *expa,
                                int N, int n, double *pik) {
    double *pro = (double *) R_alloc((size_t)N * n, sizeof(double));
    memset(pro, 0, (size_t)N * n * sizeof(double));

    pro[0 * n + (n - 1)] = 1.0;

    pik[0] = 0.0;
    for (int z = 0; z < n; z++) {
        double q_0z = compute_q_joint(w, expa, N, n, 0, z);
        pik[0] += pro[0 * n + z] * q_0z;
    }

    for (int i = 1; i < N; i++) {
        for (int z = n - 1; z >= 0; z--) {
            double q_prev = compute_q_joint(w, expa, N, n, i - 1, z);
            pro[i * n + z] += pro[(i - 1) * n + z] * (1.0 - q_prev);

            if (z + 1 < n) {
                double q_prev_zp1 = compute_q_joint(w, expa, N, n, i - 1, z + 1);
                pro[i * n + z] += pro[(i - 1) * n + (z + 1)] * q_prev_zp1;
            }
        }

        pik[i] = 0.0;
        for (int z = 0; z < n; z++) {
            double q_iz = compute_q_joint(w, expa, N, n, i, z);
            pik[i] += pro[i * n + z] * q_iz;
        }
    }
}

static void calibrate_w_from_pik(const double *pik_target, int N, int n,
                                  double *w, double eps, int max_iter) {
    double *piktilde = (double *) R_alloc(N, sizeof(double));
    double *pik_implied = (double *) R_alloc(N, sizeof(double));
    double *expa = (double *) R_alloc((size_t)N * n, sizeof(double));

    memcpy(piktilde, pik_target, N * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < N; i++) {
            double p = piktilde[i];
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            w[i] = p / (1.0 - p);
        }

        compute_expa_joint(w, N, n, expa);
        compute_pik_from_w(w, expa, N, n, pik_implied);

        double max_diff = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = pik_target[i] - pik_implied[i];
            piktilde[i] += diff;
            if (piktilde[i] < 1e-10) piktilde[i] = 1e-10;
            if (piktilde[i] > 1.0 - 1e-10) piktilde[i] = 1.0 - 1e-10;
            if (fabs(diff) > max_diff) max_diff = fabs(diff);
        }

        if (max_diff < eps) break;
    }

    for (int i = 0; i < N; i++) {
        double p = piktilde[i];
        if (p < 1e-10) p = 1e-10;
        if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
        w[i] = p / (1.0 - p);
    }
}

SEXP C_up_maxent_joint(SEXP pik_sexp, SEXP eps_sexp) {
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
    calibrate_w_from_pik(pik_valid, N_valid, n_valid, w, eps, 100);

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

SEXP C_up_systematic_joint(SEXP pik_sexp, SEXP eps_sexp) {
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

    int *M = (int *) R_alloc(N * N, sizeof(int));

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

SEXP C_up_brewer_joint(SEXP pik_sexp) {
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
