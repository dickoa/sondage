/*
 * cps_core.h - Shared routines for Conditional Poisson Sampling
 *
 * Elementary symmetric functions, conditional selection probabilities,
 * forward probability propagation, and Newton calibration.
 */

#ifndef CPS_CORE_H
#define CPS_CORE_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <math.h>

/* expa[i*n + z] = e_z(w_i, ..., w_{N-1}) */
static void cps_compute_expa(const double *w, int N, int n, double *expa) {
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

/* q(i,z) = P(select i | z still needed from i..N-1) */
static inline double cps_compute_q(const double *w, const double *expa,
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

    return w[i] * expa[(i + 1) * n + (z - 1)] / denom;
}

/* pro[i*n+z] = P(z needed at unit i), pik[i] = sum_z pro[i,z]*q(i,z) */
static void cps_compute_pik(const double *w, const double *expa,
                             int N, int n, double *pik, double *pro) {
    memset(pro, 0, (size_t)N * n * sizeof(double));

    pro[0 * n + (n - 1)] = 1.0;

    pik[0] = 0.0;
    for (int z = 0; z < n; z++) {
        double q_0z = cps_compute_q(w, expa, N, n, 0, z);
        pik[0] += pro[0 * n + z] * q_0z;
    }

    for (int i = 1; i < N; i++) {
        for (int z = n - 1; z >= 0; z--) {
            double q_prev = cps_compute_q(w, expa, N, n, i - 1, z);

            pro[i * n + z] += pro[(i - 1) * n + z] * (1.0 - q_prev);

            if (z + 1 < n) {
                double q_prev_zp1 = cps_compute_q(w, expa, N, n, i - 1, z + 1);
                pro[i * n + z] += pro[(i - 1) * n + (z + 1)] * q_prev_zp1;
            }
        }

        pik[i] = 0.0;
        for (int z = 0; z < n; z++) {
            double q_iz = cps_compute_q(w, expa, N, n, i, z);
            pik[i] += pro[i * n + z] * q_iz;
        }
    }
}

/* Calibrate w so CPS achieves pik_target. w[] and expa[] set on exit. */
static int cps_calibrate(const double *pik_target, int N, int n,
                          double *w, double *expa,
                          double eps, int max_iter) {
    double *piktilde = (double *) R_alloc(N, sizeof(double));
    double *pik_implied = (double *) R_alloc(N, sizeof(double));
    double *pro = (double *) R_alloc((size_t)N * n, sizeof(double));

    memcpy(piktilde, pik_target, N * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < N; i++) {
            double p = piktilde[i];
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            w[i] = p / (1.0 - p);
        }

        cps_compute_expa(w, N, n, expa);
        cps_compute_pik(w, expa, N, n, pik_implied, pro);

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
            cps_compute_expa(w, N, n, expa);
            return iter + 1;
        }
    }

    return max_iter;
}

static void cps_sample(const double *w, const double *expa,
                        int N, int n, int *sample_idx) {
    int remaining = n;
    int count = 0;

    for (int k = 0; k < N && remaining > 0; k++) {
        double q = cps_compute_q(w, expa, N, n, k, remaining - 1);
        if (unif_rand() < q) {
            sample_idx[count++] = k;
            remaining--;
        }
    }
}

#endif /* CPS_CORE_H */
