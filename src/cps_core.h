/*
 * cps_core.h - Shared routines for Conditional Poisson Sampling
 *
 * Poisson-binomial suffix probabilities, conditional selection
 * probabilities, forward probability propagation, and fixed-point
 * calibration.
 *
 * All quantities are kept in the probability domain. The previous
 * implementation stored raw elementary symmetric functions, which grow
 * like binomial coefficients and overflow to Inf near N ~ 1500 for
 * central sample sizes; the resulting non-finite conditionals were
 * clamped to zero, severely biasing selection toward later units.
 */

#ifndef CPS_CORE_H
#define CPS_CORE_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Number of doubles needed for the suffix probability table. */
#define CPS_TABLE_SIZE(N, n) ((size_t)((N) + 1) * (size_t)((n) + 1))

/*
 * f[i*(n+1) + z] = P(exactly z of units i, ..., N-1 are selected)
 * under independent Bernoulli(p_k) with p_k = w_k / (1 + w_k).
 * Rows run i = 0..N (row N is the empty-suffix base case); only
 * z = 0..min(n, N-i) can be nonzero. Every entry lies in [0, 1].
 *
 * The elementary symmetric functions relate to this table through
 * e_z(w_i, ...) = f[i][z] * prod_{k>=i} (1 + w_k), so conditional
 * probability ratios of ESFs can be computed from f directly without
 * ever forming the (overflowing) ESFs themselves.
 */
static void cps_compute_f(const double *w, int N, int n, double *f) {
    const int S = n + 1; /* row stride */
    memset(f, 0, CPS_TABLE_SIZE(N, n) * sizeof(double));

    f[N * S + 0] = 1.0;
    for (int i = N - 1; i >= 0; i--) {
        const double p = w[i] / (1.0 + w[i]);
        const double *fnext = f + (size_t)(i + 1) * S;
        double *fi = f + (size_t)i * S;
        int max_z = N - i;
        if (max_z > n) max_z = n;

        fi[0] = (1.0 - p) * fnext[0];
        for (int z = 1; z <= max_z; z++) {
            fi[z] = (1.0 - p) * fnext[z] + p * fnext[z - 1];
        }
    }
}

/* q(i,z) = P(select i | z+1 still needed from i..N-1) */
static inline double cps_compute_q(const double *w, const double *f,
                                    int N, int n, int i, int z) {
    const int S = n + 1;
    const int r = z + 1; /* selections still needed */

    if (r >= N - i) {
        return 1.0; /* every remaining unit must be selected */
    }

    double denom = f[(size_t)i * S + r];
    if (!(denom > 0.0)) return 0.0;

    double p = w[i] / (1.0 + w[i]);
    double q = p * f[(size_t)(i + 1) * S + (r - 1)] / denom;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    return q;
}

/* pro[i*n+z] = P(z+1 needed at unit i), pik[i] = sum_z pro[i,z]*q(i,z) */
static void cps_compute_pik(const double *w, const double *f,
                             int N, int n, double *pik, double *pro) {
    memset(pro, 0, (size_t)N * n * sizeof(double));

    pro[0 * n + (n - 1)] = 1.0;

    pik[0] = 0.0;
    for (int z = 0; z < n; z++) {
        double q_0z = cps_compute_q(w, f, N, n, 0, z);
        pik[0] += pro[0 * n + z] * q_0z;
    }

    for (int i = 1; i < N; i++) {
        for (int z = n - 1; z >= 0; z--) {
            double q_prev = cps_compute_q(w, f, N, n, i - 1, z);

            pro[i * n + z] += pro[(i - 1) * n + z] * (1.0 - q_prev);

            if (z + 1 < n) {
                double q_prev_zp1 = cps_compute_q(w, f, N, n, i - 1, z + 1);
                pro[i * n + z] += pro[(i - 1) * n + (z + 1)] * q_prev_zp1;
            }
        }

        pik[i] = 0.0;
        for (int z = 0; z < n; z++) {
            double q_iz = cps_compute_q(w, f, N, n, i, z);
            pik[i] += pro[i * n + z] * q_iz;
        }
    }
}

/*
 * Calibrate w so CPS achieves pik_target.
 *
 * Outputs:
 *   w[], f[]         - calibrated odds and suffix probability table
 *                      (set on exit; f has CPS_TABLE_SIZE(N, n) entries)
 *   *final_max_diff  - max |pik_target - pik_implied| at stop (0 on fast path)
 *   *worst_idx       - index in pik_target where the max was achieved
 *                      (-1 on fast path)
 *
 * Returns number of iterations used (1 for fast path, up to max_iter).
 *
 * The iteration is a fixed-point update: piktilde <- piktilde + diff
 * until max|diff| < eps. For inclusion probabilities close to 1 the
 * iteration asymptotes at a non-zero defect on the order of (1 - max_pik);
 * the final_max_diff / worst_idx outputs let callers emit an informative
 * warning when the target tolerance isn't reached.
 */
static int cps_calibrate(const double *pik_target, int N, int n,
                          double *w, double *f,
                          double eps, int max_iter,
                          double *final_max_diff, int *worst_idx) {
    if (final_max_diff) *final_max_diff = 0.0;
    if (worst_idx) *worst_idx = -1;

    if (N <= 0 || n <= 0) return 0;

    /*
     * Exact fast path: if all target probabilities are equal, calibration
     * has a closed-form solution with common odds w = p / (1 - p).
     */
    {
        int all_equal = 1;
        double p0 = pik_target[0];
        double tol = 16.0 * DBL_EPSILON * fmax(1.0, fabs(p0));
        for (int i = 1; i < N; i++) {
            if (fabs(pik_target[i] - p0) > tol) {
                all_equal = 0;
                break;
            }
        }
        if (all_equal) {
            double p = p0;
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            double w0 = p / (1.0 - p);
            for (int i = 0; i < N; i++) {
                w[i] = w0;
            }
            cps_compute_f(w, N, n, f);
            return 1;
        }
    }

    double *piktilde = (double *) R_alloc(N, sizeof(double));
    double *pik_implied = (double *) R_alloc(N, sizeof(double));
    double *pro = (double *) R_alloc((size_t)N * n, sizeof(double));

    memcpy(piktilde, pik_target, N * sizeof(double));

    double last_max_diff = 0.0;
    int    last_worst_idx = -1;

    for (int iter = 0; iter < max_iter; iter++) {
        R_CheckUserInterrupt();
        for (int i = 0; i < N; i++) {
            double p = piktilde[i];
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            w[i] = p / (1.0 - p);
        }

        cps_compute_f(w, N, n, f);
        cps_compute_pik(w, f, N, n, pik_implied, pro);

        double max_diff = 0.0;
        int    max_idx = -1;
        for (int i = 0; i < N; i++) {
            double diff = pik_target[i] - pik_implied[i];
            piktilde[i] += diff;

            if (piktilde[i] < 1e-10) piktilde[i] = 1e-10;
            if (piktilde[i] > 1.0 - 1e-10) piktilde[i] = 1.0 - 1e-10;

            if (fabs(diff) > max_diff) {
                max_diff = fabs(diff);
                max_idx = i;
            }
        }

        last_max_diff = max_diff;
        last_worst_idx = max_idx;

        if (max_diff < eps) {
            for (int i = 0; i < N; i++) {
                double p = piktilde[i];
                if (p < 1e-10) p = 1e-10;
                if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
                w[i] = p / (1.0 - p);
            }
            cps_compute_f(w, N, n, f);
            if (final_max_diff) *final_max_diff = max_diff;
            if (worst_idx) *worst_idx = max_idx;
            return iter + 1;
        }
    }

    if (final_max_diff) *final_max_diff = last_max_diff;
    if (worst_idx) *worst_idx = last_worst_idx;
    return max_iter;
}

/*
 * Emit an informative warning when cps_calibrate did not reach the target
 * tolerance. The message reports the final max_diff and, when available,
 * counts and flags the 1-based population indices of units whose targets
 * are within boundary_thresh of 0 or 1 (those are the usual culprits).
 *
 * idx_map maps the valid-subset index to the 0-based population index;
 * N_valid is the size of pik_valid.
 */
static void cps_warn_nonconverge(const char *ctx,
                                 int max_iter, double tol, double max_diff,
                                 const double *pik_valid, int N_valid,
                                 const int *idx_map) {
    const double boundary_thresh = 1e-3;
    int near_boundary = 0;
    int first_boundary = -1;
    if (pik_valid && idx_map && N_valid > 0) {
        for (int i = 0; i < N_valid; i++) {
            double p = pik_valid[i];
            if (p < boundary_thresh || p > 1.0 - boundary_thresh) {
                if (first_boundary < 0) first_boundary = idx_map[i] + 1;
                near_boundary++;
            }
        }
    }

    if (near_boundary > 0) {
        Rf_warning(
            "%sCPS calibration did not reach tolerance after %d iterations "
            "(max_diff = %.2e, tolerance = %.2e). %d unit(s) have pik "
            "within %.0e of 0 or 1 (first at index %d); this is the usual "
            "cause. Realized first-order inclusion probabilities will "
            "differ from the target by up to max_diff. To silence this, "
            "clip pik away from the boundary before calling.",
            ctx, max_iter, max_diff, tol,
            near_boundary, boundary_thresh, first_boundary
        );
    } else {
        Rf_warning(
            "%sCPS calibration did not reach tolerance after %d iterations "
            "(max_diff = %.2e, tolerance = %.2e). Realized first-order "
            "inclusion probabilities will differ from the target by up to "
            "max_diff.",
            ctx, max_iter, max_diff, tol
        );
    }
}

/*
 * Draw one CPS sample.
 * Returns number of selected indices written to sample_idx.
 */
static int cps_sample(const double *w, const double *f,
                       int N, int n, int *sample_idx) {
    int remaining = n;
    int count = 0;

    for (int k = 0; k < N && remaining > 0; k++) {
        int units_left = N - k;

        /*
         * Numerical guard: if remaining units equals remaining selections,
         * all tail units must be selected.
         */
        if (units_left == remaining) {
            for (int kk = k; kk < N; kk++) {
                sample_idx[count++] = kk;
            }
            remaining = 0;
            break;
        }

        double q = cps_compute_q(w, f, N, n, k, remaining - 1);

        /*
         * Clamp/repair tiny numerical deviations. This prevents NaN/Inf from
         * silently producing invalid selection paths.
         */
        if (!R_FINITE(q)) {
            q = 0.0;
        } else if (q < 0.0) {
            q = 0.0;
        } else if (q > 1.0) {
            q = 1.0;
        }

        if (unif_rand() < q) {
            sample_idx[count++] = k;
            remaining--;
        }
    }

    return count;
}

#endif /* CPS_CORE_H */
