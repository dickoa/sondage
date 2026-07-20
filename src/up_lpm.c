/*
 * up_lpm.c: Local pivotal method 2 (Grafström, Lundström & Schelin, 2012)
 *
 * Spatially balanced (well-spread) sampling with prescribed inclusion
 * probabilities. Repeatedly, a randomly chosen undecided unit and its
 * nearest undecided neighbour in the spreading space compete in a
 * pivotal step (Deville & Tillé, 1998), which resolves at least one of
 * the two to 0 or 1 while preserving E(s) = pik. Nearby units thus
 * tend to exclude each other, spreading the sample.
 *
 * Implemented from the published algorithm descriptions:
 *   Grafström, A., Lundström, N.L.P. and Schelin, L. (2012). Spatially
 *     balanced sampling through the pivotal method. Biometrics, 68(2),
 *     514-520.
 *   Deville, J.-C. and Tillé, Y. (1998). Unequal probability sampling
 *     without replacement through a splitting method. Biometrika,
 *     85(1), 89-101.
 *
 * Nearest-neighbour ties are broken uniformly at random (reservoir
 * over tied candidates), which the design needs for E(s) = pik to
 * hold on gridded coordinates where exact distance ties are common.
 *
 * The neighbour search is a linear scan over undecided units, so a
 * draw is O(N^2 * d). Fine for typical design sizes; a k-d tree would
 * take this to O(N log N) if ever needed.
 */

#include "sampling_core.h"
#include "spatial_core.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <string.h>

/* Nearest undecided neighbour of id1 in the spreading space; exact
 * distance ties are broken uniformly at random. `z` is column-major
 * N x d. Requires pool->len >= 2. */
static int lpm_nearest(const SpatialPool *pool, const double *z, int N, int d,
                       int id1) {
    int best = -1;
    int n_ties = 0;
    double best_dist = R_PosInf;

    for (int k = 0; k < pool->len; k++) {
        int id = pool->list[k];
        if (id == id1) continue;

        double dist = 0.0;
        for (int c = 0; c < d; c++) {
            double diff = z[id + c * N] - z[id1 + c * N];
            dist += diff * diff;
            if (dist > best_dist) break;
        }

        if (dist < best_dist) {
            best_dist = dist;
            best = id;
            n_ties = 1;
        } else if (dist == best_dist) {
            n_ties++;
            if (unif_rand() * n_ties < 1.0) {
                best = id;
            }
        }
    }
    return best;
}

/* Move prob[id] out of the pool once it has numerically reached a
 * boundary; snaps to exactly 0 or 1 so extraction is unambiguous. */
static void lpm_settle(SpatialPool *pool, double *prob, int id, double eps) {
    if (prob[id] <= eps) {
        prob[id] = 0.0;
        spatial_pool_remove(pool, id);
    } else if (prob[id] >= 1.0 - eps) {
        prob[id] = 1.0;
        spatial_pool_remove(pool, id);
    }
}

/* One LPM2 draw over the working probabilities in `prob`. On return
 * every prob[i] is exactly 0 or 1. */
static void lpm2_run(double *prob, const double *z, int N, int d,
                     double eps, SpatialPool *pool) {
    spatial_pool_fill(pool, prob, N);

    int steps = 0;
    while (pool->len > 1) {
        if ((++steps & 1023) == 0) R_CheckUserInterrupt();

        int k1 = (int)(unif_rand() * pool->len);
        if (k1 >= pool->len) k1 = pool->len - 1;  /* Safety clamp */
        int id1 = pool->list[k1];
        int id2 = lpm_nearest(pool, z, N, d, id1);

        /* Pivotal step (Deville & Tillé, 1998): the pair's total
         * probability mass is preserved, one unit reaches 0 or 1. */
        double p1 = prob[id1];
        double p2 = prob[id2];
        double psum = p1 + p2;

        if (psum <= 1.0) {
            /* One unit gives its mass to the other. */
            if (unif_rand() * psum < p1) {
                prob[id1] = psum;
                prob[id2] = 0.0;
            } else {
                prob[id2] = psum;
                prob[id1] = 0.0;
            }
        } else {
            /* One unit is selected, the other keeps the excess. */
            double prest = psum - 1.0;
            if (unif_rand() * (2.0 - psum) < 1.0 - p2) {
                prob[id1] = 1.0;
                prob[id2] = prest;
            } else {
                prob[id2] = 1.0;
                prob[id1] = prest;
            }
        }

        lpm_settle(pool, prob, id1, eps);
        lpm_settle(pool, prob, id2, eps);
    }

    /* A leftover unit only carries the rounding residue of an integer
     * pik total; a Bernoulli draw is exact for any residual value. */
    if (pool->len == 1) {
        int id = pool->list[0];
        prob[id] = (unif_rand() < prob[id]) ? 1.0 : 0.0;
        spatial_pool_remove(pool, id);
    }
}

/* --- R-callable entry points -------------------------------------------- */

SEXP C_lpm2(SEXP pik_sexp, SEXP z_sexp, SEXP eps_sexp) {
    const int N = length(pik_sexp);
    const int d = spatial_check_spread(z_sexp, N);

    const double eps = asReal(eps_sexp);
    const double *z = REAL(z_sexp);
    double *prob = (double *)R_alloc(N, sizeof(double));
    int *selected = (int *)R_alloc(N, sizeof(int));
    memcpy(prob, REAL(pik_sexp), (size_t)N * sizeof(double));

    SpatialPool *pool = spatial_pool_alloc(N);

    GetRNGstate();
    lpm2_run(prob, z, N, d, eps, pool);
    PutRNGstate();

    int n_selected = sampling_extract_selected(prob, N, 0.5, selected, N);

    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
    UNPROTECT(1);
    return result;
}

SEXP C_lpm2_batch(SEXP pik_sexp, SEXP z_sexp, SEXP eps_sexp,
                  SEXP nrep_sexp) {
    const int N = length(pik_sexp);
    const int d = spatial_check_spread(z_sexp, N);

    const int nrep = INTEGER(nrep_sexp)[0];
    const double eps = asReal(eps_sexp);
    const double *pik = REAL(pik_sexp);
    const double *z = REAL(z_sexp);

    double n_sum = 0.0;
    for (int i = 0; i < N; i++) n_sum += pik[i];
    const int n_target = (int)(n_sum + 0.5);

    SEXP result = PROTECT(allocMatrix(INTSXP, n_target, nrep));
    int *res = INTEGER(result);

    if (nrep <= 0) {
        UNPROTECT(1);
        return result;
    }

    double *prob = (double *)R_alloc(N, sizeof(double));
    SpatialPool *pool = spatial_pool_alloc(N);

    GetRNGstate();
    for (int s = 0; s < nrep; s++) {
        if (s % 32 == 0) R_CheckUserInterrupt();

        memcpy(prob, pik, (size_t)N * sizeof(double));
        lpm2_run(prob, z, N, d, eps, pool);

        int *col = res + (R_xlen_t)s * n_target;
        int n_selected = sampling_extract_selected(
            prob, N, 0.5, col, n_target
        );
        if (n_selected != n_target) {
            PutRNGstate();
            UNPROTECT(1);
            error(
                "lpm2 batch draw %d produced size %d, expected %d",
                s + 1,
                n_selected,
                n_target
            );
        }
    }
    PutRNGstate();

    UNPROTECT(1);
    return result;
}
