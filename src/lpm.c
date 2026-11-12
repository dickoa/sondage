/*
 * up_lpm.c: Local pivotal method 2 (Grafström, Lundström & Schelin, 2012)
 *
 * Spatially balanced (well-spread) sampling with prescribed inclusion
 * probabilities. Repeatedly, a randomly chosen undecided unit and its
 * nearest undecided neighbor in the spreading space compete in a
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
 * Nearest-neighbor ties are broken uniformly at random (reservoir
 * over tied candidates), which the design needs for E(s) = pik to
 * hold on gridded coordinates where exact distance ties are common.
 *
 * With at most LPM_KDTREE_MAX_DIM spreading variables the neighbor
 * search uses the static k-d tree of kdtree_core.h, so a draw costs
 * about O(N log N). Above that the tree prunes too little and a linear
 * scan over undecided units is used, O(N^2 * d) per draw. Both searches
 * find the same set of nearest units, and with no exact distance ties
 * neither draws a random number. With ties, the tree picks one of the
 * tied units with a single uniform draw, and the scan keeps a reservoir
 * in pool order, so a seed can select differently on the two paths,
 * always uniformly among the nearest units. Heavily tied layouts, such
 * as duplicated points, make the tree visit every tied unit, so a draw
 * there approaches O(N^2).
 */

#include "sampling_core.h"
#include "spatial_core.h"
#include "kdtree_core.h"
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <string.h>

#define LPM_KDTREE_MAX_DIM 8

/* Nearest undecided neighbor of id1 in the spreading space by linear
 * scan; exact distance ties are broken uniformly at random. `z` is
 * column-major N x d. Requires pool->len >= 2. */
static int lpm_nearest(const SpatialPool *pool, const double *z, int N, int d,
                       int id1) {
  int best = -1;
  int n_ties = 0;
  double best_dist = R_PosInf;

  for (int k = 0; k < pool->len; k++) {
    int id = pool->list[k];
    if (id == id1)
      continue;

    double dist = kd_dist(z, N, d, id, id1, best_dist);

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

/* Nearest undecided neighbor of id1 through the k-d tree. Tied units
 * are resolved with one uniform draw over the tied set. Any order of the
 * ties gives the same uniform choice, so they are not sorted: sorting
 * them by pool position was quadratic in the number of ties at each
 * pivotal step, which duplicated points make close to N. */
static int lpm_nearest_tree(KdTree *tree, const SpatialPool *pool, int id1) {
  int n_ties = kd_nearest(tree, id1);
  int *ties = tree->ties;

#ifdef SONDAGE_KDTREE_CHECK
  {
    /* Debug builds only: the tree must find exactly the units at the
     * minimum distance that a full scan finds. */
    double best = R_PosInf;
    int n_scan = 0;
    for (int k = 0; k < pool->len; k++) {
      int id = pool->list[k];
      if (id == id1)
        continue;
      double dist = kd_dist(tree->z, tree->N, tree->d, id, id1, R_PosInf);
      if (dist < best) {
        best = dist;
        n_scan = 1;
      } else if (dist == best) {
        n_scan++;
      }
    }
    int found = 0;
    for (int i = 0; i < n_ties; i++) {
      double dist =
          kd_dist(tree->z, tree->N, tree->d, ties[i], id1, R_PosInf);
      if (dist == best && pool->reverse[ties[i]] < pool->len)
        found++;
    }
    if (found != n_scan || n_ties != n_scan)
      error("k-d tree check failed for unit %d: tree %d ties, scan %d",
            id1 + 1, n_ties, n_scan);
  }
#endif

  (void)pool; /* read only by the debug check above */

  if (n_ties == 1)
    return ties[0];

  int k = (int)(unif_rand() * n_ties);
  if (k >= n_ties)
    k = n_ties - 1;
  return ties[k];
}

/* Move prob[id] out of the pool once it has numerically reached a
 * boundary; snaps to exactly 0 or 1 so extraction is unambiguous. */
static void lpm_settle(SpatialPool *pool, KdTree *tree, double *prob, int id,
                       double eps) {
  if (prob[id] <= eps) {
    prob[id] = 0.0;
  } else if (prob[id] >= 1.0 - eps) {
    prob[id] = 1.0;
  } else {
    return;
  }
  spatial_pool_remove(pool, id);
  if (tree)
    kd_remove(tree, id);
}

/* One LPM2 draw over the working probabilities in `prob`. On return
 * every prob[i] is exactly 0 or 1. `tree`, when not NULL, was built from
 * the same starting probabilities and is reset here. */
static void lpm2_run(double *prob, const double *z, int N, int d, double eps,
                     SpatialPool *pool, KdTree *tree) {
  spatial_pool_fill(pool, prob, N);
  if (tree)
    kd_reset(tree);

  int steps = 0;
  while (pool->len > 1) {
    if ((++steps & 1023) == 0)
      R_CheckUserInterrupt();

    int k1 = (int)(unif_rand() * pool->len);
    if (k1 >= pool->len)
      k1 = pool->len - 1; /* Safety clamp */
    int id1 = pool->list[k1];
    int id2 = tree ? lpm_nearest_tree(tree, pool, id1)
                   : lpm_nearest(pool, z, N, d, id1);

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

    lpm_settle(pool, tree, prob, id1, eps);
    lpm_settle(pool, tree, prob, id2, eps);
  }

  /* A leftover unit only carries the rounding residue of an integer
   * pik total; a Bernoulli draw is exact for any residual value. */
  if (pool->len == 1) {
    int id = pool->list[0];
    prob[id] = (unif_rand() < prob[id]) ? 1.0 : 0.0;
    spatial_pool_remove(pool, id);
  }
}

SEXP C_lpm2(SEXP pik_sexp, SEXP z_sexp, SEXP eps_sexp) {
  const int N = length(pik_sexp);
  const int d = spatial_check_spread(z_sexp, N);

  const double eps = asReal(eps_sexp);
  const double *z = REAL(z_sexp);
  double *prob = (double *)R_alloc(N, sizeof(double));
  int *selected = (int *)R_alloc(N, sizeof(int));
  memcpy(prob, REAL(pik_sexp), (size_t)N * sizeof(double));

  SpatialPool *pool = spatial_pool_alloc(N);
  KdTree *tree = (d <= LPM_KDTREE_MAX_DIM) ? kd_build(z, N, d, prob) : NULL;

  GetRNGstate();
  lpm2_run(prob, z, N, d, eps, pool, tree);
  PutRNGstate();

  int n_selected = sampling_extract_selected(prob, N, 0.5, selected, N);

  SEXP result = PROTECT(allocVector(INTSXP, n_selected));
  memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
  UNPROTECT(1);
  return result;
}

SEXP C_lpm2_batch(SEXP pik_sexp, SEXP z_sexp, SEXP eps_sexp, SEXP nrep_sexp) {
  const int N = length(pik_sexp);
  const int d = spatial_check_spread(z_sexp, N);

  const int nrep = INTEGER(nrep_sexp)[0];
  const double eps = asReal(eps_sexp);
  const double *pik = REAL(pik_sexp);
  const double *z = REAL(z_sexp);

  double n_sum = 0.0;
  for (int i = 0; i < N; i++)
    n_sum += pik[i];
  const int n_target = (int)(n_sum + 0.5);

  SEXP result = PROTECT(allocMatrix(INTSXP, n_target, nrep));
  int *res = INTEGER(result);

  if (nrep <= 0) {
    UNPROTECT(1);
    return result;
  }

  double *prob = (double *)R_alloc(N, sizeof(double));
  SpatialPool *pool = spatial_pool_alloc(N);
  KdTree *tree = (d <= LPM_KDTREE_MAX_DIM) ? kd_build(z, N, d, pik) : NULL;

  GetRNGstate();
  for (int s = 0; s < nrep; s++) {
    if (s % 32 == 0)
      R_CheckUserInterrupt();

    memcpy(prob, pik, (size_t)N * sizeof(double));
    lpm2_run(prob, z, N, d, eps, pool, tree);

    int *col = res + (R_xlen_t)s * n_target;
    int n_selected = sampling_extract_selected(prob, N, 0.5, col, n_target);
    if (n_selected != n_target) {
      PutRNGstate();
      UNPROTECT(1);
      error("lpm2 batch draw %d produced size %d, expected %d", s + 1,
            n_selected, n_target);
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return result;
}
