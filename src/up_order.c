/*
 * C implementations for systematic PPS, Poisson, SPS, and Pareto sampling.
 *
 * These replace the previous pure-R implementations for performance.
 * All functions select certainty units (pik exactly 1) unconditionally,
 * exclude impossible units (pik exactly 0), and run the algorithm on the
 * remaining units. Values close to but not exactly 0 or 1 are sampled as
 * given; callers validate that the pik sum is integer to floating-point
 * accuracy for fixed-size methods.
 */

#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* helpers for order sampling (SPS / Pareto) */

static void swap_pair(double *keys, int *orig, int a, int b) {
  double tk = keys[a];
  keys[a] = keys[b];
  keys[b] = tk;
  int ti = orig[a];
  orig[a] = orig[b];
  orig[b] = ti;
}

/*
 * Quickselect on orig[] within [band_lo, band_hi] (whose keys are all
 * equal): move the k smallest population indices to the front of the
 * band. orig values are unique, so this always terminates.
 */
static void tie_select_by_index(double *keys, int *orig,
                                int band_lo, int band_hi, int k) {
  if (k <= 0 || k > band_hi - band_lo) return; /* empty or full band */

  int lo = band_lo, hi = band_hi;
  const int want = band_lo + k; /* exclusive selection boundary */
  while (lo < hi) {
    int pivot = orig[lo + (hi - lo) / 2];
    int lt = lo, i = lo, gt = hi;
    while (i <= gt) {
      if (orig[i] < pivot)
        swap_pair(keys, orig, i++, lt++);
      else if (orig[i] > pivot)
        swap_pair(keys, orig, i, gt--);
      else
        i++;
    }
    if (want <= lt)
      hi = lt - 1;
    else if (want > gt + 1)
      lo = gt + 1;
    else
      return; /* boundary inside the (singleton) pivot position */
  }
}

/*
 * Partial selection on parallel arrays: rearranges keys[], orig[] so
 * that the n_select entries with the smallest keys are in positions
 * [0, n_select). Three-way quickselect: elements equal to the pivot are
 * grouped in a middle band, so fully or heavily tied keys partition in
 * linear time (a two-way Lomuto partition degrades to O(len^2) when
 * many keys are identical). Expected O(len) overall.
 *
 * Ties that straddle the selection boundary are broken toward the
 * smallest population index: deterministic and reproducible for
 * user-supplied permanent random numbers without consuming RNG values.
 *
 * Precondition: 0 <= n_select <= len.
 * When n_select == len every entry qualifies; the arrays are left
 * untouched and the caller is responsible for any ordering it needs.
 */
static void partial_sort_paired(double *keys, int *orig, int len,
                                int n_select) {
  if (n_select <= 0 || len <= 0) return;   /* nothing to do */
  if (n_select >= len) return;              /* all entries qualify */

  int lo = 0, hi = len - 1;
  while (lo < hi) {
    double pivot = keys[lo + (hi - lo) / 2];

    /* Dutch-national-flag partition:
     * [lo, lt) < pivot, [lt, gt] == pivot, (gt, hi] > pivot */
    int lt = lo, i = lo, gt = hi;
    while (i <= gt) {
      if (keys[i] < pivot)
        swap_pair(keys, orig, i++, lt++);
      else if (keys[i] > pivot)
        swap_pair(keys, orig, i, gt--);
      else
        i++;
    }

    if (n_select <= lt) {
      hi = lt - 1;
    } else if (n_select > gt + 1) {
      lo = gt + 1;
    } else {
      /* Boundary falls inside the tie band: deterministically keep the
       * tied units with the smallest population indices. */
      tie_select_by_index(keys, orig, lt, gt, n_select - lt);
      return;
    }
  }
}

/* Systematic PPS */

/*
 * C_up_systematic(pik)
 *
 * Single-pass systematic PPS: accumulate pik, select unit k when
 * the cumulative sum crosses an integer boundary (shifted by u).
 * Returns sorted 1-based indices.
 */
SEXP C_up_systematic(SEXP pik) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);

  /* Count certainty units and compute n */
  int n_certain = 0;
  double n_valid = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      n_certain++;
    else if (pk[i] > 0.0)
      n_valid += pk[i];
  }
  int n_select = (int)(n_valid + 0.5);
  int total = n_select + n_certain;

  SEXP result = PROTECT(allocVector(INTSXP, total));
  int *res = INTEGER(result);
  int j = 0;

  /* Certainty units first */
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    /*
     * Normalize the cumulative walk so it ends at exactly n_select:
     * callers guarantee |n_valid - n_select| is at floating-point scale,
     * and without this the walk could cross one extra (or one fewer)
     * integer boundary than there are slots.
     */
    const double scale = (double)n_select / n_valid;

    GetRNGstate();
    double u = unif_rand();
    PutRNGstate();

    double cs = 0.0;
    for (int i = 0; i < N; i++) {
      if (pk[i] <= 0.0 || pk[i] >= 1.0)
        continue;
      double prev = cs;
      cs += pk[i] * scale;
      if (floor(cs - u) > floor(prev - u)) {
        if (j >= total) {
          UNPROTECT(1);
          error(
            "systematic sampling produced more than the expected %d selections",
            total
          );
        }
        res[j++] = i + 1;
      }
    }
  }

  if (j != total) {
    UNPROTECT(1);
    error("systematic sampling produced %d selections, expected %d", j, total);
  }

  R_isort(res, total);

  UNPROTECT(1);
  return result;
}

/* Poisson sampling */

/*
 * C_up_poisson(pik, prn)
 *
 * prn is R_NilValue when NULL (generate uniforms internally).
 * Returns sorted 1-based indices.
 */
SEXP C_up_poisson(SEXP pik, SEXP prn) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  int has_prn = !isNull(prn);
  const double *u = has_prn ? REAL(prn) : NULL;

  /* First pass: count selections */
  int count = 0;
  int *sel = (int *)R_alloc(N, sizeof(int));

  if (has_prn) {
    for (int i = 0; i < N; i++) {
      if (u[i] < pk[i])
        sel[count++] = i + 1;
    }
  } else {
    GetRNGstate();
    for (int i = 0; i < N; i++) {
      if (unif_rand() < pk[i])
        sel[count++] = i + 1;
    }
    PutRNGstate();
  }

  SEXP result = PROTECT(allocVector(INTSXP, count));
  int *res = INTEGER(result);
  for (int i = 0; i < count; i++)
    res[i] = sel[i];

  UNPROTECT(1);
  return result;
}

/* SPS (Sequential Poisson Sampling) */

/*
 * C_up_sps(pik, prn)
 *
 * Key: xi_k = u_k / pik_k.  Select the n smallest keys among valid units.
 * prn is R_NilValue when NULL (generate uniforms internally).
 * Returns sorted 1-based indices.
 */
SEXP C_up_sps(SEXP pik, SEXP prn) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  int has_prn = !isNull(prn);
  const double *u_ext = has_prn ? REAL(prn) : NULL;

  /* Separate certainty vs valid units */
  int n_certain = 0;
  double n_sum = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      n_certain++;
    else if (pk[i] > 0.0)
      n_sum += pk[i];
  }
  int n = (int)(n_sum + n_certain + 0.5);
  int n_select = n - n_certain;

  int total = n;
  SEXP result = PROTECT(allocVector(INTSXP, total));
  int *res = INTEGER(result);
  int j = 0;

  /* Certainty units */
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    /* Build key array for valid units */
    int n_valid = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > 0.0 && pk[i] < 1.0)
        n_valid++;
    }

    if (n_select > n_valid) {
      UNPROTECT(1);
      error(
        "sps sampling needs %d units but only %d have 0 < pik < 1",
        n_select, n_valid
      );
    }

    double *keys = (double *)R_alloc(n_valid, sizeof(double));
    int *orig = (int *)R_alloc(n_valid, sizeof(int));

    if (!has_prn)
      GetRNGstate();

    int vi = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > 0.0 && pk[i] < 1.0) {
        double u_i = has_prn ? u_ext[i] : unif_rand();
        keys[vi] = u_i / pk[i];
        orig[vi] = i + 1; /* 1-based */
        vi++;
      }
    }

    if (!has_prn)
      PutRNGstate();

    /* Select the n_select smallest keys. When n_select == n_valid the
     * partition is a no-op; the loop below reads all n_valid entries and
     * the outer qsort(res, ...) restores ascending-index order. */
    partial_sort_paired(keys, orig, n_valid, n_select);

    for (int i = 0; i < n_select; i++) {
      res[j++] = orig[i];
    }
  }

  /* Sort full result (matches up_chromy.c idiom) */
  R_isort(res, total);

  UNPROTECT(1);
  return result;
}

/* Pareto sampling */

/*
 * C_up_pareto(pik, prn)
 *
 * Key: xi_k = [u_k/(1-u_k)] / [pik_k/(1-pik_k)].
 * Select the n smallest keys among valid units.
 * prn is R_NilValue when NULL (generate uniforms internally).
 * Returns sorted 1-based indices.
 */
SEXP C_up_pareto(SEXP pik, SEXP prn) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  int has_prn = !isNull(prn);
  const double *u_ext = has_prn ? REAL(prn) : NULL;

  /* Separate certainty vs valid units */
  int n_certain = 0;
  double n_sum = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      n_certain++;
    else if (pk[i] > 0.0)
      n_sum += pk[i];
  }
  int n = (int)(n_sum + n_certain + 0.5);
  int n_select = n - n_certain;

  int total = n;
  SEXP result = PROTECT(allocVector(INTSXP, total));
  int *res = INTEGER(result);
  int j = 0;

  /* Certainty units */
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    int n_valid = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > 0.0 && pk[i] < 1.0)
        n_valid++;
    }

    if (n_select > n_valid) {
      UNPROTECT(1);
      error(
        "pareto sampling needs %d units but only %d have 0 < pik < 1",
        n_select, n_valid
      );
    }

    double *keys = (double *)R_alloc(n_valid, sizeof(double));
    int *orig = (int *)R_alloc(n_valid, sizeof(int));

    if (!has_prn)
      GetRNGstate();

    int vi = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > 0.0 && pk[i] < 1.0) {
        double u_i = has_prn ? u_ext[i] : unif_rand();
        double odds_u = u_i / (1.0 - u_i);
        double odds_p = pk[i] / (1.0 - pk[i]);
        keys[vi] = odds_u / odds_p;
        orig[vi] = i + 1; /* 1-based */
        vi++;
      }
    }

    if (!has_prn)
      PutRNGstate();

    /* Select the n_select smallest keys. When n_select == n_valid the
     * partition is a no-op; the loop below reads all n_valid entries and
     * the outer qsort(res, ...) restores ascending-index order. */
    partial_sort_paired(keys, orig, n_valid, n_select);

    for (int i = 0; i < n_select; i++) {
      res[j++] = orig[i];
    }
  }

  /* Sort full result (matches up_chromy.c idiom) */
  R_isort(res, total);

  UNPROTECT(1);
  return result;
}
