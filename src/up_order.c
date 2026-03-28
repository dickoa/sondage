/*
 * C implementations for systematic PPS, Poisson, SPS, and Pareto sampling.
 *
 * These replace the previous pure-R implementations for performance.
 * All functions handle certainty units (pik >= 1-eps) by selecting
 * them unconditionally and running the algorithm on valid units only.
 */

#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* helpers for order sampling (SPS / Pareto) */

/*
 * Partial sort on parallel arrays: rearranges keys[], orig[] so that the
 * n_select entries with the smallest keys are in positions [0, n_select).
 * Uses quickselect (Lomuto partition), O(N) expected.
 */
static void partial_sort_paired(double *keys, int *orig, int len,
                                int n_select) {
  if (n_select <= 0 || len <= 0 || n_select >= len)
    return;

  int lo = 0, hi = len - 1;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    /* Move pivot to end */
    {
      double tk = keys[mid];
      keys[mid] = keys[hi];
      keys[hi] = tk;
      int ti = orig[mid];
      orig[mid] = orig[hi];
      orig[hi] = ti;
    }

    double pivot = keys[hi];
    int store = lo;
    for (int j = lo; j < hi; j++) {
      if (keys[j] < pivot) {
        double tk = keys[store];
        keys[store] = keys[j];
        keys[j] = tk;
        int ti = orig[store];
        orig[store] = orig[j];
        orig[j] = ti;
        store++;
      }
    }
    {
      double tk = keys[store];
      keys[store] = keys[hi];
      keys[hi] = tk;
      int ti = orig[store];
      orig[store] = orig[hi];
      orig[hi] = ti;
    }

    if (store == n_select - 1)
      break;
    if (store < n_select - 1)
      lo = store + 1;
    else
      hi = store - 1;
  }
}

/* Integer comparison for qsort */
static int cmp_int(const void *a, const void *b) {
  int ia = *(const int *)a;
  int ib = *(const int *)b;
  return (ia > ib) - (ia < ib);
}

/* Systematic PPS */

/*
 * C_up_systematic(pik, eps)
 *
 * Single-pass systematic PPS: accumulate pik, select unit k when
 * the cumulative sum crosses an integer boundary (shifted by u).
 * Returns sorted 1-based indices.
 */
SEXP C_up_systematic(SEXP pik, SEXP eps) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  const double epsilon = REAL(eps)[0];

  /* Count certainty units and compute n */
  int n_certain = 0;
  double n_valid = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0 - epsilon)
      n_certain++;
    else if (pk[i] > epsilon)
      n_valid += pk[i];
  }
  int n_select = (int)(n_valid + 0.5);
  int total = n_select + n_certain;

  SEXP result = PROTECT(allocVector(INTSXP, total));
  int *res = INTEGER(result);
  int j = 0;

  /* Certainty units first */
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0 - epsilon)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    GetRNGstate();
    double u = unif_rand();
    PutRNGstate();

    double cs = 0.0;
    for (int i = 0; i < N; i++) {
      if (pk[i] <= epsilon || pk[i] >= 1.0 - epsilon)
        continue;
      double prev = cs;
      cs += pk[i];
      if (floor(cs - u) > floor(prev - u)) {
        res[j++] = i + 1;
      }
    }
  }

  /* Result is already sorted (sequential scan) */
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
 * C_up_sps(pik, prn, eps)
 *
 * Key: xi_k = u_k / pik_k.  Select the n smallest keys among valid units.
 * prn is R_NilValue when NULL (generate uniforms internally).
 * Returns sorted 1-based indices.
 */
SEXP C_up_sps(SEXP pik, SEXP prn, SEXP eps) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  const double epsilon = REAL(eps)[0];
  int has_prn = !isNull(prn);
  const double *u_ext = has_prn ? REAL(prn) : NULL;

  /* Separate certainty vs valid units */
  int n_certain = 0;
  double n_sum = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0 - epsilon)
      n_certain++;
    else if (pk[i] > epsilon)
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
    if (pk[i] >= 1.0 - epsilon)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    /* Build key array for valid units */
    int n_valid = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > epsilon && pk[i] < 1.0 - epsilon)
        n_valid++;
    }

    double *keys = (double *)R_alloc(n_valid, sizeof(double));
    int *orig = (int *)R_alloc(n_valid, sizeof(int));

    if (!has_prn)
      GetRNGstate();

    int vi = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > epsilon && pk[i] < 1.0 - epsilon) {
        double u_i = has_prn ? u_ext[i] : unif_rand();
        keys[vi] = u_i / pk[i];
        orig[vi] = i + 1; /* 1-based */
        vi++;
      }
    }

    if (!has_prn)
      PutRNGstate();

    /* Select the n_select smallest keys */
    partial_sort_paired(keys, orig, n_valid, n_select);

    for (int i = 0; i < n_select; i++) {
      res[j++] = orig[i];
    }
  }

  /* Sort full result */
  qsort(res, total, sizeof(int), cmp_int);

  UNPROTECT(1);
  return result;
}

/* Pareto sampling */

/*
 * C_up_pareto(pik, prn, eps)
 *
 * Key: xi_k = [u_k/(1-u_k)] / [pik_k/(1-pik_k)].
 * Select the n smallest keys among valid units.
 * prn is R_NilValue when NULL (generate uniforms internally).
 * Returns sorted 1-based indices.
 */
SEXP C_up_pareto(SEXP pik, SEXP prn, SEXP eps) {
  const int N = LENGTH(pik);
  const double *pk = REAL(pik);
  const double epsilon = REAL(eps)[0];
  int has_prn = !isNull(prn);
  const double *u_ext = has_prn ? REAL(prn) : NULL;

  /* Separate certainty vs valid units */
  int n_certain = 0;
  double n_sum = 0.0;
  for (int i = 0; i < N; i++) {
    if (pk[i] >= 1.0 - epsilon)
      n_certain++;
    else if (pk[i] > epsilon)
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
    if (pk[i] >= 1.0 - epsilon)
      res[j++] = i + 1;
  }

  if (n_select > 0) {
    int n_valid = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > epsilon && pk[i] < 1.0 - epsilon)
        n_valid++;
    }

    double *keys = (double *)R_alloc(n_valid, sizeof(double));
    int *orig = (int *)R_alloc(n_valid, sizeof(int));

    if (!has_prn)
      GetRNGstate();

    int vi = 0;
    for (int i = 0; i < N; i++) {
      if (pk[i] > epsilon && pk[i] < 1.0 - epsilon) {
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

    /* Select the n_select smallest keys */
    partial_sort_paired(keys, orig, n_valid, n_select);

    for (int i = 0; i < n_select; i++) {
      res[j++] = orig[i];
    }
  }

  /* Sort full result */
  qsort(res, total, sizeof(int), cmp_int);

  UNPROTECT(1);
  return result;
}
