/*
 * Sampford unequal-probability sampling without replacement.
 *
 * The fast path is Sampford's original rejection construction.  To avoid
 * unbounded running time, difficult inputs fall back to Grafstrom's exact
 * non-rejective construction using the conditional-Poisson routines shared
 * with the CPS implementation.
 *
 * References
 *   Sampford, M.R. (1967), Biometrika 54, 499--513.
 *   Grafstrom, A. (2009), J. Statist. Plann. Inference 139, 2111--2114.
 */

#include "cps_core.h"
#include <stdlib.h>

#define SAMPFORD_REJECTION_TRIES 32

static int int_cmp(const void *a, const void *b) {
    const int ia = *(const int *)a;
    const int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}

/* Draw from positive weights represented by their cumulative sums. */
static int sampford_categorical(const double *cum, int N, double total) {
    const double u = unif_rand() * total;
    int lo = 0, hi = N - 1;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (u < cum[mid]) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

static int sampford_has_duplicate(const int *x, int len, int value) {
    for (int i = 0; i < len; i++) {
        if (x[i] == value) return 1;
    }
    return 0;
}

/* Partial Fisher--Yates draw for the equal-probability fast path. */
static void sampford_srs(int N, int m, int *out) {
    int *pool = (int *) R_alloc(N, sizeof(int));
    for (int i = 0; i < N; i++) pool[i] = i;
    for (int i = 0; i < m; i++) {
        int j = i + (int) floor(unif_rand() * (double)(N - i));
        if (j >= N) j = N - 1; /* guard a hypothetical rounded U == 1 */
        const int tmp = pool[i];
        pool[i] = pool[j];
        pool[j] = tmp;
        out[i] = pool[i];
    }
}

/*
 * Grafstrom's non-rejective representation.
 *
 * alpha_i is the inclusion probability of i under conditional-Poisson
 * sampling of size m-1 with Bernoulli parameters p.  Select the anchor with
 * weight p_i (1-alpha_i), then select m-1 other units by conditional Poisson
 * with odds p_i/(1-p_i).  This is exactly the Sampford design.
 */
static void sampford_nonrejective(const double *p, const double *odds,
                                  int N, int m, int *out) {
    if (m == 1) {
        double *cum = (double *) R_alloc(N, sizeof(double));
        double total = 0.0;
        for (int i = 0; i < N; i++) {
            total += p[i];
            cum[i] = total;
        }
        out[0] = sampford_categorical(cum, N, total);
        return;
    }

    const int r = m - 1;
    double *f = (double *) R_alloc(CPS_TABLE_SIZE(N, r), sizeof(double));
    double *alpha = (double *) R_alloc(N, sizeof(double));
    double *pro = (double *) R_alloc((size_t)N * r, sizeof(double));
    double *anchor_cum = (double *) R_alloc(N, sizeof(double));

    cps_compute_f(odds, N, r, f);
    cps_compute_pik(odds, f, N, r, alpha, pro);

    double anchor_total = 0.0;
    for (int i = 0; i < N; i++) {
        double a = alpha[i];
        if (a < 0.0) a = 0.0;
        if (a > 1.0) a = 1.0;
        anchor_total += p[i] * (1.0 - a);
        anchor_cum[i] = anchor_total;
    }
    if (!(anchor_total > 0.0) || !R_FINITE(anchor_total)) {
        error("Sampford non-rejective anchor probabilities are non-finite");
    }

    const int anchor = sampford_categorical(anchor_cum, N, anchor_total);
    out[0] = anchor;

    double *odds_ex = (double *) R_alloc(N - 1, sizeof(double));
    int *map = (int *) R_alloc(N - 1, sizeof(int));
    int pos = 0;
    for (int i = 0; i < N; i++) {
        if (i == anchor) continue;
        odds_ex[pos] = odds[i];
        map[pos] = i;
        pos++;
    }

    double *f_ex = (double *) R_alloc(CPS_TABLE_SIZE(N - 1, r),
                                      sizeof(double));
    int *draw_ex = (int *) R_alloc(r, sizeof(int));
    cps_compute_f(odds_ex, N - 1, r, f_ex);
    const int selected = cps_sample(odds_ex, f_ex, N - 1, r, draw_ex);
    if (selected != r) {
        error("Sampford conditional-Poisson fallback produced %d units, expected %d",
              selected, r);
    }
    for (int i = 0; i < r; i++) out[i + 1] = map[draw_ex[i]];
}

/* Draw a work-design sample, using rejection first and exact fallback. */
static void sampford_draw_work(const double *p, int N, int m, int *out) {
    if (m <= 0) return;

    int all_equal = 1;
    const double p0 = p[0];
    const double tol = 32.0 * DBL_EPSILON * fmax(1.0, fabs(p0));
    for (int i = 1; i < N; i++) {
        if (fabs(p[i] - p0) > tol) {
            all_equal = 0;
            break;
        }
    }
    if (all_equal) {
        sampford_srs(N, m, out);
        return;
    }

    double *odds = (double *) R_alloc(N, sizeof(double));
    double *cum_p = (double *) R_alloc(N, sizeof(double));
    double *cum_odds = (double *) R_alloc(N, sizeof(double));
    double p_total = 0.0, odds_total = 0.0;
    for (int i = 0; i < N; i++) {
        odds[i] = p[i] / (1.0 - p[i]);
        p_total += p[i];
        odds_total += odds[i];
        cum_p[i] = p_total;
        cum_odds[i] = odds_total;
    }
    if (!(odds_total > 0.0) || !R_FINITE(odds_total)) {
        error("Sampford odds are non-finite; inclusion probabilities are too close to 1");
    }

    if (m == 1) {
        out[0] = sampford_categorical(cum_p, N, p_total);
        return;
    }

    for (int attempt = 0; attempt < SAMPFORD_REJECTION_TRIES; attempt++) {
        out[0] = sampford_categorical(cum_p, N, p_total);
        int accepted = 1;
        for (int j = 1; j < m; j++) {
            const int unit = sampford_categorical(cum_odds, N, odds_total);
            if (sampford_has_duplicate(out, j, unit)) {
                accepted = 0;
                break;
            }
            out[j] = unit;
        }
        if (accepted) return;
    }

    sampford_nonrejective(p, odds, N, m, out);
}

SEXP C_up_sampford(SEXP pik_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik = REAL(pik_sexp);

    int N = 0, n_certain = 0;
    double n_active_sum = 0.0;
    for (int i = 0; i < N_full; i++) {
        if (pik[i] >= 1.0) n_certain++;
        else if (pik[i] > 0.0) {
            N++;
            n_active_sum += pik[i];
        }
    }
    const int n_active = (int) floor(n_active_sum + 0.5);
    const int n_total = n_certain + n_active;

    SEXP ans = PROTECT(allocVector(INTSXP, n_total));
    int *result = INTEGER(ans);
    int *active_map = (int *) R_alloc((N > 0 ? N : 1), sizeof(int));
    double *p_orig = (double *) R_alloc((N > 0 ? N : 1), sizeof(double));

    int ai = 0;
    for (int i = 0; i < N_full; i++) {
        if (pik[i] > 0.0 && pik[i] < 1.0) {
            active_map[ai] = i;
            p_orig[ai] = pik[i];
            ai++;
        }
    }

    if (N == 0 || n_active == 0) {
        int k = 0;
        for (int i = 0; i < N_full; i++) {
            if (pik[i] >= 1.0) result[k++] = i + 1;
        }
        UNPROTECT(1);
        return ans;
    }

    int use_complement = (n_active > N - n_active);
    if (use_complement) {
        /* 1-p can round to exactly 1 for subnormal p; avoid infinite odds. */
        for (int i = 0; i < N; i++) {
            if (!(1.0 - p_orig[i] < 1.0)) {
                use_complement = 0;
                break;
            }
        }
    }
    const int m = use_complement ? N - n_active : n_active;
    double *p_work = (double *) R_alloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        p_work[i] = use_complement ? 1.0 - p_orig[i] : p_orig[i];
    }
    int *work_sample = (int *) R_alloc((m > 0 ? m : 1), sizeof(int));

    GetRNGstate();
    sampford_draw_work(p_work, N, m, work_sample);
    PutRNGstate();

    unsigned char *marked = (unsigned char *) R_alloc(N, sizeof(unsigned char));
    memset(marked, 0, (size_t)N * sizeof(unsigned char));
    for (int i = 0; i < m; i++) marked[work_sample[i]] = 1;

    int out = 0;
    for (int i = 0; i < N_full; i++) {
        if (pik[i] >= 1.0) result[out++] = i + 1;
    }
    for (int i = 0; i < N; i++) {
        const int take = use_complement ? !marked[i] : marked[i];
        if (take) result[out++] = active_map[i] + 1;
    }
    if (out != n_total) {
        UNPROTECT(1);
        error("Sampford sampling produced %d units, expected %d", out, n_total);
    }
    qsort(result, (size_t)n_total, sizeof(int), int_cmp);

    UNPROTECT(1);
    return ans;
}

/* ------------------------------------------------------------------------- */
/* Exact second-order inclusion probabilities                                */

/* Positive-term ESP recurrence, excluding a pair. Used as a stable fallback. */
static double sampford_pair_direct(const double *p, const long double *q,
                                   int N, int m, int ii, int jj,
                                   long double normalizer,
                                   long double *e, long double *g) {
    const int d = m - 2;
    memset(e, 0, (size_t)(d + 1) * sizeof(long double));
    memset(g, 0, (size_t)(d + 1) * sizeof(long double));
    e[0] = 1.0L;
    int seen = 0;
    for (int k = 0; k < N; k++) {
        if (k == ii || k == jj) continue;
        int top = seen + 1;
        if (top > d) top = d;
        for (int r = top; r >= 1; r--) {
            g[r] += q[k] * (g[r - 1] + (long double)p[k] * e[r - 1]);
            e[r] += q[k] * e[r - 1];
        }
        seen++;
    }
    const long double raw = q[ii] * q[jj] *
        (((long double)m - p[ii] - p[jj]) * e[d] - g[d]);
    return (double)(raw / normalizer);
}

/* Exact pair value from a quadratic deletion of units i,j. */
static double sampford_pair_fast(const double *p, const long double *q,
                                 const long double *E, const long double *G,
                                 int N, int m, int ii, int jj,
                                 long double normalizer,
                                 long double *e, long double *g) {
    const int d = m - 2;
    e[0] = 1.0L;
    g[0] = 0.0L;
    const long double qi = q[ii], qj = q[jj];
    const long double s1 = qi + qj;
    const long double s2 = qi * qj;
    const long double h1 = (long double)p[ii] * qi +
                           (long double)p[jj] * qj;
    const long double h2 = s2 * ((long double)p[ii] + p[jj]);

    int bad = 0;
    for (int r = 1; r <= d; r++) {
        const long double em1 = e[r - 1];
        const long double em2 = (r >= 2) ? e[r - 2] : 0.0L;
        const long double gm1 = g[r - 1];
        const long double gm2 = (r >= 2) ? g[r - 2] : 0.0L;
        e[r] = E[r] - s1 * em1 - s2 * em2;
        g[r] = G[r] - s1 * gm1 - s2 * gm2 - h1 * em1 - h2 * em2;
        if (!isfinite(e[r]) || !isfinite(g[r]) ||
            e[r] < -64.0L * LDBL_EPSILON * fabsl(E[r])) {
            bad = 1;
            break;
        }
        if (e[r] < 0.0L) e[r] = 0.0L;
        if (g[r] < 0.0L &&
            g[r] > -64.0L * LDBL_EPSILON * fmaxl(1.0L, fabsl(G[r]))) {
            g[r] = 0.0L;
        }
    }

    double value = NAN;
    if (!bad) {
        const long double raw = qi * qj *
            (((long double)m - p[ii] - p[jj]) * e[d] - g[d]);
        value = (double)(raw / normalizer);
        const double lo = fmax(0.0, p[ii] + p[jj] - 1.0);
        const double hi = fmin(p[ii], p[jj]);
        const double tol = 2e-10 * fmax(1.0, hi);
        if (!R_FINITE(value) || value < lo - tol || value > hi + tol) bad = 1;
    }
    if (bad) {
        value = sampford_pair_direct(p, q, N, m, ii, jj, normalizer, e, g);
    }

    const double lo = fmax(0.0, p[ii] + p[jj] - 1.0);
    const double hi = fmin(p[ii], p[jj]);
    if (value < lo && value > lo - 2e-12) value = lo;
    if (value > hi && value < hi + 2e-12) value = hi;
    return value;
}

/*
 * Common implementation for the full matrix and a requested submatrix.
 * idx is NULL for the full matrix and otherwise contains 1-based indices.
 */
static SEXP sampford_jip_impl(SEXP pik_sexp, const int *idx, int n_out) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik = REAL(pik_sexp);
    if (idx == NULL) n_out = N_full;

    SEXP ans = PROTECT(allocMatrix(REALSXP, n_out, n_out));
    double *J = REAL(ans);
    memset(J, 0, (size_t)n_out * n_out * sizeof(double));

    int N = 0;
    double n_active_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik[k] > 0.0 && pik[k] < 1.0) {
            N++;
            n_active_sum += pik[k];
        }
    }
    const int n_active = (int) floor(n_active_sum + 0.5);

    int *pop_to_active = (int *) R_alloc(N_full, sizeof(int));
    for (int k = 0; k < N_full; k++) pop_to_active[k] = -1;
    double *p_orig = (double *) R_alloc((N > 0 ? N : 1), sizeof(double));
    int ai = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik[k] > 0.0 && pik[k] < 1.0) {
            pop_to_active[k] = ai;
            p_orig[ai++] = pik[k];
        }
    }

    /* Diagonal, certainty pairs, and zero-probability pairs. */
    for (int a = 0; a < n_out; a++) {
        const int ka = idx ? idx[a] - 1 : a;
        J[(size_t)a * n_out + a] = pik[ka];
        for (int b = 0; b < a; b++) {
            const int kb = idx ? idx[b] - 1 : b;
            double v = NAN;
            if (pik[ka] <= 0.0 || pik[kb] <= 0.0) v = 0.0;
            else if (pik[ka] >= 1.0) v = pik[kb];
            else if (pik[kb] >= 1.0) v = pik[ka];
            if (!ISNAN(v)) {
                J[(size_t)a * n_out + b] = v;
                J[(size_t)b * n_out + a] = v;
            }
        }
    }

    if (N < 2 || n_active < 1) {
        UNPROTECT(1);
        return ans;
    }

    int use_complement = (n_active > N - n_active);
    if (use_complement) {
        for (int i = 0; i < N; i++) {
            if (!(1.0 - p_orig[i] < 1.0)) {
                use_complement = 0;
                break;
            }
        }
    }
    const int m = use_complement ? N - n_active : n_active;
    double *p = (double *) R_alloc(N, sizeof(double));
    for (int i = 0; i < N; i++) p[i] = use_complement ? 1.0 - p_orig[i] : p_orig[i];

    int all_equal = 1;
    const double eqtol = 32.0 * DBL_EPSILON * fmax(1.0, fabs(p[0]));
    for (int i = 1; i < N; i++) {
        if (fabs(p[i] - p[0]) > eqtol) {
            all_equal = 0;
            break;
        }
    }

    if (all_equal) {
        const double v = (N > 1) ?
            (double)n_active * (n_active - 1) / ((double)N * (N - 1)) : 0.0;
        for (int a = 0; a < n_out; a++) {
            const int ka = idx ? idx[a] - 1 : a;
            if (pop_to_active[ka] < 0) continue;
            for (int b = 0; b < a; b++) {
                const int kb = idx ? idx[b] - 1 : b;
                if (pop_to_active[kb] < 0) continue;
                J[(size_t)a * n_out + b] = v;
                J[(size_t)b * n_out + a] = v;
            }
        }
        UNPROTECT(1);
        return ans;
    }

    if (m < 2) {
        /* Work-design pairs never co-occur. */
        for (int a = 0; a < n_out; a++) {
            const int ka = idx ? idx[a] - 1 : a;
            const int iact = pop_to_active[ka];
            if (iact < 0) continue;
            for (int b = 0; b < a; b++) {
                const int kb = idx ? idx[b] - 1 : b;
                const int jact = pop_to_active[kb];
                if (jact < 0) continue;
                double v = use_complement ? p_orig[iact] + p_orig[jact] - 1.0 : 0.0;
                if (v < 0.0 && v > -2e-12) v = 0.0;
                J[(size_t)a * n_out + b] = v;
                J[(size_t)b * n_out + a] = v;
            }
        }
        UNPROTECT(1);
        return ans;
    }

    /* Scale all odds by one common constant so sum(q) = m.  Sampford
     * probabilities are homogeneous of degree m, so this changes neither
     * the design nor any joint probability and greatly limits ESP growth. */
    double max_odds = 0.0;
    for (int i = 0; i < N; i++) {
        const double o = p[i] / (1.0 - p[i]);
        if (o > max_odds) max_odds = o;
    }
    long double scaled_sum = 0.0L;
    for (int i = 0; i < N; i++) {
        scaled_sum += (long double)(p[i] / (1.0 - p[i])) / max_odds;
    }
    long double *q = (long double *) R_alloc(N, sizeof(long double));
    for (int i = 0; i < N; i++) {
        q[i] = ((long double)(p[i] / (1.0 - p[i])) / max_odds) *
               (long double)m / scaled_sum;
    }

    long double *E = (long double *) R_alloc(m + 1, sizeof(long double));
    long double *G = (long double *) R_alloc(m + 1, sizeof(long double));
    memset(E, 0, (size_t)(m + 1) * sizeof(long double));
    memset(G, 0, (size_t)(m + 1) * sizeof(long double));
    E[0] = 1.0L;
    int seen = 0;
    for (int i = 0; i < N; i++) {
        int top = seen + 1;
        if (top > m) top = m;
        for (int r = top; r >= 1; r--) {
            G[r] += q[i] * (G[r - 1] + (long double)p[i] * E[r - 1]);
            E[r] += q[i] * E[r - 1];
        }
        seen++;
    }
    const long double normalizer = (long double)m * E[m] - G[m];
    if (!(normalizer > 0.0L) || !isfinite(normalizer)) {
        UNPROTECT(1);
        error("Sampford joint-probability normalizer is non-finite");
    }

    const int d = m - 2;
    long double *e = (long double *) R_alloc(d + 1, sizeof(long double));
    long double *g = (long double *) R_alloc(d + 1, sizeof(long double));

    for (int a = 0; a < n_out; a++) {
        if ((a & 127) == 0) R_CheckUserInterrupt();
        const int ka = idx ? idx[a] - 1 : a;
        const int iact = pop_to_active[ka];
        if (iact < 0) continue;
        for (int b = 0; b < a; b++) {
            const int kb = idx ? idx[b] - 1 : b;
            const int jact = pop_to_active[kb];
            if (jact < 0) continue;
            double work_v = sampford_pair_fast(p, q, E, G, N, m,
                                                iact, jact, normalizer, e, g);
            double v = work_v;
            if (use_complement) {
                v = p_orig[iact] + p_orig[jact] - 1.0 + work_v;
            }
            const double lo = fmax(0.0, p_orig[iact] + p_orig[jact] - 1.0);
            const double hi = fmin(p_orig[iact], p_orig[jact]);
            if (v < lo && v > lo - 2e-10) v = lo;
            if (v > hi && v < hi + 2e-10) v = hi;
            J[(size_t)a * n_out + b] = v;
            J[(size_t)b * n_out + a] = v;
        }
    }

    UNPROTECT(1);
    return ans;
}

SEXP C_sampford_jip(SEXP pik_sexp) {
    return sampford_jip_impl(pik_sexp, NULL, LENGTH(pik_sexp));
}

SEXP C_sampford_jip_sub(SEXP pik_sexp, SEXP idx_sexp) {
    return sampford_jip_impl(pik_sexp, INTEGER(idx_sexp), LENGTH(idx_sexp));
}
