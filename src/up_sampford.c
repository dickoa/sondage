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
