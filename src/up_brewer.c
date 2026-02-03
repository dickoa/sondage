#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
 * Brewer's method for unequal probability sampling
 *
 * Optimizations:
 * - Binary search for selection (O(log N) vs O(N))
 * - Compact arrays with swap-remove
 */

static inline int binary_search(const double *cumprob, int n, double u) {
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (cumprob[mid] <= u) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

SEXP C_up_brewer(SEXP pik, SEXP eps) {
    const int N_full = LENGTH(pik);
    const double *pik_ptr = REAL(pik);
    const double epsilon = REAL(eps)[0];

    int N = 0;
    int N_certain = 0;
    double n = 0.0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            N_certain++;
        } else if (pk > epsilon) {
            N++;
            n += pk;
        }
    }

    const int n_draws = (int)(n + 0.5);
    const int total_selected = n_draws + N_certain;

    SEXP s = PROTECT(allocVector(INTSXP, total_selected));
    int *s_ptr = INTEGER(s);
    int s_idx = 0;

    if (N == 0 || n_draws == 0) {
        for (int k = 0; k < N_full; k++) {
            if (pik_ptr[k] >= 1.0 - epsilon) {
                s_ptr[s_idx++] = k + 1;
            }
        }
        UNPROTECT(1);
        return s;
    }

    double *pikb = (double *) R_alloc(N, sizeof(double));
    double *cumprob = (double *) R_alloc(N, sizeof(double));
    int *orig_idx = (int *) R_alloc(N, sizeof(int));

    int j = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            s_ptr[s_idx++] = k + 1;
        } else if (pk > epsilon) {
            pikb[j] = pk;
            orig_idx[j] = k;
            j++;
        }
    }

    GetRNGstate();

    double a = 0.0;
    int n_active = N;

    for (int i = 1; i <= n_draws; i++) {
        const double m = n - a;
        const double r = (double)(n_draws - i + 1);

        double sum_p = 0.0;
        for (int j = 0; j < n_active; j++) {
            const double pk = pikb[j];
            const double denom = m - pk * r;

            if (denom > 1e-10) {
                double prob = pk * (m - pk) / denom;
                if (prob > 0.0) sum_p += prob;
            }
            cumprob[j] = sum_p;
        }

        const double u = unif_rand() * sum_p;
        int sel_j = binary_search(cumprob, n_active, u);

        s_ptr[s_idx++] = orig_idx[sel_j] + 1;
        a += pikb[sel_j];

        n_active--;
        if (sel_j < n_active) {
            pikb[sel_j] = pikb[n_active];
            orig_idx[sel_j] = orig_idx[n_active];
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return s;
}
