#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>

#define CHROMY_EPS 1e-9

/*
 * Core Chromy algorithm for fractional parts
 * pik values should be in [0,1), processes units in active[] order
 */
static int chromy_core(const double *pik, const int *active, int n_active,
                       int n_target, int *out) {
    if (n_target == 0 || n_active == 0) return 0;

    double V_prev = 0.0;
    int V_I_prev = 0;
    double V_F_prev = 0.0;
    int T = 0;
    int selected = 0;

    for (int i = 0; i < n_active; i++) {
        int k = active[i];
        double pk = pik[k];

        double V_curr = V_prev + pk;
        int V_I_curr = (int)floor(V_curr + CHROMY_EPS);
        double V_F_curr = V_curr - (double)V_I_curr;

        double prob;
        if (i == 0) {
            prob = pk;
        } else if (V_F_curr > V_F_prev + CHROMY_EPS) {
            prob = (T == V_I_prev) ? (V_F_curr - V_F_prev) / (1.0 - V_F_prev) : 0.0;
        } else {
            prob = (T == V_I_prev) ? 1.0 :
                   ((V_F_prev > CHROMY_EPS) ? V_F_curr / V_F_prev : 0.0);
        }

        if (unif_rand() < prob) {
            out[selected++] = k;
            T++;
        }

        V_prev = V_curr;
        V_I_prev = V_I_curr;
        V_F_prev = V_F_curr;
    }

    return selected;
}

static int select_start_weighted(const double *weights, const int *active,
                                 int n_active, double sum_weights) {
    if (n_active == 0) return 0;
    double u = unif_rand() * sum_weights;
    double cumsum = 0.0;
    for (int i = 0; i < n_active; i++) {
        cumsum += weights[active[i]];
        if (u <= cumsum) return i;
    }
    return n_active - 1;
}

static int chromy_precompute(
    const double *x, int N, int n, double sum_x,
    double *frac_pik,    /* out: fractional parts */
    int *floor_hits,     /* out: floor of pik */
    int *active,         /* out: indices with frac > eps */
    int *p_n_active,     /* out: count of active */
    double *p_sum_frac   /* out: sum of fractional parts */
) {
    int total_certain = 0;
    double sum_frac = 0.0;
    int n_active = 0;

    for (int k = 0; k < N; k++) {
        double pik_k = (double)n * x[k] / sum_x;
        floor_hits[k] = (int)pik_k;
        frac_pik[k] = pik_k - (double)floor_hits[k];
        total_certain += floor_hits[k];

        if (frac_pik[k] > CHROMY_EPS) {
            active[n_active++] = k;
            sum_frac += frac_pik[k];
        }
    }

    *p_n_active = n_active;
    *p_sum_frac = sum_frac;
    return total_certain;
}

static int chromy_sample_once(
    int N,
    const double *frac_pik,
    const int *floor_hits,
    const int *active,
    int n_active,
    double sum_frac,
    int n_frac,
    int *active_ord,
    int *frac_sel,
    int *sel
) {
    int nsel = 0;

    for (int k = 0; k < N; k++) {
        for (int j = 0; j < floor_hits[k]; j++) {
            sel[nsel++] = k;
        }
    }

    if (n_frac > 0 && n_active > 0) {
        int start_idx = select_start_weighted(frac_pik, active, n_active, sum_frac);

        for (int i = 0; i < n_active; i++) {
            active_ord[i] = active[(start_idx + i) % n_active];
        }

        int nf = chromy_core(frac_pik, active_ord, n_active, n_frac, frac_sel);

        for (int i = 0; i < nf; i++) {
            sel[nsel++] = frac_sel[i];
        }
    }

    return nsel;
}

SEXP C_up_chromy(SEXP r_x, SEXP r_n) {
    const int N = LENGTH(r_x);
    const double *x = REAL(r_x);
    const int n = asInteger(r_n);

    double sum_x = 0.0;
    for (int k = 0; k < N; k++) sum_x += x[k];

    double *frac_pik = (double *)R_alloc(N, sizeof(double));
    int *floor_hits = (int *)R_alloc(N, sizeof(int));
    int *active = (int *)R_alloc(N, sizeof(int));

    int n_active;
    double sum_frac;
    int total_certain = chromy_precompute(x, N, n, sum_x,
                                          frac_pik, floor_hits, active,
                                          &n_active, &sum_frac);
    int n_frac = n - total_certain;

    int *active_ord = (int *)R_alloc(n_active > 0 ? n_active : 1, sizeof(int));
    int *frac_sel = (int *)R_alloc(n_frac > 0 ? n_frac : 1, sizeof(int));
    int *sel = (int *)R_alloc(n, sizeof(int));

    GetRNGstate();
    int nsel = chromy_sample_once(N, frac_pik, floor_hits,
                                  active, n_active, sum_frac, n_frac,
                                  active_ord, frac_sel, sel);
    PutRNGstate();

    SEXP result = PROTECT(allocVector(INTSXP, nsel));
    int *out = INTEGER(result);
    for (int i = 0; i < nsel; i++) {
        out[i] = sel[i] + 1;
    }
    R_isort(out, nsel);

    UNPROTECT(1);
    return result;
}

SEXP C_up_chromy_pairexp(SEXP r_x, SEXP r_n, SEXP r_nsim) {
    const int N = LENGTH(r_x);
    const double *x = REAL(r_x);
    const int n = asInteger(r_n);
    const int nsim = asInteger(r_nsim);

    double sum_x = 0.0;
    for (int k = 0; k < N; k++) sum_x += x[k];

    double *frac_pik = (double *)R_alloc(N, sizeof(double));
    int *floor_hits = (int *)R_alloc(N, sizeof(int));
    int *active = (int *)R_alloc(N, sizeof(int));

    int n_active;
    double sum_frac;
    int total_certain = chromy_precompute(x, N, n, sum_x,
                                          frac_pik, floor_hits, active,
                                          &n_active, &sum_frac);
    int n_frac = n - total_certain;

    SEXP result = PROTECT(allocMatrix(REALSXP, N, N));
    double *joint = REAL(result);
    memset(joint, 0, (size_t)N * N * sizeof(double));

    int *active_ord = (int *)R_alloc(n_active > 0 ? n_active : 1, sizeof(int));
    int *frac_sel = (int *)R_alloc(n_frac > 0 ? n_frac : 1, sizeof(int));
    int *sel = (int *)R_alloc(n, sizeof(int));

    GetRNGstate();

    for (int sim = 0; sim < nsim; sim++) {
        int nsel = chromy_sample_once(N, frac_pik, floor_hits,
                                      active, n_active, sum_frac, n_frac,
                                      active_ord, frac_sel, sel);

        for (int i = 0; i < nsel; i++) {
            int si = sel[i];
            joint[si + si * N] += 1.0;
            for (int j = i + 1; j < nsel; j++) {
                int sj = sel[j];
                joint[si + sj * N] += 1.0;
                joint[sj + si * N] += 1.0;
            }
        }
    }

    PutRNGstate();

    double inv_nsim = 1.0 / (double)nsim;
    for (int i = 0; i < N * N; i++) {
        joint[i] *= inv_nsim;
    }

    UNPROTECT(1);
    return result;
}
