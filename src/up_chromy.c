#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define CHROMY_EPS 1e-9

/*
 * Core Chromy algorithm for fractional parts
 * pik values should be in [0,1), units with pik=0 are skipped
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

SEXP C_up_chromy(SEXP r_x, SEXP r_n) {
    const int N = LENGTH(r_x);
    const double *x = REAL(r_x);
    const int n = asInteger(r_n);

    double sum_x = 0.0;
    for (int k = 0; k < N; k++) {
        sum_x += x[k];
    }

    double *pik = (double *)R_alloc(N, sizeof(double));
    int *floor_hits = (int *)R_alloc(N, sizeof(int));
    double *frac_pik = (double *)R_alloc(N, sizeof(double));
    int *active = (int *)R_alloc(N, sizeof(int));

    int total_certain = 0;
    double sum_frac = 0.0;
    int n_active = 0;

    for (int k = 0; k < N; k++) {
        pik[k] = (double)n * x[k] / sum_x;
        floor_hits[k] = (int)pik[k];  /* truncate, don't add eps */
        frac_pik[k] = pik[k] - (double)floor_hits[k];
        total_certain += floor_hits[k];

        if (frac_pik[k] > CHROMY_EPS) {
            active[n_active++] = k;
            sum_frac += frac_pik[k];
        }
    }

    int n_frac = n - total_certain;

    GetRNGstate();

    int start_idx = 0;
    if (n_frac > 0 && n_active > 0) {
        start_idx = select_start_weighted(frac_pik, active, n_active, sum_frac);
    }

    int *active_ord = (int *)R_alloc(n_active, sizeof(int));
    for (int i = 0; i < n_active; i++) {
        active_ord[i] = active[(start_idx + i) % n_active];
    }

    int *frac_sel = (int *)R_alloc(n_frac > 0 ? n_frac : 1, sizeof(int));
    int n_frac_actual = chromy_core(frac_pik, active_ord, n_active, n_frac, frac_sel);

    PutRNGstate();

    SEXP result = PROTECT(allocVector(INTSXP, total_certain + n_frac_actual));
    int *out = INTEGER(result);
    int out_idx = 0;

    for (int k = 0; k < N; k++) {
        for (int j = 0; j < floor_hits[k]; j++) {
            out[out_idx++] = k + 1;
        }
    }

    for (int i = 0; i < n_frac_actual; i++) {
        out[out_idx++] = frac_sel[i] + 1;
    }

    R_isort(out, total_certain + n_frac_actual);
    UNPROTECT(1);
    return result;
}
