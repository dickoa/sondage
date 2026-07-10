/*
 * Joint Inclusion Probabilities
 *
 * - C_cps_jip: Exact CPS joint probabilities (Aires' formula)
 * - C_up_systematic_jip: Exact systematic sampling joint probabilities
 * - C_high_entropy_jip: High entropy approximation (equation 18, Brewer & Donadio 2003)
 */

#include "cps_core.h"
#include <float.h>
#include <R_ext/Utils.h>

/*
 * Exact CPS joint probability for one pair of valid units (i, j),
 * computed in the probability domain:
 *
 *   pi_ij = p_i * p_j * P(W \ {i,j} selects n-2) / P(W selects n)
 *
 * with p_k = w_k / (1 + w_k). The numerator is a fresh forward
 * Poisson-binomial recurrence excluding units i and j: O(N * n) time,
 * O(n) workspace, and no subtraction of near-equal quantities. (The
 * previous leave-one/two-out polynomial division cancelled
 * catastrophically for repeated weights and produced NaN at N ~ 800.)
 *
 * buf must hold n + 1 doubles. f_n = P(W selects n) from the
 * calibration table; it is never small because calibration centers the
 * unconditional Poisson-binomial distribution on n.
 */
static double cps_pair_exact(const double *w, int N, int n,
                             int i, int j, double f_n, double *buf) {
    if (n < 2 || !(f_n > 0.0)) return 0.0;
    const int m = n - 2;

    memset(buf, 0, (size_t)(m + 1) * sizeof(double));
    buf[0] = 1.0;
    int zmax = 0;

    for (int k = 0; k < N; k++) {
        if (k == i || k == j) continue;
        const double p = w[k] / (1.0 + w[k]);
        const int zt = (zmax < m) ? zmax + 1 : m;
        for (int z = zt; z >= 1; z--) {
            buf[z] = (1.0 - p) * buf[z] + p * buf[z - 1];
        }
        buf[0] *= (1.0 - p);
        zmax = zt;
    }

    const double p_i = w[i] / (1.0 + w[i]);
    const double p_j = w[j] / (1.0 + w[j]);
    return p_i * p_j * buf[m] / f_n;
}

/*
 * Group boundaries for weights sorted ascending: a new group starts
 * whenever the gap to the previous weight exceeds a few ulps. Units
 * with equal calibrated weights are exchangeable, so every pair drawn
 * from a fixed pair of groups shares one joint probability.
 *
 * Fills gstart[0..G] (G+1 entries, gstart[G] = len) and returns G.
 */
static int cps_weight_groups(const double *w_sorted, int len, int *gstart) {
    const double tie = 64.0 * DBL_EPSILON;
    int G = 1;
    gstart[0] = 0;
    for (int t = 1; t < len; t++) {
        if (w_sorted[t] - w_sorted[t - 1] > tie * fmax(1.0, w_sorted[t])) {
            gstart[G++] = t;
        }
    }
    gstart[G] = len;
    return G;
}

/*
 * Transform a working-design pair probability back to the original
 * design and clamp it into its Frechet bounds. When the calibration ran
 * on the complement design (n > N - n), v is the probability that both
 * units are excluded, and pi_ij = pik_i + pik_j - 1 + v.
 */
static double cps_pair_finalize(double v, double pik_i, double pik_j,
                                int use_complement) {
    double pi_ij = use_complement ? (pik_i + pik_j - 1.0 + v) : v;
    double lo = pik_i + pik_j - 1.0;
    if (lo < 0.0) lo = 0.0;
    if (pi_ij < lo) pi_ij = lo;
    if (pi_ij > pik_i) pi_ij = pik_i;
    if (pi_ij > pik_j) pi_ij = pik_j;
    return pi_ij;
}

SEXP C_cps_jip(SEXP pik_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);

    double n_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        n_sum += pik_full[k];
    }
    const int n = (int)(n_sum + 0.5);

    SEXP result = PROTECT(allocMatrix(REALSXP, N_full, N_full));
    double *pikl = REAL(result);
    memset(pikl, 0, (size_t)N_full * N_full * sizeof(double));

    for (int k = 0; k < N_full; k++) {
        pikl[(size_t)k * N_full + k] = pik_full[k];
    }

    if (n < 2) {
        UNPROTECT(1);
        return result;
    }

    /* Certainty units: pi_kj = pik_j for every j */
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] >= 1.0) {
            for (int l = 0; l < N_full; l++) {
                pikl[(size_t)k * N_full + l] = pik_full[l];
                pikl[(size_t)l * N_full + k] = pik_full[l];
            }
            pikl[(size_t)k * N_full + k] = pik_full[k];
        }
    }

    int N_valid = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            N_valid++;
        }
    }

    if (N_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    int *valid_idx = (int *) R_alloc(N_valid, sizeof(int));
    double *pik_valid = (double *) R_alloc(N_valid, sizeof(double));

    int j = 0;
    double n_valid_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            valid_idx[j] = k;
            pik_valid[j] = pik_full[k];
            n_valid_sum += pik_full[k];
            j++;
        }
    }
    const int n_valid = (int)(n_valid_sum + 0.5);

    if (n_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    /* Calibrate in the smaller of the design and its complement: the
     * complement of a maximum-entropy design is maximum entropy with
     * pik_c = 1 - pik, and q_ij (both excluded) maps back through
     * pi_ij = pik_i + pik_j - 1 + q_ij. This bounds the table at
     * O(N * min(n, N - n)). */
    const int use_complement = (n_valid > N_valid - n_valid);
    const int n_work = use_complement ? (N_valid - n_valid) : n_valid;

    double *pik_work = (double *) R_alloc(N_valid, sizeof(double));
    for (int i = 0; i < N_valid; i++) {
        pik_work[i] = use_complement ? (1.0 - pik_valid[i]) : pik_valid[i];
    }

    double *w = (double *) R_alloc(N_valid, sizeof(double));
    double *f = (double *) R_alloc(CPS_TABLE_SIZE(N_valid, n_work),
                                   sizeof(double));
    double cal_max_diff = 0.0;
    int    cal_worst_idx = -1;
    int cal_iters = cps_calibrate(pik_work, N_valid, n_work, w, f, 1e-9,
                                  500, &cal_max_diff, &cal_worst_idx);
    if (cal_iters >= 500) {
        cps_warn_nonconverge("joint_inclusion_prob: ", 500, 1e-9,
                             cal_max_diff, pik_valid, N_valid, valid_idx);
    }
    (void) cal_worst_idx;

    const double f_n = f[n_work]; /* row 0, column n_work */

    /* Relative near-equality threshold for the Aires difference formula */
    const double tau = 16.0 * sqrt(DBL_EPSILON);

    /* Check if all weights are approximately equal */
    double w_min = w[0], w_max = w[0];
    for (int i = 1; i < N_valid; i++) {
        if (w[i] < w_min) w_min = w[i];
        if (w[i] > w_max) w_max = w[i];
    }
    double w_scale = fmax(1.0, w_max);
    int all_equal = (w_max - w_min) <= tau * w_scale;

    if (all_equal) {
        /* All weights equal => CPS reduces to SRS on valid subset */
        double pi_ij = (double)n_valid * (n_valid - 1) /
                       ((double)N_valid * (N_valid - 1));
        for (int i = 0; i < N_valid; i++) {
            for (int jj = i + 1; jj < N_valid; jj++) {
                int k = valid_idx[i];
                int l = valid_idx[jj];
                pikl[(size_t)l * N_full + k] = pi_ij;
                pikl[(size_t)k * N_full + l] = pi_ij;
            }
        }
        UNPROTECT(1);
        return result;
    }

    /* Sort valid units by calibrated weight and group exact ties. */
    double *w_sorted = (double *) R_alloc(N_valid, sizeof(double));
    int *ord = (int *) R_alloc(N_valid, sizeof(int));
    memcpy(w_sorted, w, N_valid * sizeof(double));
    for (int i = 0; i < N_valid; i++) ord[i] = i;
    rsort_with_index(w_sorted, ord, N_valid);

    int *gstart = (int *) R_alloc(N_valid + 1, sizeof(int));
    const int G = cps_weight_groups(w_sorted, N_valid, gstart);

    double *buf = (double *) R_alloc(n_work + 1, sizeof(double));

    for (int g = 0; g < G; g++) {
        R_CheckUserInterrupt();
        for (int h = g; h < G; h++) {
            const int gs = gstart[g], ge = gstart[g + 1];
            const int hs = gstart[h], he = gstart[h + 1];
            if (g == h && ge - gs < 2) continue;

            const double w_lo = w_sorted[gs];
            const double w_hi = w_sorted[hs];
            const double scale = fmax(1.0, w_hi);
            const int near = (w_hi - w_lo) <= tau * scale;

            double v_block = 0.0;
            if (near) {
                /* One exact evaluation covers every pair in the block:
                 * units within a group are exchangeable. */
                const int rep_i = ord[gs];
                const int rep_j = (g == h) ? ord[gs + 1] : ord[hs];
                v_block = cps_pair_exact(w, N_valid, n_work,
                                         rep_i, rep_j, f_n, buf);
            }

            for (int a = gs; a < ge; a++) {
                const int b0 = (g == h) ? a + 1 : hs;
                for (int b = b0; b < he; b++) {
                    const int vi = ord[a];
                    const int vj = ord[b];
                    double v;
                    if (near) {
                        v = v_block;
                    } else {
                        /* Aires' identity; the block gap guarantees a
                         * well-separated denominator. */
                        const double r_i = w[vi];
                        const double r_j = w[vj];
                        v = (r_j * pik_work[vi] - r_i * pik_work[vj]) /
                            (r_j - r_i);
                    }
                    const double pi_ij = cps_pair_finalize(
                        v, pik_valid[vi], pik_valid[vj], use_complement);

                    const int k = valid_idx[vi];
                    const int l = valid_idx[vj];
                    pikl[(size_t)l * N_full + k] = pi_ij;
                    pikl[(size_t)k * N_full + l] = pi_ij;
                }
            }
        }
    }

    UNPROTECT(1);
    return result;
}

/*
 * CPS JIP submatrix for sampled units only.
 * idx_sexp: 1-based integer vector of sampled population indices.
 * Returns n_s x n_s matrix instead of N x N.
 *
 * Calibration (the dominant cost) still runs on the full valid
 * population, but pair values are only computed for sampled units,
 * using the same group-block structure as the full matrix.
 */
SEXP C_cps_jip_sub(SEXP pik_sexp, SEXP idx_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);
    const int n_s = LENGTH(idx_sexp);
    const int *idx_r = INTEGER(idx_sexp); /* 1-based R indices */

    double n_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        n_sum += pik_full[k];
    }
    const int n = (int)(n_sum + 0.5);

    /* Allocate n_s x n_s output */
    SEXP result = PROTECT(allocMatrix(REALSXP, n_s, n_s));
    double *pikl = REAL(result);
    memset(pikl, 0, (size_t)n_s * n_s * sizeof(double));

    for (int i = 0; i < n_s; i++) {
        pikl[(size_t)i * n_s + i] = pik_full[idx_r[i] - 1];
    }

    if (n < 2) {
        UNPROTECT(1);
        return result;
    }

    int N_valid = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            N_valid++;
        }
    }

    /* Handle certainty units in the output:
     * certainty-any pair => pi_ij = pik[other] */
    for (int i = 0; i < n_s; i++) {
        int ki = idx_r[i] - 1;
        if (pik_full[ki] >= 1.0) {
            for (int j = 0; j < n_s; j++) {
                int kj = idx_r[j] - 1;
                pikl[(size_t)j * n_s + i] = pik_full[kj];
                pikl[(size_t)i * n_s + j] = pik_full[kj];
            }
            pikl[(size_t)i * n_s + i] = pik_full[ki];
        }
    }

    if (N_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    /* Build valid arrays and reverse lookup */
    double *pik_valid = (double *) R_alloc(N_valid, sizeof(double));
    int *pop_to_valid = (int *) R_alloc(N_full, sizeof(int));
    for (int k = 0; k < N_full; k++) pop_to_valid[k] = -1;

    int vi = 0;
    double n_valid_sum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            pik_valid[vi] = pik_full[k];
            pop_to_valid[k] = vi;
            n_valid_sum += pik_full[k];
            vi++;
        }
    }
    const int n_valid = (int)(n_valid_sum + 0.5);

    if (n_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    /* Work in the smaller of the design and its complement (see
     * C_cps_jip). */
    const int use_complement = (n_valid > N_valid - n_valid);
    const int n_work = use_complement ? (N_valid - n_valid) : n_valid;

    double *pik_work = (double *) R_alloc(N_valid, sizeof(double));
    for (int i = 0; i < N_valid; i++) {
        pik_work[i] = use_complement ? (1.0 - pik_valid[i]) : pik_valid[i];
    }

    double *w = (double *) R_alloc(N_valid, sizeof(double));
    double *f = (double *) R_alloc(CPS_TABLE_SIZE(N_valid, n_work),
                                   sizeof(double));
    double cal_max_diff = 0.0;
    int    cal_worst_idx = -1;
    int cal_iters = cps_calibrate(pik_work, N_valid, n_work, w, f, 1e-9,
                                  500, &cal_max_diff, &cal_worst_idx);
    if (cal_iters >= 500) {
        /*
         * valid_idx isn't built here (we only have pop_to_valid). Build a
         * valid-to-pop map inline so the warning can flag 1-based population
         * indices of near-boundary units.
         */
        int *valid_to_pop = (int *) R_alloc(N_valid, sizeof(int));
        for (int k = 0; k < N_full; k++) {
            if (pop_to_valid[k] >= 0) valid_to_pop[pop_to_valid[k]] = k;
        }
        cps_warn_nonconverge("joint_inclusion_prob: ", 500, 1e-9,
                             cal_max_diff, pik_valid, N_valid, valid_to_pop);
    }
    (void) cal_worst_idx;

    const double f_n = f[n_work]; /* row 0, column n_work */

    /* Identify target units: sampled AND valid */
    int n_target = 0;
    for (int i = 0; i < n_s; i++) {
        if (pop_to_valid[idx_r[i] - 1] >= 0) n_target++;
    }

    if (n_target < 2) {
        UNPROTECT(1);
        return result;
    }

    int *tgt_out = (int *) R_alloc(n_target, sizeof(int));
    int *tgt_vi = (int *) R_alloc(n_target, sizeof(int));

    int ti = 0;
    for (int i = 0; i < n_s; i++) {
        int v = pop_to_valid[idx_r[i] - 1];
        if (v >= 0) {
            tgt_out[ti] = i;
            tgt_vi[ti] = v;
            ti++;
        }
    }

    /* Relative near-equality threshold for the Aires difference formula */
    const double tau = 16.0 * sqrt(DBL_EPSILON);

    /* Check if all weights are approximately equal */
    double w_min = w[0], w_max = w[0];
    for (int i = 1; i < N_valid; i++) {
        if (w[i] < w_min) w_min = w[i];
        if (w[i] > w_max) w_max = w[i];
    }
    double w_scale = fmax(1.0, w_max);
    int all_equal = (w_max - w_min) <= tau * w_scale;

    if (all_equal) {
        /* All weights equal => CPS reduces to SRS on valid subset */
        double pi_ij = (double)n_valid * (n_valid - 1) /
                       ((double)N_valid * (N_valid - 1));
        for (int i = 0; i < n_target; i++) {
            for (int jj = i + 1; jj < n_target; jj++) {
                int si = tgt_out[i];
                int sj = tgt_out[jj];
                pikl[(size_t)sj * n_s + si] = pi_ij;
                pikl[(size_t)si * n_s + sj] = pi_ij;
            }
        }
        UNPROTECT(1);
        return result;
    }

    /* Sort target units by calibrated weight and group exact ties.
     * Weights come from the full-population calibration; the exact
     * pair evaluation also runs over the full valid population. */
    double *wt = (double *) R_alloc(n_target, sizeof(double));
    int *ordt = (int *) R_alloc(n_target, sizeof(int));
    for (int t = 0; t < n_target; t++) {
        wt[t] = w[tgt_vi[t]];
        ordt[t] = t;
    }
    rsort_with_index(wt, ordt, n_target);

    int *gstart = (int *) R_alloc(n_target + 1, sizeof(int));
    const int G = cps_weight_groups(wt, n_target, gstart);

    double *buf = (double *) R_alloc(n_work + 1, sizeof(double));

    for (int g = 0; g < G; g++) {
        R_CheckUserInterrupt();
        for (int h = g; h < G; h++) {
            const int gs = gstart[g], ge = gstart[g + 1];
            const int hs = gstart[h], he = gstart[h + 1];
            if (g == h && ge - gs < 2) continue;

            const double w_lo = wt[gs];
            const double w_hi = wt[hs];
            const double scale = fmax(1.0, w_hi);
            const int near = (w_hi - w_lo) <= tau * scale;

            double v_block = 0.0;
            if (near) {
                const int rep_i = tgt_vi[ordt[gs]];
                const int rep_j = (g == h) ? tgt_vi[ordt[gs + 1]]
                                           : tgt_vi[ordt[hs]];
                v_block = cps_pair_exact(w, N_valid, n_work,
                                         rep_i, rep_j, f_n, buf);
            }

            for (int a = gs; a < ge; a++) {
                const int b0 = (g == h) ? a + 1 : hs;
                for (int b = b0; b < he; b++) {
                    const int ta = ordt[a];
                    const int tb = ordt[b];
                    const int via = tgt_vi[ta];
                    const int vib = tgt_vi[tb];
                    double v;
                    if (near) {
                        v = v_block;
                    } else {
                        const double r_i = w[via];
                        const double r_j = w[vib];
                        v = (r_j * pik_work[via] - r_i * pik_work[vib]) /
                            (r_j - r_i);
                    }
                    const double pi_ij = cps_pair_finalize(
                        v, pik_valid[via], pik_valid[vib], use_complement);

                    const int si = tgt_out[ta];
                    const int sj = tgt_out[tb];
                    pikl[(size_t)sj * n_s + si] = pi_ij;
                    pikl[(size_t)si * n_s + sj] = pi_ij;
                }
            }
        }
    }

    UNPROTECT(1);
    return result;
}

/*
 * Total overlap length of the circular arcs [a1, a1+l1) and [a2, a2+l2)
 * on the unit circle (0 <= a < 1, 0 <= l < 1): sum the linear overlaps
 * over integer shifts of the second arc.
 */
static double circ_overlap(double a1, double l1, double a2, double l2) {
    double tot = 0.0;
    for (int k = -1; k <= 1; k++) {
        double lo = fmax(a1, a2 + (double)k);
        double hi = fmin(a1 + l1, a2 + (double)k + l2);
        if (hi > lo) tot += hi - lo;
    }
    return tot;
}

/*
 * Exact systematic PPS joint probabilities.
 *
 * Unit k (in the order supplied) occupies the circular arc
 * [V_{k-1} mod 1, V_{k-1} mod 1 + pik_k) where V is the cumulative sum
 * of the valid pik; the unit is selected when the random start falls in
 * its arc. Two units are jointly selected exactly when the start falls
 * in both arcs, so pi_ij is the overlap length of the two arcs -- an
 * O(1) computation per pair. Total complexity O(N^2) with O(N) extra
 * memory (the previous implementation built an N x N interval indicator
 * matrix and was O(N^3)).
 */
SEXP C_up_systematic_jip(SEXP pik_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);

    SEXP result = PROTECT(allocMatrix(REALSXP, N_full, N_full));
    double *pikl = REAL(result);

    /* Defaults: independence for pairs involving certainty/zero units,
     * marginals on the diagonal. Valid-valid pairs are overwritten. */
    for (int i = 0; i < N_full; i++) {
        for (int j = 0; j < N_full; j++) {
            pikl[(size_t)j * N_full + i] = pik_full[i] * pik_full[j];
        }
        pikl[(size_t)i * N_full + i] = pik_full[i];
    }

    int N = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            N++;
        }
    }

    if (N < 2) {
        UNPROTECT(1);
        return result;
    }

    int *valid_idx = (int *) R_alloc(N, sizeof(int));
    double *len = (double *) R_alloc(N, sizeof(double));
    double *start = (double *) R_alloc(N, sizeof(double));

    int j = 0;
    double cumsum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            valid_idx[j] = k;
            len[j] = pik_full[k];
            start[j] = fmod(cumsum, 1.0);
            cumsum += pik_full[k];
            j++;
        }
    }

    for (int i = 0; i < N; i++) {
        R_CheckUserInterrupt();
        for (int jj = i + 1; jj < N; jj++) {
            double pi_ij = circ_overlap(start[i], len[i],
                                        start[jj], len[jj]);

            int k = valid_idx[i];
            int l = valid_idx[jj];
            pikl[(size_t)l * N_full + k] = pi_ij;
            pikl[(size_t)k * N_full + l] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}

/*
 * Systematic JIP submatrix for sampled units only.
 * idx_sexp: 1-based integer vector of sampled population indices.
 * Returns n_s x n_s matrix instead of N x N.
 *
 * The arc positions require one O(N) pass over the full population;
 * the pair loop is O(n_s^2) with O(1) work per pair.
 */
SEXP C_up_systematic_jip_sub(SEXP pik_sexp, SEXP idx_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);
    const int n_s = LENGTH(idx_sexp);
    const int *idx_r = INTEGER(idx_sexp); /* 1-based R indices */

    /* Allocate n_s x n_s output */
    SEXP result = PROTECT(allocMatrix(REALSXP, n_s, n_s));
    double *pikl = REAL(result);

    /* Initialize: pikl[i,j] = pik[i]*pik[j], diagonal = pik[i]
     * This handles certainty-certainty (1*1), certainty-valid (1*pik_j),
     * and zero-anything (0) pairs correctly as defaults. */
    for (int i = 0; i < n_s; i++) {
        double pi_i = pik_full[idx_r[i] - 1];
        for (int j = 0; j < n_s; j++) {
            pikl[(size_t)j * n_s + i] = pi_i * pik_full[idx_r[j] - 1];
        }
        pikl[(size_t)i * n_s + i] = pi_i;
    }

    /* Arc start of every valid unit (prefix sums over full population) */
    int N_valid = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) N_valid++;
    }

    if (N_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    double *start_pop = (double *) R_alloc(N_full, sizeof(double));
    double cumsum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            start_pop[k] = fmod(cumsum, 1.0);
            cumsum += pik_full[k];
        } else {
            start_pop[k] = -1.0; /* not a valid unit */
        }
    }

    for (int i = 0; i < n_s; i++) {
        R_CheckUserInterrupt();
        int ki = idx_r[i] - 1;
        if (start_pop[ki] < 0.0) continue;

        for (int jj = i + 1; jj < n_s; jj++) {
            int kj = idx_r[jj] - 1;
            if (start_pop[kj] < 0.0) continue;

            double pi_ij = circ_overlap(start_pop[ki], pik_full[ki],
                                        start_pop[kj], pik_full[kj]);

            pikl[(size_t)jj * n_s + i] = pi_ij;
            pikl[(size_t)i * n_s + jj] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_high_entropy_jip(SEXP pik_sexp, SEXP eps_sexp) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik_full = REAL(pik_sexp);
    const double eps = REAL(eps_sexp)[0];

    SEXP result = PROTECT(allocMatrix(REALSXP, N_full, N_full));
    double *pikl = REAL(result);
    memset(pikl, 0, (size_t)N_full * N_full * sizeof(double));

    /* Diagonal = first-order inclusion probabilities */
    for (int k = 0; k < N_full; k++) {
        pikl[k * N_full + k] = pik_full[k];
    }

    /* Handle certainty units: pi_ij = pik[other] */
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] >= 1.0 - eps) {
            for (int l = 0; l < N_full; l++) {
                pikl[k * N_full + l] = pik_full[l];
                pikl[l * N_full + k] = pik_full[l];
            }
            pikl[k * N_full + k] = pik_full[k];
        }
    }

    /* Count valid (non-certainty, non-zero) units */
    int N_valid = 0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            N_valid++;
        }
    }

    if (N_valid < 2) {
        UNPROTECT(1);
        return result;
    }

    /* Extract valid subset */
    int *valid_idx = (int *) R_alloc(N_valid, sizeof(int));
    double *pik_valid = (double *) R_alloc(N_valid, sizeof(double));

    int j = 0;
    double n = 0.0;
    double sum_pik2 = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > eps && pik_full[k] < 1.0 - eps) {
            valid_idx[j] = k;
            pik_valid[j] = pik_full[k];
            n += pik_full[k];
            sum_pik2 += pik_full[k] * pik_full[k];
            j++;
        }
    }

    if (n <= 1.0 + eps) {
        UNPROTECT(1);
        return result;
    }

    /* Compute high-entropy coefficients on valid subset */
    double *c = (double *) R_alloc(N_valid, sizeof(double));
    const double nm1 = n - 1.0;
    const double coef1 = (2.0 * n - 1.0) / nm1;
    const double coef2 = sum_pik2 / nm1;

    for (int k = 0; k < N_valid; k++) {
        double denom = n - coef1 * pik_valid[k] + coef2;
        if (fabs(denom) < 1e-15) denom = 1e-15;
        c[k] = nm1 / denom;
    }

    /* Compute joint probabilities for valid pairs, with clamping */
    for (int i = 0; i < N_valid; i++) {
        R_CheckUserInterrupt();
        for (int jj = i + 1; jj < N_valid; jj++) {
            double pi_ij = pik_valid[i] * pik_valid[jj] *
                           (c[i] + c[jj]) / 2.0;

            if (pi_ij < 0.0) pi_ij = 0.0;
            if (pi_ij > pik_valid[i]) pi_ij = pik_valid[i];
            if (pi_ij > pik_valid[jj]) pi_ij = pik_valid[jj];

            int k = valid_idx[i];
            int l = valid_idx[jj];
            pikl[l * N_full + k] = pi_ij;
            pikl[k * N_full + l] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}
