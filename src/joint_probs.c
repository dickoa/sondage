/*
 * Joint Inclusion Probabilities
 *
 * - C_cps_jip: Exact CPS joint probabilities (Aires' formula)
 * - C_sampford_jip: Exact Sampford joint probabilities
 * - C_up_systematic_jip: Exact systematic sampling joint probabilities
 * - C_wr_systematic_jeh: Exact joint expected hits of systematic_mr
 * - C_high_entropy_jip: High entropy approximation (equation 18, Brewer & Donadio 2003)
 */

#include "cps_core.h"
#include "systematic_core.h"
#include <float.h>
#include <R_ext/Utils.h>

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

/*
 * Aires' identity for one pair, with the marginals implied by the odds.
 * The operands are put in one order, the smaller odds first, so a pair
 * gets the same value wherever it is computed. Where the compiler fuses
 * a * b - c * d into a fused multiply-add (arm64, for one), swapping the
 * two units changes the rounding, and the full and sampled-only matrices
 * differed in the twelfth digit.
 */
static double cps_aires(const double *w, const double *marg, int i, int j) {
    if (w[j] < w[i] || (w[j] == w[i] && j < i)) {
        const int t = i;
        i = j;
        j = t;
    }
    return (w[j] * marg[i] - w[i] * marg[j]) / (w[j] - w[i]);
}

/* Odds w_lo <= w_hi closer than the relative gap. */
static int cps_odds_near(double w_lo, double w_hi, double gap) {
    return (w_hi - w_lo) <= gap * w_hi;
}

/* Add one unit with selection probability p to a Poisson-binomial table
 * of the counts 0..m, truncated at m. */
static void cps_pb_add(double *buf, int m, double p) {
    for (int z = m; z >= 1; z--) {
        buf[z] = (1.0 - p) * buf[z] + p * buf[z - 1];
    }
    buf[0] *= (1.0 - p);
}

/*
 * Working-design pi of the units at odds positions a < b, given
 * acc[z] = P(z of the units before a and between a and b) and the suffix
 * table f (stride m + 1) of the positions in odds order.
 */
static double cps_sweep_pair(const double *acc, const double *f, int m,
                             const double *ws, int a, int b, double f_n) {
    const double *suf = f + (size_t)(b + 1) * (m + 1);
    double v = 0.0;
    for (int z = 0; z <= m; z++) v += acc[z] * suf[m - z];
    const double pa = ws[a] / (1.0 + ws[a]);
    const double pb = ws[b] / (1.0 + ws[b]);
    return pa * pb * v / f_n;
}

/*
 * Exact joint probabilities of the working CPS design for a set of target
 * units, written to an ld x ld output through cps_pair_finalize().
 *
 * Pairs whose calibrated odds are well apart use Aires' identity
 *
 *   pi_ij = (w_j pi_i - w_i pi_j) / (w_j - w_i)
 *
 * with the marginals pi implied by the odds. The identity is exact for
 * those, not for the targets, which calibration only meets to its
 * tolerance, and the division by the gap amplified that residue. Its
 * rounding error grows like DBL_EPSILON over the relative gap, so pairs
 * closer than cbrt(DBL_EPSILON) take the subtraction-free form
 *
 *   pi_ij = p_i p_j P(n - 2 of U without i, j) / P(n of U),
 *
 * with p = w / (1 + w). With the units in odds order, U without i and j
 * splits into the units before i, those between i and j, and those after
 * j. A running prefix table covers the first, a suffix table the last,
 * and stepping from i to j adds the middle, so a pair costs O(n) after
 * an O(N n) pass, where a fresh recurrence per pair cost O(N n). Units
 * with tied odds are exchangeable, so one value covers every pair of two
 * tie groups.
 *
 * w and marg hold the N valid units in index order. f must hold
 * CPS_TABLE_SIZE(N, n) doubles and is overwritten. f_n = P(n of U).
 * Target t is valid unit tgt[t], written at output index row[t].
 */
static void cps_working_pairs(const double *w, double *f, int N, int n,
                              const double *marg, double f_n,
                              const int *tgt, int nt, const int *row,
                              int ld, double *out,
                              const double *pik_valid, int use_complement) {
    if (nt < 2) return;
    const double gap = cbrt(DBL_EPSILON);

#define CPS_WRITE_PAIR(ta, tb, v) do {                                    \
        const int va_ = tgt[ta], vb_ = tgt[tb];                           \
        const double pi_ = cps_pair_finalize((v), pik_valid[va_],         \
                                             pik_valid[vb_],              \
                                             use_complement);             \
        out[(size_t)row[tb] * ld + row[ta]] = pi_;                        \
        out[(size_t)row[ta] * ld + row[tb]] = pi_;                        \
    } while (0)


    /* Population in odds order, its tie groups, and target positions */
    double *ws = (double *) R_alloc(N, sizeof(double));
    int *ord = (int *) R_alloc(N, sizeof(int));
    memcpy(ws, w, (size_t)N * sizeof(double));
    for (int k = 0; k < N; k++) ord[k] = k;
    rsort_with_index(ws, ord, N);
    int *pos = (int *) R_alloc(N, sizeof(int));
    for (int k = 0; k < N; k++) pos[ord[k]] = k;

    double *tp = (double *) R_alloc(nt, sizeof(double));
    int *tpo = (int *) R_alloc(nt, sizeof(int));
    for (int t = 0; t < nt; t++) {
        tp[t] = (double) pos[tgt[t]];
        tpo[t] = t;
    }
    rsort_with_index(tp, tpo, nt);

    int any_near = 0;
    if (n >= 2) {
        for (int t = 1; t < nt && !any_near; t++) {
            any_near = cps_odds_near(w[tgt[tpo[t - 1]]], w[tgt[tpo[t]]], gap);
        }
    }

    if (!any_near) {
        for (int a = 0; a < nt; a++) {
            R_CheckUserInterrupt();
            for (int b = a + 1; b < nt; b++) {
                /* No two units are drawn together when n < 2. */
                const double v =
                    (n < 2) ? 0.0 : cps_aires(w, marg, tgt[a], tgt[b]);
                CPS_WRITE_PAIR(a, b, v);
            }
        }
        return;
    }

    int *gstart = (int *) R_alloc(N + 1, sizeof(int));
    const int G = cps_weight_groups(ws, N, gstart);
    int *grp = (int *) R_alloc(N, sizeof(int));
    for (int g = 0; g < G; g++) {
        for (int k = gstart[g]; k < gstart[g + 1]; k++) grp[k] = g;
    }

    /* Suffix table in odds order: f[k * S + z] = P(z of positions k..N-1) */
    const int m = n - 2;
    const int S = m + 1;
    cps_compute_f(ws, N, m, f);

    double *pre = (double *) R_alloc(S, sizeof(double));
    double *acc = (double *) R_alloc(S, sizeof(double));
    double *val = (double *) R_alloc(G, sizeof(double));
    memset(pre, 0, (size_t)S * sizeof(double));
    memset(val, 0, (size_t)G * sizeof(double));
    pre[0] = 1.0;

    int tcur = 0;
    for (int g = 0; g < G; g++) {
        R_CheckUserInterrupt();
        const int gs = gstart[g], ge = gstart[g + 1];
        const int t0 = tcur;
        while (tcur < nt && (int) tp[tcur] < ge) tcur++;
        const int t1 = tcur;

        if (t1 > t0) {
            /* One value per near group pair (g, h) */
            int hmax = g;
            memcpy(acc, pre, (size_t)S * sizeof(double));
            if (ge - gs >= 2) {
                val[0] = cps_sweep_pair(acc, f, m, ws, gs, gs + 1, f_n);
            }
            int cur = gs + 1;
            for (int h = g + 1;
                 h < G && cps_odds_near(ws[gs], ws[gstart[h]], gap); h++) {
                const int hs = gstart[h];
                while (cur < hs) {
                    cps_pb_add(acc, m, ws[cur] / (1.0 + ws[cur]));
                    cur++;
                }
                val[h - g] = cps_sweep_pair(acc, f, m, ws, gs, hs, f_n);
                hmax = h;
            }

            for (int a = t0; a < t1; a++) {
                const int ta = tpo[a];
                for (int b = a + 1; b < nt; b++) {
                    const int tb = tpo[b];
                    const int hb = grp[pos[tgt[tb]]];
                    const double v = (hb <= hmax)
                        ? val[hb - g]
                        : cps_aires(w, marg, tgt[ta], tgt[tb]);
                    CPS_WRITE_PAIR(ta, tb, v);
                }
            }
        }

        for (int k = gs; k < ge; k++) {
            cps_pb_add(pre, m, ws[k] / (1.0 + ws[k]));
        }
    }

#undef CPS_WRITE_PAIR
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

    /* The active units sum to their count within the input tolerance, so
     * the complement design is empty and there is nothing to calibrate. */
    if (n_work == 0) {
        for (int i = 0; i < N_valid; i++) {
            for (int jj = i + 1; jj < N_valid; jj++) {
                const double pi_ij = cps_pair_finalize(
                    0.0, pik_valid[i], pik_valid[jj], use_complement);
                int k = valid_idx[i];
                int l = valid_idx[jj];
                pikl[(size_t)l * N_full + k] = pi_ij;
                pikl[(size_t)k * N_full + l] = pi_ij;
            }
        }
        UNPROTECT(1);
        return result;
    }

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

    /* The table and the marginals the calibrated odds imply, which
     * Aires' identity needs (see cps_working_pairs). */
    cps_compute_f(w, N_valid, n_work, f);
    const double f_n = f[n_work]; /* row 0, column n_work */
    double *marg = (double *) R_alloc(N_valid, sizeof(double));
    double *pro = (double *) R_alloc(2 * (size_t)n_work, sizeof(double));
    cps_compute_pik(w, f, N_valid, n_work, marg, pro);

    /* Exactly tied odds: CPS reduces to SRS on the valid units */
    double w_min = w[0], w_max = w[0];
    for (int i = 1; i < N_valid; i++) {
        if (w[i] < w_min) w_min = w[i];
        if (w[i] > w_max) w_max = w[i];
    }
    int all_equal = (w_max - w_min) <= 64.0 * DBL_EPSILON * fmax(1.0, w_max);

    if (all_equal) {
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

    int *tgt = (int *) R_alloc(N_valid, sizeof(int));
    for (int i = 0; i < N_valid; i++) tgt[i] = i;
    cps_working_pairs(w, f, N_valid, n_work, marg, f_n, tgt, N_valid,
                      valid_idx, N_full, pikl, pik_valid, use_complement);

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

    /* Empty complement design (see C_cps_jip). */
    if (n_work == 0) {
        for (int i = 0; i < n_s; i++) {
            const int vi_i = pop_to_valid[idx_r[i] - 1];
            if (vi_i < 0) continue;
            for (int j = i + 1; j < n_s; j++) {
                const int vi_j = pop_to_valid[idx_r[j] - 1];
                if (vi_j < 0) continue;
                const double pi_ij = cps_pair_finalize(
                    0.0, pik_valid[vi_i], pik_valid[vi_j], use_complement);
                pikl[(size_t)j * n_s + i] = pi_ij;
                pikl[(size_t)i * n_s + j] = pi_ij;
            }
        }
        UNPROTECT(1);
        return result;
    }

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

    /* The table and the marginals the calibrated odds imply (see
     * cps_working_pairs). */
    cps_compute_f(w, N_valid, n_work, f);
    const double f_n = f[n_work]; /* row 0, column n_work */
    double *marg = (double *) R_alloc(N_valid, sizeof(double));
    double *pro = (double *) R_alloc(2 * (size_t)n_work, sizeof(double));
    cps_compute_pik(w, f, N_valid, n_work, marg, pro);

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

    /* Exactly tied odds: CPS reduces to SRS on the valid units */
    double w_min = w[0], w_max = w[0];
    for (int i = 1; i < N_valid; i++) {
        if (w[i] < w_min) w_min = w[i];
        if (w[i] > w_max) w_max = w[i];
    }
    int all_equal = (w_max - w_min) <= 64.0 * DBL_EPSILON * fmax(1.0, w_max);

    if (all_equal) {
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

    /* Odds come from the full-population calibration, and the exact
     * near-pair form also runs over the full valid population. */
    cps_working_pairs(w, f, N_valid, n_work, marg, f_n, tgt_vi, n_target,
                      tgt_out, n_s, pikl, pik_valid, use_complement);

    UNPROTECT(1);
    return result;
}

/* Snap a pair value that rounding left just outside its valid bounds. */
static double sampford_clamp_pair(double value, double pi, double pj) {
    const double lo = fmax(0.0, pi + pj - 1.0);
    const double hi = fmin(pi, pj);
    if (value < lo && value > lo - 2e-12) return lo;
    if (value > hi && value < hi + 2e-12) return hi;
    return value;
}

/*
 * Exact Sampford joint probability for one pair, computed with an
 * elementary-symmetric-polynomial recurrence excluding that pair.
 */
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
    return sampford_clamp_pair((double)(raw / normalizer), p[ii], p[jj]);
}

/*
 * Build the polynomials E and G of all active units up to degree top, where
 * E[r] sums prod(q) and G[r] sums prod(q) * sum(p) over subsets of size r.
 * Returns 0 when a coefficient leaves the floating-point range.
 */
static int sampford_build_poly(const double *p, const long double *q, int N,
                               int top, long double *E, long double *G) {
    memset(E, 0, (size_t)(top + 1) * sizeof(long double));
    memset(G, 0, (size_t)(top + 1) * sizeof(long double));
    E[0] = 1.0L;
    for (int i = 0; i < N; i++) {
        const int t = (i + 1 < top) ? i + 1 : top;
        for (int r = t; r >= 1; r--) {
            G[r] += q[i] * (G[r - 1] + (long double)p[i] * E[r - 1]);
            E[r] += q[i] * E[r - 1];
        }
    }
    for (int r = 0; r <= top; r++) {
        if (!isfinite(E[r]) || !isfinite(G[r])) return 0;
    }
    return 1;
}

/*
 * Remove the factor of one unit from (E, G), returning degrees 0..upto.
 * Deleting forward subtracts q * e[r - 1] from E[r], which cancels once that
 * term dominates E[r], so forward deletion stops there. The remaining degrees
 * come from degree top downward, where the same step does not cancel. A
 * negative top means the polynomials stop short of their full degree, and the
 * function then returns 0 if forward deletion could not reach upto.
 */
static int sampford_deflate(const long double *E, const long double *G,
                            int top, long double qk, long double pk,
                            long double *e, long double *g, int upto) {
    e[0] = 1.0L;
    g[0] = 0.0L;
    int c = 0;
    while (c < upto) {
        const int r = c + 1;
        const long double sub = qk * e[r - 1];
        if (sub > 0.5L * E[r]) break;
        e[r] = E[r] - sub;
        g[r] = G[r] - qk * (g[r - 1] + pk * e[r - 1]);
        c = r;
    }
    if (c < upto) {
        if (top < 0) return 0;
        long double en = 0.0L, gn = 0.0L;
        for (int r = top; r >= c + 2; r--) {
            const long double ep = (E[r] - en) / qk;
            const long double gp = (G[r] - gn) / qk - pk * ep;
            if (r - 1 <= upto) {
                e[r - 1] = ep;
                g[r - 1] = gp;
            }
            en = ep;
            gn = gp;
        }
    }
    return 1;
}

/*
 * Pair value from (E1, G1), the polynomials that already exclude unit ii.
 * top1 is their full degree, or negative when they stop short of it. Returns
 * NAN when the deletion could not finish or the value is outside its bounds,
 * and the caller then uses a slower route.
 */
static double sampford_pair_fast(const double *p, const long double *q,
                                 const long double *E1, const long double *G1,
                                 int top1, int m, int ii, int jj,
                                 long double normalizer,
                                 long double *e, long double *g) {
    const int d = m - 2;
    if (!sampford_deflate(E1, G1, top1, q[jj], (long double)p[jj], e, g, d)) {
        return NAN;
    }
    const long double raw = q[ii] * q[jj] *
        (((long double)m - p[ii] - p[jj]) * e[d] - g[d]);
    const double value = (double)(raw / normalizer);
    const double lo = fmax(0.0, p[ii] + p[jj] - 1.0);
    const double hi = fmin(p[ii], p[jj]);
    const double tol = 2e-10 * fmax(1.0, hi);
    if (!R_FINITE(value) || value < lo - tol || value > hi + tol) return NAN;
    return sampford_clamp_pair(value, p[ii], p[jj]);
}

/*
 * Common implementation for the full matrix and a requested submatrix.
 * idx is NULL for the full matrix and otherwise contains 1-based indices.
 */
static SEXP sampford_jip_impl(SEXP pik_sexp, const int *idx, int n_out) {
    const int N_full = LENGTH(pik_sexp);
    const double *pik = REAL(pik_sexp);
    if (idx == NULL) n_out = N_full;

    SEXP result = PROTECT(allocMatrix(REALSXP, n_out, n_out));
    double *pikl = REAL(result);
    memset(pikl, 0, (size_t)n_out * n_out * sizeof(double));

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
        pikl[(size_t)a * n_out + a] = pik[ka];
        for (int b = 0; b < a; b++) {
            const int kb = idx ? idx[b] - 1 : b;
            double v = NAN;
            if (pik[ka] <= 0.0 || pik[kb] <= 0.0) v = 0.0;
            else if (pik[ka] >= 1.0) v = pik[kb];
            else if (pik[kb] >= 1.0) v = pik[ka];
            if (!ISNAN(v)) {
                pikl[(size_t)a * n_out + b] = v;
                pikl[(size_t)b * n_out + a] = v;
            }
        }
    }

    if (N < 2 || n_active < 1) {
        UNPROTECT(1);
        return result;
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
    for (int i = 0; i < N; i++) {
        p[i] = use_complement ? 1.0 - p_orig[i] : p_orig[i];
    }

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
                pikl[(size_t)a * n_out + b] = v;
                pikl[(size_t)b * n_out + a] = v;
            }
        }
        UNPROTECT(1);
        return result;
    }

    if (m < 2) {
        for (int a = 0; a < n_out; a++) {
            const int ka = idx ? idx[a] - 1 : a;
            const int iact = pop_to_active[ka];
            if (iact < 0) continue;
            for (int b = 0; b < a; b++) {
                const int kb = idx ? idx[b] - 1 : b;
                const int jact = pop_to_active[kb];
                if (jact < 0) continue;
                double v = use_complement ?
                    p_orig[iact] + p_orig[jact] - 1.0 : 0.0;
                if (v < 0.0 && v > -2e-12) v = 0.0;
                pikl[(size_t)a * n_out + b] = v;
                pikl[(size_t)b * n_out + a] = v;
            }
        }
        UNPROTECT(1);
        return result;
    }

    /* Scale the odds to limit growth in the ESP recurrence. */
    double max_odds = 0.0;
    for (int i = 0; i < N; i++) {
        const double o = p[i] / (1.0 - p[i]);
        if (o > max_odds) max_odds = o;
    }
    long double scaled_sum = 0.0L;
    for (int i = 0; i < N; i++) {
        scaled_sum += (long double)(p[i] / (1.0 - p[i])) / max_odds;
    }
    /* R_allocLD, not R_alloc: R_alloc only guarantees alignment for double,
       but long double requires 16-byte alignment on x86-64. */
    long double *q = R_allocLD((size_t) N);
    for (int i = 0; i < N; i++) {
        q[i] = ((long double)(p[i] / (1.0 - p[i])) / max_odds) *
               (long double)m / scaled_sum;
    }

    /*
     * Degree m is enough while forward deletion stays stable, which keeps the
     * cost at O(N * m) for a small sample from a large population. Degree N is
     * built once, the first time a pair needs the backward path, unless the
     * O(N^2) build would cost more than recomputing every pair directly.
     */
    const int d = m - 2;
    const int can_extend = (double)N <= 0.5 * (double)n_out * n_out * d;
    const size_t poly_len = (size_t)(can_extend ? N : m) + 1;
    long double *E = R_allocLD(poly_len);
    long double *G = R_allocLD(poly_len);
    if (!sampford_build_poly(p, q, N, m, E, G)) {
        UNPROTECT(1);
        error("Sampford joint probabilities exceed the floating-point range "
              "of this platform");
    }
    const long double normalizer = (long double)m * E[m] - G[m];
    if (!(normalizer > 0.0L) || !isfinite(normalizer)) {
        UNPROTECT(1);
        error("Sampford joint-probability normalizer is non-finite");
    }

    int full = 0;
    long double *e = R_allocLD(poly_len);
    long double *g = R_allocLD(poly_len);
    long double *E1 = R_allocLD(poly_len);
    long double *G1 = R_allocLD(poly_len);

    for (int a = 0; a < n_out; a++) {
        if ((a & 127) == 0) R_CheckUserInterrupt();
        const int ka = idx ? idx[a] - 1 : a;
        const int iact = pop_to_active[ka];
        if (iact < 0) continue;
        int have1 = sampford_deflate(E, G, full ? N : -1, q[iact],
                                     (long double)p[iact], E1, G1,
                                     full ? N - 1 : d);
        for (int b = 0; b < a; b++) {
            const int kb = idx ? idx[b] - 1 : b;
            const int jact = pop_to_active[kb];
            if (jact < 0) continue;
            double work_v = have1 ?
                sampford_pair_fast(p, q, E1, G1, full ? N - 1 : -1, m,
                                   iact, jact, normalizer, e, g) : NAN;
            if (ISNAN(work_v) && !full && can_extend) {
                if (!sampford_build_poly(p, q, N, N, E, G)) {
                    UNPROTECT(1);
                    error("Sampford joint probabilities exceed the "
                          "floating-point range of this platform");
                }
                full = 1;
                have1 = sampford_deflate(E, G, N, q[iact],
                                         (long double)p[iact], E1, G1, N - 1);
                work_v = sampford_pair_fast(p, q, E1, G1, N - 1, m,
                                            iact, jact, normalizer, e, g);
            }
            if (ISNAN(work_v)) {
                work_v = sampford_pair_direct(p, q, N, m, iact, jact,
                                              normalizer, e, g);
            }
            double v = work_v;
            if (use_complement) {
                v = p_orig[iact] + p_orig[jact] - 1.0 + work_v;
            }
            const double lo = fmax(0.0, p_orig[iact] + p_orig[jact] - 1.0);
            const double hi = fmin(p_orig[iact], p_orig[jact]);
            if (v < lo && v > lo - 2e-10) v = lo;
            if (v > hi && v < hi + 2e-10) v = hi;
            pikl[(size_t)a * n_out + b] = v;
            pikl[(size_t)b * n_out + a] = v;
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_sampford_jip(SEXP pik_sexp) {
    return sampford_jip_impl(pik_sexp, NULL, LENGTH(pik_sexp));
}

SEXP C_sampford_jip_sub(SEXP pik_sexp, SEXP idx_sexp) {
    return sampford_jip_impl(pik_sexp, INTEGER(idx_sexp), LENGTH(idx_sexp));
}

/*
 * Total overlap length of the raw cumulative intervals [a1, b1) and
 * [a2, b2) on the unit circle. Translate the second interval next to the
 * first, then sum the linear overlaps over the three relevant integer shifts.
 * Shared cumulative endpoints remain bit-identical, so exactly abutting arcs
 * have exactly zero overlap.
 */
static double circ_overlap(double a1, double b1, double a2, double b2) {
    double tot = 0.0;
    const double base_shift = nearbyint(a1 - a2);
    for (int k = -1; k <= 1; k++) {
        const double shift = base_shift + (double)k;
        double lo = fmax(a1, a2 + shift);
        double hi = fmin(b1, b2 + shift);
        if (hi > lo) tot += hi - lo;
    }
    return tot;
}

/*
 * Exact systematic PPS joint probabilities.
 *
 * Unit k (in the order supplied) occupies the circular projection of the
 * raw cumulative interval [V_{k-1}, V_k), where V is the cumulative sum
 * of the valid pik. The unit is selected when the random start falls in
 * its projection. Two units are jointly selected exactly when the start
 * falls in both projections, so pi_ij is their overlap length -- an
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
    double *start = (double *) R_alloc(N, sizeof(double));
    double *end = (double *) R_alloc(N, sizeof(double));

    int j = 0;
    double cumsum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            valid_idx[j] = k;
            start[j] = cumsum;
            cumsum += pik_full[k];
            end[j] = cumsum;
            j++;
        }
    }

    for (int i = 0; i < N; i++) {
        R_CheckUserInterrupt();
        for (int jj = i + 1; jj < N; jj++) {
            double pi_ij = circ_overlap(start[i], end[i],
                                        start[jj], end[jj]);

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
    double *end_pop = (double *) R_alloc(N_full, sizeof(double));
    double cumsum = 0.0;
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] > 0.0 && pik_full[k] < 1.0) {
            start_pop[k] = cumsum;
            cumsum += pik_full[k];
            end_pop[k] = cumsum;
        } else {
            start_pop[k] = -1.0; /* not a valid unit */
            end_pop[k] = -1.0;
        }
    }

    for (int i = 0; i < n_s; i++) {
        R_CheckUserInterrupt();
        int ki = idx_r[i] - 1;
        if (start_pop[ki] < 0.0) continue;

        for (int jj = i + 1; jj < n_s; jj++) {
            int kj = idx_r[jj] - 1;
            if (start_pop[kj] < 0.0) continue;

            double pi_ij = circ_overlap(start_pop[ki], end_pop[ki],
                                        start_pop[kj], end_pop[kj]);

            pikl[(size_t)jj * n_s + i] = pi_ij;
            pikl[(size_t)i * n_s + jj] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}

/*
 * Exact joint expected hits of systematic PPS with minimum replacement.
 *
 * With hits_k = m_k + I_k, where I is the walk over the scaled
 * fractional parts, E(I_i) is the length g_i of unit i's arc and
 * E(I_i I_j) the overlap pi_ij of the two arcs. Then
 *
 *   E(hits_i hits_j) = m_i m_j + m_i g_j + m_j g_i + pi_ij   (i != j)
 *   E(hits_i^2)      = m_i^2 + 2 m_i g_i + g_i.
 *
 * The expanded form matters. Writing e_i e_j + pi_ij - g_i g_j instead
 * subtracts two equal products when m_i = m_j = 0 and leaves residues of
 * order 1e-18 where the answer is exactly 0, which hides structural
 * zeros from sampling_cov(weighted = TRUE).
 *
 * The split, the scaling and the arc ends come from systematic_core.h,
 * as in the draw, so both agree on every m_k and every arc. idx (1-based,
 * length n_out) selects the sampled units, or is NULL for all N.
 * O(N + n_out^2) time.
 */
static SEXP wr_systematic_jeh_impl(SEXP e_sexp, SEXP n_sexp,
                                   const int *idx, int n_out) {
    const int N = LENGTH(e_sexp);
    const double *e = REAL(e_sexp);
    const int n = asInteger(n_sexp);

    int *m = (int *) R_alloc(N, sizeof(int));
    double *f = (double *) R_alloc(N, sizeof(double));
    double *lo = (double *) R_alloc(N, sizeof(double));
    double *hi = (double *) R_alloc(N, sizeof(double));

    const int n_frac = sys_split_hits(e, N, n, m, f);
    if (n_frac < 0)
        error("systematic_mr: expected hits do not add up to %d", n);

    /* Arc of every active unit, built as the walk builds it */
    const int last = sys_last_active(f, N);
    double cs = 0.0;
    for (int k = 0; k < N; k++) {
        if (f[k] > 0.0) {
            lo[k] = cs;
            cs = sys_arc_end(cs, f, k, last, n_frac);
            hi[k] = cs;
        } else {
            lo[k] = 0.0;
            hi[k] = 0.0;
        }
    }

    SEXP result = PROTECT(allocMatrix(REALSXP, n_out, n_out));
    double *out = REAL(result);

    for (int a = 0; a < n_out; a++) {
        R_CheckUserInterrupt();
        const int i = idx ? idx[a] - 1 : a;
        const double mi = (double) m[i];
        const double gi = hi[i] - lo[i];

        out[(size_t) a * n_out + a] = mi * mi + 2.0 * mi * gi + gi;

        for (int b = a + 1; b < n_out; b++) {
            const int j = idx ? idx[b] - 1 : b;
            const double mj = (double) m[j];
            const double gj = hi[j] - lo[j];
            const double pij = (gi > 0.0 && gj > 0.0)
                ? circ_overlap(lo[i], hi[i], lo[j], hi[j])
                : 0.0;
            const double v = mi * mj + mi * gj + mj * gi + pij;
            out[(size_t) b * n_out + a] = v;
            out[(size_t) a * n_out + b] = v;
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_wr_systematic_jeh(SEXP e_sexp, SEXP n_sexp) {
    return wr_systematic_jeh_impl(e_sexp, n_sexp, NULL, LENGTH(e_sexp));
}

SEXP C_wr_systematic_jeh_sub(SEXP e_sexp, SEXP n_sexp, SEXP idx_sexp) {
    return wr_systematic_jeh_impl(e_sexp, n_sexp, INTEGER(idx_sexp),
                                  LENGTH(idx_sexp));
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
        pikl[(size_t)k * N_full + k] = pik_full[k];
    }

    /* Certainty units: pi_ij is the smaller probability of the pair, so
     * two near-certainty units do not depend on the processing order. */
    for (int k = 0; k < N_full; k++) {
        if (pik_full[k] >= 1.0 - eps) {
            for (int l = 0; l < N_full; l++) {
                const double v = fmin(pik_full[k], pik_full[l]);
                pikl[(size_t)k * N_full + l] = v;
                pikl[(size_t)l * N_full + k] = v;
            }
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
            pikl[(size_t)l * N_full + k] = pi_ij;
            pikl[(size_t)k * N_full + l] = pi_ij;
        }
    }

    UNPROTECT(1);
    return result;
}
