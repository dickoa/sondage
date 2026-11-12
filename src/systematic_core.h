/*
 * Shared steps of the systematic PPS walk.
 *
 * Used by the without-replacement draw (C_up_systematic), the
 * minimum-replacement draw (C_up_wr_systematic and its batch twin) and
 * the exact joint expected hits of the latter (C_wr_systematic_jeh*).
 * Keeping the arithmetic in one place makes the draw and its moments
 * agree on every unit's whole and fractional part, and makes the two
 * draws select the same units when every expected hit is below 1.
 */

#ifndef SONDAGE_SYSTEMATIC_CORE_H
#define SONDAGE_SYSTEMATIC_CORE_H

#include <float.h>
#include <limits.h>
#include <math.h>

/*
 * An expected hit within this relative distance of a whole number of at
 * least 1 is treated as that number. It keeps an exact whole number
 * whole after the n * prob round trip, where it can come back a few ulps
 * low. Values near 0 are left alone, since 0 * prob is exactly 0.
 */
#define SYS_SNAP_REL (64.0 * DBL_EPSILON)

/*
 * Scale the positive entries of v, all below 1, so that they sum to
 * n_sel. An entry that reaches 1 after scaling is certain. It is set to
 * 0, bump[k] is incremented, n_sel drops by one, and the remaining
 * entries are scaled again. Without such an entry, each value is the
 * input times n_sel / sum, the arithmetic of the original walk.
 *
 * Returns the final n_sel (0 when nothing is left to draw), or -1 when
 * n_sel > 0 and no positive entry remains.
 */
static inline int sys_rescale(double *v, int N, int n_sel, int *bump) {
    while (n_sel > 0) {
        double s = 0.0;
        for (int k = 0; k < N; k++) {
            if (v[k] > 0.0) s += v[k];
        }
        if (!(s > 0.0)) return -1;

        const double scale = (double)n_sel / s;
        int promoted = 0;
        for (int k = 0; k < N; k++) {
            if (v[k] > 0.0) {
                v[k] *= scale;
                if (v[k] >= 1.0) {
                    v[k] = 0.0;
                    bump[k]++;
                    promoted++;
                }
            }
        }
        if (promoted == 0) return n_sel;
        n_sel -= promoted;
    }
    for (int k = 0; k < N; k++) v[k] = 0.0;
    return 0;
}

/* Index of the last unit with a positive value, or -1. */
static inline int sys_last_active(const double *v, int N) {
    for (int k = N - 1; k >= 0; k--) {
        if (v[k] > 0.0) return k;
    }
    return -1;
}

/*
 * End of unit k's arc on the walk: the running sum plus v[k], except
 * that the last active unit ends at exactly n_sel. The pin keeps
 * rounding in the running sum from dropping or adding a selection at
 * the end of the walk.
 */
static inline double sys_arc_end(double cs, const double *v, int k,
                                 int last, int n_sel) {
    return (k == last) ? (double)n_sel : cs + v[k];
}

/*
 * One systematic walk over the positive entries of v (scaled to sum to
 * n_sel by sys_rescale) from the start u in (0, 1). Unit k is selected
 * when its arc contains a point u + j, which is the crossing test
 * floor(end - u) > floor(start - u). Sets sel[k] to 0 or 1 and returns
 * the number of units selected.
 */
static inline int sys_walk(const double *v, int N, int n_sel, int last,
                           double u, int *sel) {
    int count = 0;
    double cs = 0.0;
    for (int k = 0; k < N; k++) {
        sel[k] = 0;
        if (!(v[k] > 0.0)) continue;
        const double prev = cs;
        cs = sys_arc_end(cs, v, k, last, n_sel);
        if (floor(cs - u) > floor(prev - u)) {
            sel[k] = 1;
            count++;
        }
    }
    return count;
}

/*
 * Split expected hits e (summing to n within the input tolerance) into
 * whole parts m and fractional parts f, then scale f with sys_rescale so
 * that it sums to n_frac = n - sum(m). A fractional part that rescaling
 * pushes to 1 moves into m.
 *
 * Returns n_frac (0 when the draw is deterministic), or -1 when the
 * whole parts exceed n or no fractional mass is left to carry n_frac.
 */
static inline int sys_split_hits(const double *e, int N, int n, int *m,
                                 double *f) {
    double sum_m = 0.0;
    for (int k = 0; k < N; k++) {
        double x = e[k];
        const double r = nearbyint(x);
        if (r >= 1.0 && fabs(x - r) <= SYS_SNAP_REL * x) x = r;
        /* The R callers refuse such inputs, so this guards direct calls:
         * the cast below is undefined above INT_MAX and for NaN. */
        if (!(x >= 0.0 && x <= (double)INT_MAX)) return -1;
        const double whole = floor(x);
        m[k] = (int)whole;
        f[k] = x - whole;
        sum_m += whole;
    }
    if (sum_m > (double)n) return -1;

    const int n_frac = n - (int)sum_m;
    if (n_frac == 0) {
        for (int k = 0; k < N; k++) f[k] = 0.0;
        return 0;
    }
    return sys_rescale(f, N, n_frac, m);
}

#endif /* SONDAGE_SYSTEMATIC_CORE_H */
