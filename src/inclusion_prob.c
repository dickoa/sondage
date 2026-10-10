#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <math.h>

SEXP C_inclusion_prob(SEXP a, SEXP n) {
    const int len = length(a);
    const double n_val = asReal(n);
    const double *a_ptr = REAL(a);

    if (len == 0) {
        SEXP res = PROTECT(allocVector(REALSXP, 0));
        UNPROTECT(1);
        return res;
    }
    if (ISNA(n_val) || ISNAN(n_val)) {
        error("n must not be NA or NaN");
    }
    if (n_val < 0) {
        error("n must be non-negative");
    }

    SEXP pik = PROTECT(allocVector(REALSXP, len));
    double *pik_ptr = REAL(pik);

    /*
     * The R wrapper (inclusion_prob.default) validates that x contains no
     * NA / NaN / Inf before calling us, so we never see non-finite input
     * here and can operate on plain doubles.
     *
     * Sizes are divided by their maximum before summation so the result
     * is invariant to rescaling of x: raw sums of very large sizes would
     * overflow to Inf, and sums of subnormal sizes lose all precision.
     *
     * n may be fractional (the expected size of a random-size design).
     * The domain is 0 <= n <= n_pos, checked exactly: any tolerance would
     * let a tiny target through on all-zero sizes.
     */
    double max_a = 0.0;
    for (int i = 0; i < len; i++) {
        if (a_ptr[i] > max_a) max_a = a_ptr[i];
    }

    double sum_a = 0.0;
    int n_pos = 0;
    for (int i = 0; i < len; i++) {
        double val = a_ptr[i];
        if (val <= 0.0) {
            pik_ptr[i] = 0.0;
        } else {
            pik_ptr[i] = val / max_a;
            sum_a += pik_ptr[i];
            n_pos++;
        }
    }

    if (n_val > (double)n_pos) {
        UNPROTECT(1);
        error("'n' (%.15g) exceeds the number of units with positive size (%d)",
              n_val, n_pos);
    }
    if (n_pos == 0) {
        /* n == 0 here, and every pik is already 0 */
        UNPROTECT(1);
        return pik;
    }

    const double scale = n_val / sum_a;
    int n_capped = 0;
    double sum_uncapped = 0.0;
    /*
     * Uncapped values are a / base times a factor common to all of them.
     * max_u is the largest raw size among the uncapped units.
     */
    double base = max_a;
    double max_u = 0.0;

    for (int i = 0; i < len; i++) {
        double pi = pik_ptr[i] * scale;
        if (pi >= 1.0) {
            pik_ptr[i] = 1.0;
            n_capped++;
        } else {
            pik_ptr[i] = pi;
            sum_uncapped += pi;
            if (a_ptr[i] > max_u) max_u = a_ptr[i];
        }
    }

    /*
     * A unit is uncapped while its size is positive and its pik below 1.
     * Its pik can be 0 when its size underflowed against the base.
     */
    while (n_capped > 0 && max_u > 0.0) {
        R_CheckUserInterrupt();
        double remaining_n = n_val - (double)n_capped;

        if (max_u < 1e-200 * base) {
            /*
             * Sizes span 1e200 or more, and the units that set the base
             * are all certain. The remaining values were divided by a
             * base so large that they may be subnormal or 0, with few or
             * no significant bits, and scaling them up would amplify the
             * loss or overflow to Inf. Recompute them from the raw sizes
             * against the largest uncapped size, halved so that they stay
             * below the capped value 1. A unit still subnormal against
             * the new base is at least 1e108 times smaller than the
             * largest uncapped unit, so its probability, and its error,
             * are below 1e-108. Inputs with a smaller span never reach
             * this branch.
             */
            base = max_u;
            sum_uncapped = 0.0;
            for (int i = 0; i < len; i++) {
                if (a_ptr[i] > 0.0 && pik_ptr[i] < 1.0) {
                    pik_ptr[i] = a_ptr[i] / base / 2.0;
                    sum_uncapped += pik_ptr[i];
                }
            }
        } else if (sum_uncapped == 0.0) {
            break;
        }
        double rescale = remaining_n / sum_uncapped;

        int new_capped = 0;
        double new_sum_uncapped = 0.0;
        double new_max_u = 0.0;

        for (int i = 0; i < len; i++) {
            if (a_ptr[i] > 0.0 && pik_ptr[i] < 1.0) {
                double pi = pik_ptr[i] * rescale;
                if (pi >= 1.0) {
                    pik_ptr[i] = 1.0;
                    new_capped++;
                } else {
                    pik_ptr[i] = pi;
                    new_sum_uncapped += pi;
                    if (a_ptr[i] > new_max_u) new_max_u = a_ptr[i];
                }
            }
        }

        if (new_capped == 0) {
            break;
        }

        n_capped += new_capped;
        sum_uncapped = new_sum_uncapped;
        max_u = new_max_u;
    }

    /*
     * With n <= n_pos the capping always has a solution, and the
     * recomputation above keeps the factor finite, so this is not expected
     * to fire. It keeps a numerical failure from being returned as a
     * design. Written as !(<=) so that a NaN total fails.
     */
    double total = 0.0;
    for (int i = 0; i < len; i++) total += pik_ptr[i];
    if (!(fabs(total - n_val) <= 64.0 * DBL_EPSILON * (double)len * n_val)) {
        UNPROTECT(1);
        error("inclusion probabilities sum to %.17g, not n = %.17g "
              "(numerical failure)", total, n_val);
    }

    UNPROTECT(1);
    return pik;
}
