#include <R.h>
#include <Rinternals.h>
#include <math.h>

SEXP C_inclusion_prob(SEXP a, SEXP n) {
    const int len = length(a);
    const double n_val = asReal(n);
    const double *a_ptr = REAL(a);

    if (len == 0) {
        return allocVector(REALSXP, 0);
    }
    if (ISNA(n_val) || ISNAN(n_val)) {
        error("n must not be NA or NaN");
    }
    if (n_val < 0) {
        error("n must be non-negative");
    }

    SEXP pik = PROTECT(allocVector(REALSXP, len));
    double *pik_ptr = REAL(pik);

    double sum_a = 0.0;
    for (int i = 0; i < len; i++) {
        double val = a_ptr[i];
        if (ISNA(val) || ISNAN(val)) {
            pik_ptr[i] = NA_REAL;
        } else if (val <= 0.0) {
            pik_ptr[i] = 0.0;
        } else {
            pik_ptr[i] = val;
            sum_a += val;
        }
    }

    if (sum_a == 0.0) {
        for (int i = 0; i < len; i++) {
            if (!ISNA(pik_ptr[i])) {
                pik_ptr[i] = 0.0;
            }
        }
        UNPROTECT(1);
        return pik;
    }

    const double scale = n_val / sum_a;
    int n_capped = 0;
    double sum_uncapped = 0.0;

    for (int i = 0; i < len; i++) {
        if (!ISNA(pik_ptr[i])) {
            double pi = pik_ptr[i] * scale;
            if (pi >= 1.0) {
                pik_ptr[i] = 1.0;
                n_capped++;
            } else {
                pik_ptr[i] = pi;
                sum_uncapped += pi;
            }
        }
    }

    while (n_capped > 0 && sum_uncapped > 0.0) {
        double remaining_n = n_val - (double)n_capped;
        double rescale = remaining_n / sum_uncapped;

        int new_capped = 0;
        double new_sum_uncapped = 0.0;

        for (int i = 0; i < len; i++) {
            if (!ISNA(pik_ptr[i]) && pik_ptr[i] < 1.0) {
                double pi = pik_ptr[i] * rescale;
                if (pi >= 1.0) {
                    pik_ptr[i] = 1.0;
                    new_capped++;
                } else {
                    pik_ptr[i] = pi;
                    new_sum_uncapped += pi;
                }
            }
        }

        if (new_capped == 0) {
            break;
        }

        n_capped += new_capped;
        sum_uncapped = new_sum_uncapped;
    }

    UNPROTECT(1);
    return pik;
}
