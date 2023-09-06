#include <R.h>
#include <Rinternals.h>

SEXP inclusion_probs(SEXP a, SEXP n) {
    int i, l, l1;
    double sum_a = 0.0;
    int len = length(a);
    double n_val = asReal(n);
    double* a_ptr = REAL(a);
    SEXP pik1;
    PROTECT(pik1 = allocVector(REALSXP, len));
    double* pik1_ptr = REAL(pik1);
    // Calculate sum of a and correct negative values
    for (i = 0; i < len; i++) {
        if (a_ptr[i] < 0) {
            a_ptr[i] = 0;
        }
        sum_a += a_ptr[i];
    }
    // Initialize pik1
    for (i = 0; i < len; i++) {
        pik1_ptr[i] = (sum_a == 0) ? 0 : n_val * a_ptr[i] / sum_a;
    }
    // Count and adjust inclusion probabilities greater than or equal to 1
    l = 0;
    for (i = 0; i < len; i++) {
        if (pik1_ptr[i] >= 1) {
            l++;
        }
    }
    if (l > 0) {
        l1 = 0;
        while (l != l1) {
            double temp_sum = 0;
            for (i = 0; i < len; i++) {
                if (pik1_ptr[i] < 1) {
                    temp_sum += pik1_ptr[i];
                }
            }
            for (i = 0; i < len; i++) {
                pik1_ptr[i] = (pik1_ptr[i] < 1) ? (n_val - l) * (pik1_ptr[i] / temp_sum) : 1;
            }

            l1 = l;
            l = 0;
            for (i = 0; i < len; i++) {
                if (pik1_ptr[i] >= 1) {
                    l++;
                }
            }
        }
    }
    UNPROTECT(1);
    return pik1;
}
