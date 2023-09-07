#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <stdlib.h>

SEXP up_brewer(SEXP pik, SEXP eps) {
    // Extract data from SEXP arguments
    double *pik_ptr = REAL(pik);
    double epsilon = REAL(eps)[0];
    int N_val = LENGTH(pik);

    // Initialize variables
    double size = 0.0;
    double a;
    double u;

    int ans_idx = 0;
    SEXP ans;
    PROTECT(ans = allocVector(INTSXP, N_val));  // Allocating maximum possible size
    int *ans_ptr = INTEGER(ans);

    // Random seed
    GetRNGstate();

    // Count and allocate space for pik values between eps and 1 - eps
    int count = 0;
    for (int k = 0; k < N_val; k++) {
        if (pik_ptr[k] > epsilon && pik_ptr[k] < 1 - epsilon) {
            count++;
        }
    }

    double *filtered_pik = malloc(count * sizeof(double));
    int *original_idx = malloc(count * sizeof(int));
    int *sb = calloc(count, sizeof(int));
    double *p = malloc(count * sizeof(double));
    int filtered_idx = 0;

    if (filtered_pik == NULL || original_idx == NULL || sb == NULL || p == NULL) {
        error("Memory allocation failed");
    }

    /* Pre-loop to handle pik >= 1 - eps and pik <= eps, and populate filtered_pik */
    for (int k = 0; k < N_val; k++) {
        if (pik_ptr[k] >= 1 - epsilon) {
            ans_ptr[ans_idx++] = k + 1;
        } else if (pik_ptr[k] > epsilon) {
            filtered_pik[filtered_idx] = pik_ptr[k];
            original_idx[filtered_idx] = k;
            size += pik_ptr[k];
            filtered_idx++;
        }
    }

    int n = fround(size, 0);

    /* Main loop  */
    for (int i = 1; i <= n; i++) {
        double p_sum = 0.0;
        a = 0.0;
	double n_minus_i_plus_1 = n - i + 1;

        for (int k = 0; k < count; k++) {
            a += filtered_pik[k] * sb[k];
        }

        /* Reduce repeated calculations */
        double n_minus_a = n - a;

        /* Calculate p */
        for (int k = 0; k < count; k++) {
            p[k] = ((1 - sb[k]) * filtered_pik[k] * (n_minus_a - filtered_pik[k])) / (n_minus_a - filtered_pik[k] * n_minus_i_plus_1);
            p_sum += p[k];
        }

        /* Division costs more than multiplication */
        double inv_p_sum = 1 / p_sum;

        /* Generate uniform distribution */
        u = unif_rand();

        /* Find the first index j such that p[j] - u > 0 */
        double cumprob = 0.0;
        int j;
        for (j = 0; j < count; j++) {
            cumprob += p[j] * inv_p_sum;
            if (u < cumprob) {
                break;
            }
        }

        /* Update sb and ans */
        sb[j] = 1;
        ans_ptr[ans_idx++] = original_idx[j] + 1;
    }

    /* Resize ans to actual size  */
    SEXP ans_resized;
    PROTECT(ans_resized = allocVector(INTSXP, ans_idx));
    memcpy(INTEGER(ans_resized), ans_ptr, ans_idx * sizeof(int));

    /* Free the dynamically allocated memory */
    free(filtered_pik);
    free(original_idx);
    free(sb);
    free(p);

    PutRNGstate();
    UNPROTECT(2);
    return ans_resized;
}
