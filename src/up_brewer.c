#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
 * up_brewer - Highly optimized Brewer's method
 * 
 * Implements Algorithm 6.10 from Tillé's "Sampling Algorithms"
 * Returns integer vector of selected indices (1-indexed for R)
 * 
 * Optimizations:
 * 1. Packed array: only iterate over unselected units
 * 2. Swap-to-end: O(1) removal of selected units
 * 3. Incremental 'a' update: O(1) instead of O(N)
 * 4. Fused probability + cumsum in single pass
 * 5. No branch in inner probability loop
 */

SEXP C_up_brewer(SEXP pik, SEXP eps) {
    const int N_full = LENGTH(pik);
    const double *pik_ptr = REAL(pik);
    const double epsilon = REAL(eps)[0];
    
    /* Single pass: count valid units, certainty selections, and compute n */
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
    
    /* Allocate output - indices */
    SEXP s;
    PROTECT(s = allocVector(INTSXP, total_selected));
    int *s_ptr = INTEGER(s);
    int s_idx = 0;
    
    /* Working arrays */
    double *pikb = (double *) R_alloc(N, sizeof(double));     /* Inclusion probs */
    double *cumprob = (double *) R_alloc(N, sizeof(double));  /* Cumulative probs */
    int *orig_idx = (int *) R_alloc(N, sizeof(int));          /* Original indices */
    int *active = (int *) R_alloc(N, sizeof(int));            /* Active unit indices (0..N-1) */
    
    /* Extract valid units and certainty selections */
    int j = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            /* Certainty selection */
            s_ptr[s_idx++] = k + 1;  /* 1-indexed for R */
        } else if (pk > epsilon) {
            pikb[j] = pk;
            orig_idx[j] = k;
            active[j] = j;
            j++;
        }
    }
    
    if (n_draws == 0) {
        UNPROTECT(1);
        return s;
    }
    
    GetRNGstate();
    
    double a = 0.0;        /* Sum of pik for selected units */
    int n_active = N;      /* Number of active (unselected) units */
    
    for (int i = 1; i <= n_draws; i++) {
        const double m = n - a;                    /* Remaining mass */
        const double r = (double)(n_draws - i + 1); /* Remaining draws */
        
        /* Compute cumulative probabilities for active units only */
        double sum_p = 0.0;
        for (int j = 0; j < n_active; j++) {
            const int k = active[j];
            const double pk = pikb[k];
            const double denom = m - pk * r;
            
            /* Brewer's formula: p_k ∝ π_k * (m - π_k) / (m - π_k * r) */
            if (denom > 1e-10) {
                double prob = pk * (m - pk) / denom;
                if (prob > 0.0) sum_p += prob;
            }
            cumprob[j] = sum_p;
        }
        
        /* Select unit: find first j where cumprob[j] > u */
        const double u = unif_rand() * sum_p;
        int sel_j = n_active - 1;  /* Default: last active */
        
        for (int j = 0; j < n_active; j++) {
            if (cumprob[j] > u) {
                sel_j = j;
                break;
            }
        }
        
        /* Get the actual unit index */
        const int sel_k = active[sel_j];
        
        /* Store selected index (1-indexed for R) */
        s_ptr[s_idx++] = orig_idx[sel_k] + 1;
        
        /* Update a */
        a += pikb[sel_k];
        
        /* Remove from active set by swapping with last active unit */
        n_active--;
        if (sel_j < n_active) {
            active[sel_j] = active[n_active];
        }
    }
    
    PutRNGstate();
    
    UNPROTECT(1);
    return s;
}
