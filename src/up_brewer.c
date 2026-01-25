#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP C_up_brewer(SEXP pik, SEXP eps) {
    const int N_full = LENGTH(pik);
    const double *pik_ptr = REAL(pik);
    const double epsilon = REAL(eps)[0];
    
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
    
    SEXP s = PROTECT(allocVector(INTSXP, total_selected));
    int *s_ptr = INTEGER(s);
    int s_idx = 0;
    
    double *pikb = (double *) R_alloc(N, sizeof(double));
    double *cumprob = (double *) R_alloc(N, sizeof(double));
    int *orig_idx = (int *) R_alloc(N, sizeof(int));
    int *active = (int *) R_alloc(N, sizeof(int));
    
    int j = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            s_ptr[s_idx++] = k + 1;
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
    
    double a = 0.0;
    int n_active = N;
    
    for (int i = 1; i <= n_draws; i++) {
        const double m = n - a;
        const double r = (double)(n_draws - i + 1);
        
        double sum_p = 0.0;
        for (int j = 0; j < n_active; j++) {
            const int k = active[j];
            const double pk = pikb[k];
            const double denom = m - pk * r;
            
            if (denom > 1e-10) {
                double prob = pk * (m - pk) / denom;
                if (prob > 0.0) sum_p += prob;
            }
            cumprob[j] = sum_p;
        }
        
        const double u = unif_rand() * sum_p;
        int sel_j = n_active - 1;
        
        for (int j = 0; j < n_active; j++) {
            if (cumprob[j] > u) {
                sel_j = j;
                break;
            }
        }
        
        const int sel_k = active[sel_j];
        s_ptr[s_idx++] = orig_idx[sel_k] + 1;
        a += pikb[sel_k];
        
        n_active--;
        if (sel_j < n_active) {
            active[sel_j] = active[n_active];
        }
    }
    
    PutRNGstate();
    UNPROTECT(1);
    return s;
}
