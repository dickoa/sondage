#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>

/*
 * Maximum Entropy Design Creation
 * 
 * Creates design objects for maximum entropy sampling (also known as
 * Conditional Poisson Sampling). The design object contains precomputed
 * lambda parameters that enable fast repeated sampling.
 * 
 * Key optimizations:
 * 1. Cache exp(lambda) per Newton iteration
 * 2. Precompute logit(pi_target) once
 * 3. Reuse final exp(lambda) for design object
 */

/* ============================================================
 * Compute π(λ, Sn) from λ with CACHED exp(λ)
 * 
 * exp_lambda must be precomputed: exp_lambda[k] = exp(lambda[k])
 * This eliminates (n-1) × N exp() calls compared to recomputing each iteration!
 * ============================================================ */
static void compute_pi_from_lambda_cached(const double *exp_lambda, int N, int n, 
                                          double *pi_out, double *pi_prev) {
    /* Initialize: π(λ, S0) = 0 */
    memset(pi_prev, 0, N * sizeof(double));
    
    for (int j = 1; j <= n; j++) {
        /* Compute denominator: Σ exp(λ_k) × (1 - π_k^{j-1}) */
        double denom = 0.0;
        for (int k = 0; k < N; k++) {
            denom += exp_lambda[k] * (1.0 - pi_prev[k]);
        }
        
        /* Compute π^j */
        const double scale = (double)j / denom;
        for (int k = 0; k < N; k++) {
            pi_out[k] = exp_lambda[k] * (1.0 - pi_prev[k]) * scale;
        }
        
        /* Copy to prev for next iteration (memcpy is SIMD-optimized, very fast) */
        if (j < n) {
            memcpy(pi_prev, pi_out, N * sizeof(double));
        }
    }
}

/* ============================================================
 * Optimized Newton method with cached exp(λ) and precomputed logit(π_target)
 * 
 * Returns: exp_lambda (caller provides buffer, gets exp(λ) at convergence)
 * ============================================================ */
static int compute_lambda_newton_v2(const double *pi_target, int N, int n,
                                    double *lambda, double *exp_lambda,
                                    double tol, int max_iter) {
    /* Allocate working arrays */
    double *pi_current = (double *) R_alloc(N, sizeof(double));
    double *pi_work = (double *) R_alloc(N, sizeof(double));
    double *logit_target = (double *) R_alloc(N, sizeof(double));
    
    /* Precompute logit(π_target) - done ONCE, not per iteration */
    for (int k = 0; k < N; k++) {
        double p = pi_target[k];
        /* Clamp for numerical stability */
        if (p < 1e-10) p = 1e-10;
        if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
        logit_target[k] = log(p / (1.0 - p));
        /* Initialize λ(0) = logit(π) */
        lambda[k] = logit_target[k];
    }
    
    /* Newton iterations */
    for (int iter = 0; iter < max_iter; iter++) {
        /* Compute exp(λ) once per iteration */
        for (int k = 0; k < N; k++) {
            exp_lambda[k] = exp(lambda[k]);
        }
        
        /* Compute π(λ, Sn) using cached exp(λ) */
        compute_pi_from_lambda_cached(exp_lambda, N, n, pi_current, pi_work);
        
        /* Update: λ_new = λ + logit(π_target) - logit(π_current) */
        double max_diff = 0.0;
        for (int k = 0; k < N; k++) {
            double p_current = pi_current[k];
            
            /* Clamp for numerical stability */
            if (p_current < 1e-10) p_current = 1e-10;
            if (p_current > 1.0 - 1e-10) p_current = 1.0 - 1e-10;
            
            double logit_current = log(p_current / (1.0 - p_current));
            double diff = logit_target[k] - logit_current;
            lambda[k] += diff;
            
            double abs_diff = fabs(diff);
            if (abs_diff > max_diff) max_diff = abs_diff;
        }
        
        /* Check convergence */
        if (max_diff < tol) {
            /* Recompute exp(λ) for final converged λ */
            for (int k = 0; k < N; k++) {
                exp_lambda[k] = exp(lambda[k]);
            }
            return iter + 1;
        }
    }
    
    /* Final exp(λ) if we hit max_iter */
    for (int k = 0; k < N; k++) {
        exp_lambda[k] = exp(lambda[k]);
    }
    
    return max_iter;
}

/* ============================================================
 * Create maximum entropy design object
 * 
 * Computes lambda parameters from inclusion probabilities using
 * Newton's method. Handles certainty selections (pik >= 1-eps).
 * ============================================================ */
SEXP C_maxent_design_create(SEXP pik, SEXP eps) {
    const double *pik_ptr = REAL(pik);
    const double epsilon = REAL(eps)[0];
    const int N_full = LENGTH(pik);
    
    /* Count valid units, certainty selections, and compute n */
    int N = 0;           /* Units with eps < pik < 1-eps */
    int N_certain = 0;   /* Units with pik >= 1-eps */
    double n_sum = 0.0;
    
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            N_certain++;
        } else if (pk > epsilon) {
            N++;
            n_sum += pk;
        }
    }
    
    /* Sample size for random units */
    const int n = (int)(n_sum + 0.5);
    
    /* Extract valid units and certainty selections */
    double *pikb = NULL;
    int *idx_arr = NULL;
    int *certain_idx = NULL;
    
    if (N > 0) {
        pikb = (double *) R_alloc(N, sizeof(double));
        idx_arr = (int *) R_alloc(N, sizeof(int));
    }
    if (N_certain > 0) {
        certain_idx = (int *) R_alloc(N_certain, sizeof(int));
    }
    
    int j = 0, c = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            certain_idx[c++] = k;
        } else if (pk > epsilon) {
            pikb[j] = pk;
            idx_arr[j] = k;
            j++;
        }
    }
    
    /* Allocate λ and exp(λ) */
    double *lambda = NULL;
    double *exp_lambda_arr = NULL;
    int newton_iters = 0;
    
    if (N > 0 && n > 0) {
        lambda = (double *) R_alloc(N, sizeof(double));
        exp_lambda_arr = (double *) R_alloc(N, sizeof(double));
        
        /* Compute λ from π using optimized Newton */
        newton_iters = compute_lambda_newton_v2(pikb, N, n, lambda, exp_lambda_arr, 
                                                 1e-10, 100);
        
        if (newton_iters >= 100) {
            warning("Newton method did not converge in 100 iterations");
        }
    }
    
    /* Create output list - now with 11 elements */
    SEXP result = PROTECT(allocVector(VECSXP, 11));
    SEXP names = PROTECT(allocVector(STRSXP, 11));
    
    /* pik_valid */
    SEXP pik_valid_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(pik_valid_sexp), pikb, N * sizeof(double));
    SET_VECTOR_ELT(result, 0, pik_valid_sexp);
    SET_STRING_ELT(names, 0, mkChar("pik_valid"));
    
    /* lambda */
    SEXP lambda_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(lambda_sexp), lambda, N * sizeof(double));
    SET_VECTOR_ELT(result, 1, lambda_sexp);
    SET_STRING_ELT(names, 1, mkChar("lambda"));
    
    /* exp_lambda */
    SEXP exp_lambda = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(exp_lambda), exp_lambda_arr, N * sizeof(double));
    SET_VECTOR_ELT(result, 2, exp_lambda);
    SET_STRING_ELT(names, 2, mkChar("exp_lambda"));
    
    /* pi_poisson = exp(λ) / (1 + exp(λ)) */
    SEXP pi_poisson = PROTECT(allocVector(REALSXP, N));
    double *pi_poisson_ptr = REAL(pi_poisson);
    for (int k = 0; k < N; k++) {
        pi_poisson_ptr[k] = exp_lambda_arr[k] / (1.0 + exp_lambda_arr[k]);
    }
    SET_VECTOR_ELT(result, 3, pi_poisson);
    SET_STRING_ELT(names, 3, mkChar("pi_poisson"));
    
    /* N (random units) */
    SET_VECTOR_ELT(result, 4, ScalarInteger(N));
    SET_STRING_ELT(names, 4, mkChar("N"));
    
    /* n (sample size for random units) */
    SET_VECTOR_ELT(result, 5, ScalarInteger(n));
    SET_STRING_ELT(names, 5, mkChar("n"));
    
    /* idx (1-indexed for R) */
    SEXP idx = PROTECT(allocVector(INTSXP, N));
    int *idx_ptr = INTEGER(idx);
    for (int k = 0; k < N; k++) {
        idx_ptr[k] = idx_arr[k] + 1;
    }
    SET_VECTOR_ELT(result, 6, idx);
    SET_STRING_ELT(names, 6, mkChar("idx"));
    
    /* N_full */
    SET_VECTOR_ELT(result, 7, ScalarInteger(N_full));
    SET_STRING_ELT(names, 7, mkChar("N_full"));
    
    /* newton_iters */
    SET_VECTOR_ELT(result, 8, ScalarInteger(newton_iters));
    SET_STRING_ELT(names, 8, mkChar("newton_iters"));
    
    /* N_certain (new) */
    SET_VECTOR_ELT(result, 9, ScalarInteger(N_certain));
    SET_STRING_ELT(names, 9, mkChar("N_certain"));
    
    /* certain_idx (new, 1-indexed for R) */
    SEXP certain_idx_sexp = PROTECT(allocVector(INTSXP, N_certain));
    int *certain_idx_ptr = INTEGER(certain_idx_sexp);
    for (int k = 0; k < N_certain; k++) {
        certain_idx_ptr[k] = certain_idx[k] + 1;
    }
    SET_VECTOR_ELT(result, 10, certain_idx_sexp);
    SET_STRING_ELT(names, 10, mkChar("certain_idx"));
    
    setAttrib(result, R_NamesSymbol, names);
    setAttrib(result, R_ClassSymbol, mkString("maxent_design"));
    
    UNPROTECT(8);
    return result;
}
