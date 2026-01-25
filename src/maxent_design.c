#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include "maxent.h"

static void compute_pi_from_lambda(const double *exp_lambda, int N, int n, 
                                   double *pi_out, double *pi_prev) {
    memset(pi_prev, 0, N * sizeof(double));
    
    for (int j = 1; j <= n; j++) {
        double denom = 0.0;
        for (int k = 0; k < N; k++) {
            denom += exp_lambda[k] * (1.0 - pi_prev[k]);
        }
        
        const double scale = (double)j / denom;
        for (int k = 0; k < N; k++) {
            pi_out[k] = exp_lambda[k] * (1.0 - pi_prev[k]) * scale;
        }
        
        if (j < n) {
            memcpy(pi_prev, pi_out, N * sizeof(double));
        }
    }
}

static int compute_lambda_newton(const double *pi_target, int N, int n,
                                 double *lambda, double *exp_lambda,
                                 double tol, int max_iter) {
    double *pi_current = (double *) R_alloc(N, sizeof(double));
    double *pi_work = (double *) R_alloc(N, sizeof(double));
    double *logit_target = (double *) R_alloc(N, sizeof(double));
    
    for (int k = 0; k < N; k++) {
        double p = pi_target[k];
        if (p < 1e-10) p = 1e-10;
        if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
        logit_target[k] = log(p / (1.0 - p));
        lambda[k] = logit_target[k];
    }
    
    for (int iter = 0; iter < max_iter; iter++) {
        for (int k = 0; k < N; k++) {
            exp_lambda[k] = exp(lambda[k]);
        }
        
        compute_pi_from_lambda(exp_lambda, N, n, pi_current, pi_work);
        
        double max_diff = 0.0;
        for (int k = 0; k < N; k++) {
            double p = pi_current[k];
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            
            double diff = logit_target[k] - log(p / (1.0 - p));
            lambda[k] += diff;
            
            if (fabs(diff) > max_diff) max_diff = fabs(diff);
        }
        
        if (max_diff < tol) {
            for (int k = 0; k < N; k++) {
                exp_lambda[k] = exp(lambda[k]);
            }
            return iter + 1;
        }
    }
    
    for (int k = 0; k < N; k++) {
        exp_lambda[k] = exp(lambda[k]);
    }
    
    return max_iter;
}

SEXP C_maxent_design_create(SEXP pik, SEXP eps) {
    const double *pik_ptr = REAL(pik);
    const double epsilon = REAL(eps)[0];
    const int N_full = LENGTH(pik);
    
    int N = 0;
    int N_certain = 0;
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
    
    const int n = (int)(n_sum + 0.5);
    
    double *pikb = NULL;
    int *idx_arr = NULL;
    int *certain_arr = NULL;
    
    if (N > 0) {
        pikb = (double *) R_alloc(N, sizeof(double));
        idx_arr = (int *) R_alloc(N, sizeof(int));
    }
    if (N_certain > 0) {
        certain_arr = (int *) R_alloc(N_certain, sizeof(int));
    }
    
    int j = 0, c = 0;
    for (int k = 0; k < N_full; k++) {
        double pk = pik_ptr[k];
        if (pk >= 1.0 - epsilon) {
            certain_arr[c++] = k;
        } else if (pk > epsilon) {
            pikb[j] = pk;
            idx_arr[j] = k;
            j++;
        }
    }
    
    double *lambda = NULL;
    double *exp_lambda_arr = NULL;
    int newton_iters = 0;
    
    if (N > 0 && n > 0) {
        lambda = (double *) R_alloc(N, sizeof(double));
        exp_lambda_arr = (double *) R_alloc(N, sizeof(double));
        newton_iters = compute_lambda_newton(pikb, N, n, lambda, exp_lambda_arr, 1e-10, 100);
        
        if (newton_iters >= 100) {
            warning("Newton method did not converge in 100 iterations");
        }
    }
    
    SEXP result = PROTECT(allocVector(VECSXP, DESIGN_LENGTH));
    SEXP names = PROTECT(allocVector(STRSXP, DESIGN_LENGTH));
    
    SEXP pik_valid_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(pik_valid_sexp), pikb, N * sizeof(double));
    SET_VECTOR_ELT(result, DESIGN_PIK_VALID, pik_valid_sexp);
    SET_STRING_ELT(names, DESIGN_PIK_VALID, mkChar("pik_valid"));
    
    SEXP lambda_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(lambda_sexp), lambda, N * sizeof(double));
    SET_VECTOR_ELT(result, DESIGN_LAMBDA, lambda_sexp);
    SET_STRING_ELT(names, DESIGN_LAMBDA, mkChar("lambda"));
    
    SEXP exp_lambda_sexp = PROTECT(allocVector(REALSXP, N));
    if (N > 0) memcpy(REAL(exp_lambda_sexp), exp_lambda_arr, N * sizeof(double));
    SET_VECTOR_ELT(result, DESIGN_EXP_LAMBDA, exp_lambda_sexp);
    SET_STRING_ELT(names, DESIGN_EXP_LAMBDA, mkChar("exp_lambda"));
    
    SEXP pi_poisson_sexp = PROTECT(allocVector(REALSXP, N));
    double *pi_poisson_ptr = REAL(pi_poisson_sexp);
    for (int k = 0; k < N; k++) {
        pi_poisson_ptr[k] = exp_lambda_arr[k] / (1.0 + exp_lambda_arr[k]);
    }
    SET_VECTOR_ELT(result, DESIGN_PI_POISSON, pi_poisson_sexp);
    SET_STRING_ELT(names, DESIGN_PI_POISSON, mkChar("pi_poisson"));
    
    SET_VECTOR_ELT(result, DESIGN_N, ScalarInteger(N));
    SET_STRING_ELT(names, DESIGN_N, mkChar("N"));
    
    SET_VECTOR_ELT(result, DESIGN_N_SAMPLE, ScalarInteger(n));
    SET_STRING_ELT(names, DESIGN_N_SAMPLE, mkChar("n"));
    
    SEXP idx_sexp = PROTECT(allocVector(INTSXP, N));
    int *idx_ptr = INTEGER(idx_sexp);
    for (int k = 0; k < N; k++) {
        idx_ptr[k] = idx_arr[k] + 1;
    }
    SET_VECTOR_ELT(result, DESIGN_IDX, idx_sexp);
    SET_STRING_ELT(names, DESIGN_IDX, mkChar("idx"));
    
    SET_VECTOR_ELT(result, DESIGN_N_FULL, ScalarInteger(N_full));
    SET_STRING_ELT(names, DESIGN_N_FULL, mkChar("N_full"));
    
    SET_VECTOR_ELT(result, DESIGN_NEWTON_ITERS, ScalarInteger(newton_iters));
    SET_STRING_ELT(names, DESIGN_NEWTON_ITERS, mkChar("newton_iters"));
    
    SET_VECTOR_ELT(result, DESIGN_N_CERTAIN, ScalarInteger(N_certain));
    SET_STRING_ELT(names, DESIGN_N_CERTAIN, mkChar("N_certain"));
    
    SEXP certain_sexp = PROTECT(allocVector(INTSXP, N_certain));
    int *certain_ptr = INTEGER(certain_sexp);
    for (int k = 0; k < N_certain; k++) {
        certain_ptr[k] = certain_arr[k] + 1;
    }
    SET_VECTOR_ELT(result, DESIGN_CERTAIN_IDX, certain_sexp);
    SET_STRING_ELT(names, DESIGN_CERTAIN_IDX, mkChar("certain_idx"));
    
    setAttrib(result, R_NamesSymbol, names);
    setAttrib(result, R_ClassSymbol, mkString("maxent_design"));
    
    UNPROTECT(8);
    return result;
}
