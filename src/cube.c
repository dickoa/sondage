/*
 * cube.c - Cube method for balanced sampling
 *
 * C port of BalancedSampling cube method for the sondage package.
 * Original C++ implementation by Wilmer Prentius (GPL >= 2).
 * https://CRAN.R-project.org/package=BalancedSampling
 *
 * References:
 * - Deville, J.C. and Tillé, Y. (2004). Efficient balanced sampling:
 *   the cube method. Biometrika, 91(4), 893-912.
 * - Chauvet, G. and Tillé, Y. (2006). A fast algorithm for balanced
 *   sampling. Computational Statistics, 21(1), 53-62.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* =========================================================================
 * Utility macros
 * ========================================================================= */

/* Row-major indexing (for RREF working matrix) */
#define IDX_RM(row, col, ncols) ((row) * (ncols) + (col))

/* Column-major indexing (for R matrices) */
#define IDX_CM(row, col, nrows) ((col) * (nrows) + (row))

/* Check if probability is effectively 0 or 1 */
#define PROB_IS_ZERO(p, eps) ((p) < (eps))
#define PROB_IS_ONE(p, eps)  ((p) > 1.0 - (eps))
#define PROB_IS_INT(p, eps)  (PROB_IS_ZERO(p, eps) || PROB_IS_ONE(p, eps))

/* =========================================================================
 * Internal types
 * ========================================================================= */

typedef struct {
    int *list;      /* Active indices: list[0..len-1] are valid */
    int *reverse;   /* reverse[id] = position in list, or capacity if removed */
    int len;        /* Current number of active indices */
    int capacity;   /* Maximum capacity (= N) */
} IndexList;

typedef struct {
    int N;              /* Population size */
    int p;              /* Number of balancing variables */
    double eps;         /* Numerical tolerance */
    double *prob;       /* Working probabilities [N] */
    double *amat;       /* A = X / pi, column-major [N x p] */
    double *bmat;       /* Working matrix, row-major [(p+1) x p] for RREF */
    double *uvec;       /* Null space vector [p+1] */
    int *candidates;    /* Current candidate indices [p+1] */
    IndexList *idx;     /* Active index list */
} CubeWorkspace;

static IndexList *indexlist_alloc(int capacity) {
    IndexList *il = (IndexList *)R_alloc(1, sizeof(IndexList));
    il->list = (int *)R_alloc(capacity, sizeof(int));
    il->reverse = (int *)R_alloc(capacity, sizeof(int));
    il->len = 0;
    il->capacity = capacity;
    return il;
}

static void indexlist_clear(IndexList *il) {
    il->len = 0;
    for (int i = 0; i < il->capacity; i++) {
        il->reverse[i] = il->capacity;  /* Mark all as removed */
    }
}

static void indexlist_add(IndexList *il, int id) {
    il->list[il->len] = id;
    il->reverse[id] = il->len;
    il->len++;
}

static void indexlist_remove(IndexList *il, int id) {
    int k = il->reverse[id];
    if (k >= il->len) return;  /* Already removed */

    il->len--;
    il->reverse[id] = il->capacity;  /* Mark as removed */

    if (k < il->len) {
        int last_id = il->list[il->len];
        il->list[k] = last_id;
        il->reverse[last_id] = k;
    }
}

static int indexlist_get(IndexList *il, int k) {
    return il->list[k];
}

static void indexlist_shuffle(IndexList *il) {
    for (int i = 0; i < il->len - 1; i++) {
        int j = i + (int)(unif_rand() * (il->len - i));
        if (j >= il->len) j = il->len - 1;  /* Clamp */
        if (i != j) {
            int tmp = il->list[i];
            il->list[i] = il->list[j];
            il->list[j] = tmp;
            il->reverse[il->list[i]] = i;
            il->reverse[il->list[j]] = j;
        }
    }
}

static void cube_rref(double *mat, int nrows, int ncols) {
    int lead = 0;

    for (int r = 0; r < nrows; r++) {
        if (lead >= ncols)
            return;

        int i = r;

        /* Find non-zero pivot in column 'lead' */
        while (mat[IDX_RM(i, lead, ncols)] == 0.0) {
            i++;
            if (i == nrows) {
                i = r;
                lead++;
                if (lead >= ncols)
                    return;
            }
        }

        if (i != r) {
            for (int k = 0; k < ncols; k++) {
                double tmp = mat[IDX_RM(i, k, ncols)];
                mat[IDX_RM(i, k, ncols)] = mat[IDX_RM(r, k, ncols)];
                mat[IDX_RM(r, k, ncols)] = tmp;
            }
        }

        double pivot = mat[IDX_RM(r, lead, ncols)];
        mat[IDX_RM(r, lead, ncols)] = 1.0;
        for (int k = lead + 1; k < ncols; k++) {
            mat[IDX_RM(r, k, ncols)] /= pivot;
        }

        for (int j = 0; j < nrows; j++) {
            if (j == r) continue;

            double factor = mat[IDX_RM(j, lead, ncols)];
            mat[IDX_RM(j, lead, ncols)] = 0.0;
            for (int k = lead + 1; k < ncols; k++) {
                mat[IDX_RM(j, k, ncols)] -= mat[IDX_RM(r, k, ncols)] * factor;
            }
        }

        lead++;
    }
}

static void cube_null_vector(double *uvec, const double *mat, int ncols) {
    int nrows = ncols - 1;

    if (mat[IDX_RM(nrows - 1, nrows - 1, ncols)] == 1.0) {
        uvec[ncols - 1] = 1.0;
        for (int i = 0; i < nrows; i++) {
            uvec[i] = -mat[IDX_RM(i, ncols - 1, ncols)];
        }
        return;
    }

    for (int k = 1; k < ncols; k++) {
        uvec[k] = (k % 2 == 0) ? -1.0 : 1.0;
    }

    for (int i = 0; i < nrows; i++) {
        int lead = 0;

        for (; lead < ncols; lead++) {
            if (mat[IDX_RM(i, lead, ncols)] == 1.0) break;
        }

        if (lead >= ncols) continue;

        uvec[lead] = 0.0;
        for (int k = lead + 1; k < ncols; k++) {
            uvec[lead] -= uvec[k] * mat[IDX_RM(i, k, ncols)];
        }
    }
}

CubeWorkspace *cube_workspace_alloc(int N, int p, double eps) {
    CubeWorkspace *ws = (CubeWorkspace *)R_alloc(1, sizeof(CubeWorkspace));

    ws->N = N;
    ws->p = p;
    ws->eps = eps;

    ws->prob = (double *)R_alloc(N, sizeof(double));
    ws->amat = (double *)R_alloc(N * p, sizeof(double));
    ws->bmat = (double *)R_alloc((p + 1) * p, sizeof(double));
    ws->uvec = (double *)R_alloc(p + 1, sizeof(double));
    ws->candidates = (int *)R_alloc(p + 1, sizeof(int));

    ws->idx = indexlist_alloc(N);

    return ws;
}

void cube_workspace_free(CubeWorkspace *ws) {
    /* R_alloc memory is automatically freed, nothing to do */
    (void)ws;
}

void cube_workspace_init(CubeWorkspace *ws, const double *prob, const double *X) {
    int N = ws->N;
    int p = ws->p;
    double eps = ws->eps;

    memcpy(ws->prob, prob, N * sizeof(double));

    for (int j = 0; j < p; j++) {
        for (int i = 0; i < N; i++) {
            double pi = ws->prob[i];
            if (pi > eps) {
                ws->amat[IDX_CM(i, j, N)] = X[IDX_CM(i, j, N)] / pi;
            } else {
                ws->amat[IDX_CM(i, j, N)] = 0.0;
            }
        }
    }

    indexlist_clear(ws->idx);
    for (int i = 0; i < N; i++) {
        if (!PROB_IS_INT(ws->prob[i], eps)) {
            indexlist_add(ws->idx, i);
        }
    }

    indexlist_shuffle(ws->idx);
}

static void cube_update(CubeWorkspace *ws, int n_cand) {
    double eps = ws->eps;

    int nrows = n_cand - 1;
    int ncols = n_cand;

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            int id = ws->candidates[j];
            ws->bmat[IDX_RM(i, j, ncols)] = ws->amat[IDX_CM(id, i, ws->N)];
        }
    }

    cube_rref(ws->bmat, nrows, ncols);
    cube_null_vector(ws->uvec, ws->bmat, ncols);

    double lambda1 = DBL_MAX;
    double lambda2 = DBL_MAX;

    for (int i = 0; i < n_cand; i++) {
        int id = ws->candidates[i];
        double pi = ws->prob[id];
        double u = ws->uvec[i];

        double lval1 = fabs(pi / u);
        double lval2 = fabs((1.0 - pi) / u);

        if (u >= 0.0) {
            if (lambda1 > lval2) lambda1 = lval2;
            if (lambda2 > lval1) lambda2 = lval1;
        } else {
            if (lambda1 > lval1) lambda1 = lval1;
            if (lambda2 > lval2) lambda2 = lval2;
        }
    }

    double lambda = (unif_rand() * (lambda1 + lambda2) < lambda2) ? lambda1 : -lambda2;

    for (int i = 0; i < n_cand; i++) {
        int id = ws->candidates[i];
        ws->prob[id] += lambda * ws->uvec[i];

        if (ws->prob[id] < eps) ws->prob[id] = 0.0;
        if (ws->prob[id] > 1.0 - eps) ws->prob[id] = 1.0;
    }
}

void cube_flight(CubeWorkspace *ws) {
    int p = ws->p;
    double eps = ws->eps;
    int max_cand = p + 1;

    while (ws->idx->len >= max_cand) {
        for (int i = 0; i < max_cand; i++) {
            ws->candidates[i] = indexlist_get(ws->idx, i);
        }

        cube_update(ws, max_cand);

        for (int i = 0; i < max_cand; i++) {
            int id = ws->candidates[i];
            if (PROB_IS_INT(ws->prob[id], eps)) {
                indexlist_remove(ws->idx, id);
            }
        }
    }
}

void cube_landing(CubeWorkspace *ws) {
    double eps = ws->eps;

    while (ws->idx->len > 1) {
        int n_cand = ws->idx->len;

        for (int i = 0; i < n_cand; i++) {
            ws->candidates[i] = indexlist_get(ws->idx, i);
        }

        cube_update(ws, n_cand);

        for (int i = 0; i < n_cand; i++) {
            int id = ws->candidates[i];
            if (PROB_IS_INT(ws->prob[id], eps)) {
                indexlist_remove(ws->idx, id);
            }
        }
    }

    if (ws->idx->len == 1) {
        int id = indexlist_get(ws->idx, 0);
        if (unif_rand() < ws->prob[id]) {
            ws->prob[id] = 1.0;
        } else {
            ws->prob[id] = 0.0;
        }
        indexlist_remove(ws->idx, id);
    }
}

int *cube_sample(CubeWorkspace *ws, int *sample_size) {
    int N = ws->N;
    double eps = ws->eps;

    int n = 0;
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(ws->prob[i], eps)) n++;
    }

    int *sample = (int *)R_alloc(n, sizeof(int));
    int k = 0;
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(ws->prob[i], eps)) {
            sample[k++] = i + 1;  /* R is 1-indexed */
        }
    }

    *sample_size = n;
    return sample;
}

static void build_amat_stratum(CubeWorkspace *ws, const double *prob,
                                const double *X, const int *subset, int n_subset) {
    int N = ws->N;
    int p = ws->p;  /* This is p+1 (including the 1 column) */
    double eps = ws->eps;

    /* First column: 1 for sample size constraint */
    for (int i = 0; i < n_subset; i++) {
        int id = subset[i];
        ws->amat[IDX_CM(id, 0, N)] = 1.0;
    }

    /* Remaining columns: X/π */
    for (int j = 1; j < p; j++) {
        for (int i = 0; i < n_subset; i++) {
            int id = subset[i];
            double pi = prob[id];
            if (pi > eps) {
                ws->amat[IDX_CM(id, j, N)] = X[IDX_CM(id, j - 1, N)] / pi;
            } else {
                ws->amat[IDX_CM(id, j, N)] = 0.0;
            }
        }
    }
}

static void build_amat_global(CubeWorkspace *ws, const double *prob,
                               const double *X, int p_orig, int n_strata,
                               const int *strata, const int *strata_levels) {
    int N = ws->N;
    double eps = ws->eps;

    for (int i = 0; i < ws->idx->len; i++) {
        int id = indexlist_get(ws->idx, i);

        for (int h = 0; h < n_strata; h++) {
            ws->amat[IDX_CM(id, h, N)] = (strata[id] == strata_levels[h]) ? 1.0 : 0.0;
        }

        double pi = prob[id];
        for (int j = 0; j < p_orig; j++) {
            if (pi > eps) {
                ws->amat[IDX_CM(id, n_strata + j, N)] = X[IDX_CM(id, j, N)] / pi;
            } else {
                ws->amat[IDX_CM(id, n_strata + j, N)] = 0.0;
            }
        }
    }
}

SEXP C_cube_stratified(SEXP prob_sexp, SEXP X_sexp, SEXP strata_sexp, SEXP eps_sexp) {
    int N = length(prob_sexp);
    SEXP dim = getAttrib(X_sexp, R_DimSymbol);
    int p_orig = isNull(dim) ? 1 : INTEGER(dim)[1];
    double eps = asReal(eps_sexp);

    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);
    const int *strata = INTEGER(strata_sexp);

    int max_strata = 0;
    for (int i = 0; i < N; i++) {
        if (strata[i] > max_strata) max_strata = strata[i];
    }
    int n_strata = max_strata;

    int *strata_count = (int *)R_alloc(n_strata, sizeof(int));
    int *strata_levels = (int *)R_alloc(n_strata, sizeof(int));
    memset(strata_count, 0, n_strata * sizeof(int));

    for (int h = 0; h < n_strata; h++) {
        strata_levels[h] = h + 1;
    }

    double *work_prob = (double *)R_alloc(N, sizeof(double));
    memcpy(work_prob, prob, N * sizeof(double));

    int **stratum_units = (int **)R_alloc(n_strata, sizeof(int *));
    int *stratum_sizes = (int *)R_alloc(n_strata, sizeof(int));

    for (int h = 0; h < n_strata; h++) {
        stratum_units[h] = (int *)R_alloc(N, sizeof(int));
        stratum_sizes[h] = 0;
    }

    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(work_prob[i], eps)) {
            work_prob[i] = 1.0;
        } else if (PROB_IS_ZERO(work_prob[i], eps)) {
            work_prob[i] = 0.0;
        } else {
            int h = strata[i] - 1;  /* 0-indexed */
            stratum_units[h][stratum_sizes[h]++] = i;
        }
    }

    GetRNGstate();

    /* PHASE 1: Flight per stratum */
    int p_stratum = p_orig + 1;  /* 1 column + p_orig columns */
    int max_cand_stratum = p_stratum + 1;

    CubeWorkspace *ws = cube_workspace_alloc(N, p_stratum, eps);
    memcpy(ws->prob, work_prob, N * sizeof(double));

    for (int h = 0; h < n_strata; h++) {
        if (stratum_sizes[h] < max_cand_stratum) continue;

        indexlist_clear(ws->idx);
        for (int i = 0; i < stratum_sizes[h]; i++) {
            indexlist_add(ws->idx, stratum_units[h][i]);
        }
        indexlist_shuffle(ws->idx);

        build_amat_stratum(ws, prob, X, stratum_units[h], stratum_sizes[h]);

        cube_flight(ws);

        int new_size = 0;
        for (int i = 0; i < stratum_sizes[h]; i++) {
            int id = stratum_units[h][i];
            if (!PROB_IS_INT(ws->prob[id], eps)) {
                stratum_units[h][new_size++] = id;
            }
        }
        stratum_sizes[h] = new_size;
    }

    memcpy(work_prob, ws->prob, N * sizeof(double));

    /* PHASE 2: Flight on full */
    int total_remaining = 0;
    for (int h = 0; h < n_strata; h++) {
        total_remaining += stratum_sizes[h];
    }

    int active_strata = 0;
    for (int h = 0; h < n_strata; h++) {
        if (stratum_sizes[h] > 0) active_strata++;
    }

    if (total_remaining > 0 && active_strata > 0) {
        int p_global = active_strata + p_orig;
        int max_cand_global = p_global + 1;

        if (total_remaining >= max_cand_global) {
            CubeWorkspace *ws_global = cube_workspace_alloc(N, p_global, eps);
            memcpy(ws_global->prob, work_prob, N * sizeof(double));

            int *active_strata_levels = (int *)R_alloc(active_strata, sizeof(int));
            int ai = 0;
            for (int h = 0; h < n_strata; h++) {
                if (stratum_sizes[h] > 0) {
                    active_strata_levels[ai++] = h + 1;
                }
            }

            indexlist_clear(ws_global->idx);
            for (int h = 0; h < n_strata; h++) {
                for (int i = 0; i < stratum_sizes[h]; i++) {
                    indexlist_add(ws_global->idx, stratum_units[h][i]);
                }
            }
            indexlist_shuffle(ws_global->idx);

            build_amat_global(ws_global, prob, X, p_orig, active_strata,
                             strata, active_strata_levels);

            cube_flight(ws_global);

            memcpy(work_prob, ws_global->prob, N * sizeof(double));
            for (int h = 0; h < n_strata; h++) {
                int new_size = 0;
                for (int i = 0; i < stratum_sizes[h]; i++) {
                    int id = stratum_units[h][i];
                    if (!PROB_IS_INT(ws_global->prob[id], eps)) {
                        stratum_units[h][new_size++] = id;
                    }
                }
                stratum_sizes[h] = new_size;
            }
        }
    }

    /* PHASE 3: Landing per stratum */
    memcpy(ws->prob, work_prob, N * sizeof(double));

    for (int h = 0; h < n_strata; h++) {
        if (stratum_sizes[h] == 0) continue;

        indexlist_clear(ws->idx);
        for (int i = 0; i < stratum_sizes[h]; i++) {
            indexlist_add(ws->idx, stratum_units[h][i]);
        }

        build_amat_stratum(ws, prob, X, stratum_units[h], stratum_sizes[h]);

        cube_landing(ws);
    }

    PutRNGstate();

    /* Extract sample */
    int n_selected = 0;
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(ws->prob[i], eps)) n_selected++;
    }

    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    int *res = INTEGER(result);
    int k = 0;
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(ws->prob[i], eps)) {
            res[k++] = i + 1;  /* 1-indexed */
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP C_cube(SEXP prob_sexp, SEXP X_sexp, SEXP eps_sexp) {
    int N = length(prob_sexp);
    SEXP dim = getAttrib(X_sexp, R_DimSymbol);
    int p = isNull(dim) ? 1 : INTEGER(dim)[1];
    double eps = asReal(eps_sexp);

    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);

    CubeWorkspace *ws = cube_workspace_alloc(N, p, eps);

    GetRNGstate();
    cube_workspace_init(ws, prob, X);
    cube_flight(ws);
    cube_landing(ws);
    PutRNGstate();

    int sample_size;
    int *sample = cube_sample(ws, &sample_size);

    SEXP result = PROTECT(allocVector(INTSXP, sample_size));
    memcpy(INTEGER(result), sample, sample_size * sizeof(int));

    UNPROTECT(1);
    return result;
}

SEXP C_cube_batch(SEXP prob_sexp, SEXP X_sexp, SEXP eps_sexp, SEXP n_samples_sexp) {
    int N = length(prob_sexp);
    SEXP dim = getAttrib(X_sexp, R_DimSymbol);
    int p = isNull(dim) ? 1 : INTEGER(dim)[1];
    double eps = asReal(eps_sexp);
    int nrep = asInteger(n_samples_sexp);

    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);

    double sum_prob = 0.0;
    for (int i = 0; i < N; i++) {
        sum_prob += prob[i];
    }
    int n = (int)(sum_prob + 0.5);

    CubeWorkspace *ws = cube_workspace_alloc(N, p, eps);

    SEXP result = PROTECT(allocMatrix(INTSXP, n, nrep));
    int *result_ptr = INTEGER(result);

    GetRNGstate();

    for (int rep = 0; rep < nrep; rep++) {
        cube_workspace_init(ws, prob, X);

        cube_flight(ws);
        cube_landing(ws);

        int *col = result_ptr + rep * n;
        int k = 0;
        for (int i = 0; i < N && k < n; i++) {
            if (PROB_IS_ONE(ws->prob[i], eps)) {
                col[k++] = i + 1;  /* 1-indexed */
            }
        }

        while (k < n) {
            col[k++] = NA_INTEGER;
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return result;
}
