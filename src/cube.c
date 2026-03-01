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
 * - Chauvet, G. (2009). Stratified balanced sampling. Survey
 *   Methodology, 35, 115-119.
 *
 * NOTE on auxiliary variable ordering:
 *   The landing phase relaxes constraints starting from the LAST variable.
 *   Users must order auxiliary variables by importance (most important first).
 *   For the stratified cube, the within-stratum sample size constraint is
 *   automatically placed in position 0 (always preserved).
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/* Row-major indexing (for RREF working matrix) */
#define IDX_RM(row, col, ncols) ((row) * (ncols) + (col))

/* Column-major indexing (for R matrices) */
#define IDX_CM(row, col, nrows) ((col) * (nrows) + (row))

/* Check if probability is effectively 0 or 1 */
#define PROB_IS_ZERO(p, eps) ((p) < (eps))
#define PROB_IS_ONE(p, eps)  ((p) > 1.0 - (eps))
#define PROB_IS_INT(p, eps)  (PROB_IS_ZERO(p, eps) || PROB_IS_ONE(p, eps))

/* Tolerance for RREF pivoting (relative to eps) */
#define RREF_TOL(eps) ((eps) * 0.01)

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
    double *bmat;       /* Working matrix, row-major [(p+1) x (p+1)] for RREF */
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

static void indexlist_shuffle(IndexList *il) {
    for (int i = 0; i < il->len - 1; i++) {
        int j = i + (int)(unif_rand() * (il->len - i));
        if (j >= il->len) j = il->len - 1;  /* Safety clamp */
        if (i != j) {
            int tmp = il->list[i];
            il->list[i] = il->list[j];
            il->list[j] = tmp;
            il->reverse[il->list[i]] = i;
            il->reverse[il->list[j]] = j;
        }
    }
}

static void cube_rref(double *mat, int nrows, int ncols, double tol) {
    int lead = 0;

    for (int r = 0; r < nrows; r++) {
        if (lead >= ncols)
            return;

        int best = r;
        double best_val = fabs(mat[IDX_RM(r, lead, ncols)]);
        for (int i = r + 1; i < nrows; i++) {
            double v = fabs(mat[IDX_RM(i, lead, ncols)]);
            if (v > best_val) {
                best = i;
                best_val = v;
            }
        }

        if (best_val < tol) {
            /* No valid pivot in this column, skip it */
            r--;
            lead++;
            continue;
        }

        if (best != r) {
            for (int k = 0; k < ncols; k++) {
                double tmp = mat[IDX_RM(best, k, ncols)];
                mat[IDX_RM(best, k, ncols)] = mat[IDX_RM(r, k, ncols)];
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
            if (fabs(factor) < tol) {
                mat[IDX_RM(j, lead, ncols)] = 0.0;
                continue;
            }
            mat[IDX_RM(j, lead, ncols)] = 0.0;
            for (int k = lead + 1; k < ncols; k++) {
                mat[IDX_RM(j, k, ncols)] -= mat[IDX_RM(r, k, ncols)] * factor;
            }
        }

        lead++;
    }
}

static void cube_null_vector(double *uvec, const double *mat, int ncols,
                             double tol) {
    int nrows = ncols - 1;

    if (nrows <= 0) {
        uvec[0] = 1.0;
        return;
    }

    if (fabs(mat[IDX_RM(nrows - 1, nrows - 1, ncols)] - 1.0) < tol) {
        uvec[ncols - 1] = 1.0;
        for (int i = 0; i < nrows; i++) {
            uvec[i] = -mat[IDX_RM(i, ncols - 1, ncols)];
        }
        return;
    }

    for (int k = 0; k < ncols; k++) {
        uvec[k] = (k % 2 == 0) ? -1.0 : 1.0;
    }

    for (int i = 0; i < nrows; i++) {
        int lead = 0;

        for (; lead < ncols; lead++) {
            if (fabs(mat[IDX_RM(i, lead, ncols)] - 1.0) < tol) break;
        }

        if (lead >= ncols) continue;

        uvec[lead] = 0.0;
        for (int k = lead + 1; k < ncols; k++) {
            uvec[lead] -= uvec[k] * mat[IDX_RM(i, k, ncols)];
        }
    }
}

static CubeWorkspace *cube_workspace_alloc(int N, int p, double eps) {
    CubeWorkspace *ws = (CubeWorkspace *)R_alloc(1, sizeof(CubeWorkspace));

    ws->N = N;
    ws->p = p;
    ws->eps = eps;

    ws->prob = (double *)R_alloc(N, sizeof(double));
    ws->amat = (double *)R_alloc((size_t)N * p, sizeof(double));
    ws->bmat = (double *)R_alloc((size_t)(p + 1) * (p + 1), sizeof(double));
    ws->uvec = (double *)R_alloc(p + 1, sizeof(double));
    ws->candidates = (int *)R_alloc(p + 1, sizeof(int));

    ws->idx = indexlist_alloc(N);

    return ws;
}

static void cube_workspace_build_amat(CubeWorkspace *ws, const double *prob,
                                      const double *X) {
    int N = ws->N;
    int p = ws->p;
    double eps = ws->eps;

    for (int j = 0; j < p; j++) {
        for (int i = 0; i < N; i++) {
            double pi = prob[i];
            if (pi > eps) {
                ws->amat[IDX_CM(i, j, N)] = X[IDX_CM(i, j, N)] / pi;
            } else {
                ws->amat[IDX_CM(i, j, N)] = 0.0;
            }
        }
    }
}

static void cube_workspace_reset(CubeWorkspace *ws, const double *prob) {
    int N = ws->N;
    double eps = ws->eps;

    memcpy(ws->prob, prob, N * sizeof(double));

    indexlist_clear(ws->idx);
    for (int i = 0; i < N; i++) {
        if (!PROB_IS_INT(ws->prob[i], eps)) {
            indexlist_add(ws->idx, i);
        }
    }

    indexlist_shuffle(ws->idx);
}

static void cube_workspace_init(CubeWorkspace *ws, const double *prob,
                                const double *X) {
    cube_workspace_build_amat(ws, prob, X);
    cube_workspace_reset(ws, prob);
}

static void cube_update(CubeWorkspace *ws, int n_cand) {
    double eps = ws->eps;
    double tol = RREF_TOL(eps);

    int nrows = n_cand - 1;
    int ncols = n_cand;

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            int id = ws->candidates[j];
            ws->bmat[IDX_RM(i, j, ncols)] = ws->amat[IDX_CM(id, i, ws->N)];
        }
    }

    cube_rref(ws->bmat, nrows, ncols, tol);
    cube_null_vector(ws->uvec, ws->bmat, ncols, tol);

    double lambda1 = DBL_MAX;
    double lambda2 = DBL_MAX;

    for (int i = 0; i < n_cand; i++) {
        int id = ws->candidates[i];
        double pi = ws->prob[id];
        double u = ws->uvec[i];

        if (fabs(u) < tol) continue;

        double lval1 = fabs(pi / u);
        double lval2 = fabs((1.0 - pi) / u);

        if (u > 0.0) {
            if (lambda1 > lval2) lambda1 = lval2;
            if (lambda2 > lval1) lambda2 = lval1;
        } else {
            if (lambda1 > lval1) lambda1 = lval1;
            if (lambda2 > lval2) lambda2 = lval2;
        }
    }

    if (lambda1 + lambda2 < tol) return;

    double lambda = (unif_rand() * (lambda1 + lambda2) < lambda2)
                    ? lambda1 : -lambda2;

    for (int i = 0; i < n_cand; i++) {
        int id = ws->candidates[i];
        ws->prob[id] += lambda * ws->uvec[i];

        if (ws->prob[id] < eps) ws->prob[id] = 0.0;
        if (ws->prob[id] > 1.0 - eps) ws->prob[id] = 1.0;
    }
}

static void cube_flight(CubeWorkspace *ws) {
    int max_cand = ws->p + 1;
    double eps = ws->eps;
    int max_iter = (ws->N > 46340) ? INT_MAX : ws->N * ws->N;  /* Safety limit */
    int iter = 0;

    while (ws->idx->len >= max_cand) {
        R_CheckUserInterrupt();
        if (++iter > max_iter) {
            error("cube flight phase did not converge after %d iterations. "
                  "This usually indicates collinear balancing variables.",
                  max_iter);
        }

        for (int i = 0; i < max_cand; i++) {
            ws->candidates[i] = ws->idx->list[i];
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

static void cube_landing(CubeWorkspace *ws) {
    double eps = ws->eps;
    int max_iter = (ws->N > 46340) ? INT_MAX : ws->N * ws->N;  /* Safety limit */
    int iter = 0;

    while (ws->idx->len > 1) {
        R_CheckUserInterrupt();
        if (++iter > max_iter) {
            error("cube landing phase did not converge after %d iterations.",
                  max_iter);
        }

        int n_cand = ws->idx->len;

        for (int i = 0; i < n_cand; i++) {
            ws->candidates[i] = ws->idx->list[i];
        }

        cube_update(ws, n_cand);

        for (int i = 0; i < n_cand; i++) {
            int id = ws->candidates[i];
            if (PROB_IS_INT(ws->prob[id], eps)) {
                indexlist_remove(ws->idx, id);
            }
        }
    }

    /* Single remaining unit: resolve by coin flip */
    if (ws->idx->len == 1) {
        int id = ws->idx->list[0];
        ws->prob[id] = (unif_rand() < ws->prob[id]) ? 1.0 : 0.0;
        indexlist_remove(ws->idx, id);
    }
}

static void build_amat_stratum(CubeWorkspace *ws, const double *prob,
                               const double *X, const int *subset,
                               int n_subset) {
    int N = ws->N;
    int p = ws->p;
    double eps = ws->eps;

    for (int i = 0; i < n_subset; i++) {
        int id = subset[i];
        ws->amat[IDX_CM(id, 0, N)] = 1.0;
    }

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

/* Build A matrix for global flight (phase 2 of stratified cube):
 * Columns 0..H-1: stratum indicator dummies (per-stratum size constraints)
 * Columns H..H+p-1: X[,j] / pi (balancing variables) */
static void build_amat_global(CubeWorkspace *ws, const double *prob,
                              const double *X, int p_orig, int n_strata,
                              const int *strata, const int *strata_levels) {
    int N = ws->N;
    double eps = ws->eps;

    for (int i = 0; i < ws->idx->len; i++) {
        int id = ws->idx->list[i];

        for (int h = 0; h < n_strata; h++) {
            ws->amat[IDX_CM(id, h, N)] =
                (strata[id] == strata_levels[h]) ? 1.0 : 0.0;
        }

        double pi = prob[id];
        for (int j = 0; j < p_orig; j++) {
            if (pi > eps) {
                ws->amat[IDX_CM(id, n_strata + j, N)] =
                    X[IDX_CM(id, j, N)] / pi;
            } else {
                ws->amat[IDX_CM(id, n_strata + j, N)] = 0.0;
            }
        }
    }
}

/* Precompute X/pi for original auxiliary variables (N x p_orig, column-major).
 * Used by stratified batch to avoid repeated division across replicates. */
static void build_x_over_pi(double *Xpi, const double *prob, const double *X,
                            int N, int p_orig, double eps) {
    for (int j = 0; j < p_orig; j++) {
        for (int i = 0; i < N; i++) {
            double pi = prob[i];
            if (pi > eps) {
                Xpi[IDX_CM(i, j, N)] = X[IDX_CM(i, j, N)] / pi;
            } else {
                Xpi[IDX_CM(i, j, N)] = 0.0;
            }
        }
    }
}

/* Like build_amat_stratum, but reads from precomputed Xpi instead of dividing. */
static void build_amat_stratum_xpi(CubeWorkspace *ws, const double *Xpi,
                                   const int *subset, int n_subset) {
    int N = ws->N;
    int p = ws->p;

    for (int i = 0; i < n_subset; i++) {
        int id = subset[i];
        ws->amat[IDX_CM(id, 0, N)] = 1.0;
    }

    for (int j = 1; j < p; j++) {
        for (int i = 0; i < n_subset; i++) {
            int id = subset[i];
            ws->amat[IDX_CM(id, j, N)] = Xpi[IDX_CM(id, j - 1, N)];
        }
    }
}

/* Like build_amat_global, but reads from precomputed Xpi instead of dividing. */
static void build_amat_global_xpi(CubeWorkspace *ws, const double *Xpi,
                                  int p_orig, int n_strata,
                                  const int *strata, const int *strata_levels) {
    int N = ws->N;

    for (int i = 0; i < ws->idx->len; i++) {
        int id = ws->idx->list[i];

        for (int h = 0; h < n_strata; h++) {
            ws->amat[IDX_CM(id, h, N)] =
                (strata[id] == strata_levels[h]) ? 1.0 : 0.0;
        }

        for (int j = 0; j < p_orig; j++) {
            ws->amat[IDX_CM(id, n_strata + j, N)] = Xpi[IDX_CM(id, j, N)];
        }
    }
}

/* --- Shared utilities --------------------------------------------------- */

static int cube_extract_selected(const double *prob, int N, double eps,
                                 int *out_idx) {
    int k = 0;
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(prob[i], eps)) {
            out_idx[k++] = i + 1;
        }
    }
    return k;
}

static int cube_strata_nlevels(const int *strata, int N) {
    int n_strata = 0;
    for (int i = 0; i < N; i++) {
        if (strata[i] < 1) {
            error("strata values must be positive integers, got %d", strata[i]);
        }
        if (strata[i] > n_strata) n_strata = strata[i];
    }
    return n_strata;
}

static void cube_check_X_dim(SEXP X_sexp, int N, int *p_out) {
    SEXP dim = getAttrib(X_sexp, R_DimSymbol);
    int p = isNull(dim) ? 1 : INTEGER(dim)[1];
    int X_nrow = isNull(dim) ? length(X_sexp) : INTEGER(dim)[0];

    if (X_nrow != N) {
        error("nrow(X) = %d does not match length(prob) = %d", X_nrow, N);
    }
    *p_out = p;
}

/* --- Internal draw helpers (no RNG state management) -------------------- */

static int cube_draw_unstratified(const double *prob, const double *X,
                                  int N, int p, double eps, int *out_idx) {
    CubeWorkspace *ws = cube_workspace_alloc(N, p, eps);
    cube_workspace_init(ws, prob, X);
    cube_flight(ws);
    cube_landing(ws);
    return cube_extract_selected(ws->prob, N, eps, out_idx);
}

static int cube_draw_stratified(const double *prob, const double *X,
                                const int *strata, int N, int p_orig,
                                int n_strata, double eps, int *out_idx) {
    int *strata_levels = (int *)R_alloc(n_strata, sizeof(int));
    for (int h = 0; h < n_strata; h++) {
        strata_levels[h] = h + 1;
    }

    double *work_prob = (double *)R_alloc(N, sizeof(double));
    memcpy(work_prob, prob, N * sizeof(double));

    int *stratum_sizes = (int *)R_alloc(n_strata, sizeof(int));
    memset(stratum_sizes, 0, n_strata * sizeof(int));

    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(work_prob[i], eps)) {
            work_prob[i] = 1.0;
        } else if (PROB_IS_ZERO(work_prob[i], eps)) {
            work_prob[i] = 0.0;
        } else {
            stratum_sizes[strata[i] - 1]++;
        }
    }

    int *offsets = (int *)R_alloc(n_strata, sizeof(int));
    offsets[0] = 0;
    for (int h = 1; h < n_strata; h++) {
        offsets[h] = offsets[h - 1] + stratum_sizes[h - 1];
    }
    int total_active = offsets[n_strata - 1] + stratum_sizes[n_strata - 1];

    int *all_units = (int *)R_alloc(total_active > 0 ? total_active : 1,
                                    sizeof(int));
    int **stratum_units = (int **)R_alloc(n_strata, sizeof(int *));
    for (int h = 0; h < n_strata; h++) {
        stratum_units[h] = all_units + offsets[h];
    }

    int *fill_pos = (int *)R_alloc(n_strata, sizeof(int));
    memset(fill_pos, 0, n_strata * sizeof(int));
    for (int i = 0; i < N; i++) {
        if (!PROB_IS_INT(work_prob[i], eps)) {
            int h = strata[i] - 1;
            stratum_units[h][fill_pos[h]++] = i;
        }
    }

    /* PHASE 1: Flight per stratum */
    int p_stratum = p_orig + 1;  /* 1 (size constraint) + p_orig */
    int max_cand_stratum = p_stratum + 1;

    CubeWorkspace *ws = cube_workspace_alloc(N, p_stratum, eps);
    memcpy(ws->prob, work_prob, N * sizeof(double));

    for (int h = 0; h < n_strata; h++) {
        R_CheckUserInterrupt();
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

    /* PHASE 2: Global flight across strata */
    int total_remaining = 0;
    int active_strata = 0;
    for (int h = 0; h < n_strata; h++) {
        total_remaining += stratum_sizes[h];
        if (stratum_sizes[h] > 0) active_strata++;
    }

    if (total_remaining > 0 && active_strata > 0) {
        int p_global = active_strata + p_orig;
        int max_cand_global = p_global + 1;

        if (total_remaining >= max_cand_global) {
            CubeWorkspace *ws_global =
                cube_workspace_alloc(N, p_global, eps);
            memcpy(ws_global->prob, work_prob, N * sizeof(double));

            int *active_levels = (int *)R_alloc(active_strata, sizeof(int));
            int ai = 0;
            for (int h = 0; h < n_strata; h++) {
                if (stratum_sizes[h] > 0) {
                    active_levels[ai++] = h + 1;
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
                              strata, active_levels);

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
        R_CheckUserInterrupt();
        if (stratum_sizes[h] == 0) continue;

        indexlist_clear(ws->idx);
        for (int i = 0; i < stratum_sizes[h]; i++) {
            indexlist_add(ws->idx, stratum_units[h][i]);
        }

        build_amat_stratum(ws, prob, X, stratum_units[h], stratum_sizes[h]);

        cube_landing(ws);
    }

    return cube_extract_selected(ws->prob, N, eps, out_idx);
}

/* --- R-callable entry points -------------------------------------------- */

SEXP C_cube(SEXP prob_sexp, SEXP X_sexp, SEXP eps_sexp) {
    const int N = length(prob_sexp);
    int p = 0;
    cube_check_X_dim(X_sexp, N, &p);

    const double eps = asReal(eps_sexp);
    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);
    int *selected = (int *)R_alloc(N, sizeof(int));

    GetRNGstate();
    int n_selected = cube_draw_unstratified(prob, X, N, p, eps, selected);
    PutRNGstate();

    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
    UNPROTECT(1);
    return result;
}

SEXP C_cube_batch(SEXP prob_sexp, SEXP X_sexp, SEXP eps_sexp,
                  SEXP nrep_sexp) {
    const int N = length(prob_sexp);
    int p = 0;
    cube_check_X_dim(X_sexp, N, &p);

    const int nrep = INTEGER(nrep_sexp)[0];
    const double eps = asReal(eps_sexp);
    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);

    double n_sum = 0.0;
    for (int i = 0; i < N; i++) n_sum += prob[i];
    const int n_target = (int)(n_sum + 0.5);

    SEXP result = PROTECT(allocMatrix(INTSXP, n_target, nrep));
    int *res = INTEGER(result);

    if (nrep <= 0) {
        UNPROTECT(1);
        return result;
    }

    /* Allocate workspace once; build A = X/pi once (constant across reps) */
    CubeWorkspace *ws = cube_workspace_alloc(N, p, eps);
    cube_workspace_build_amat(ws, prob, X);

    GetRNGstate();
    for (int s = 0; s < nrep; s++) {
        if (s % 32 == 0) R_CheckUserInterrupt();

        cube_workspace_reset(ws, prob);
        cube_flight(ws);
        cube_landing(ws);

        int *col = res + s * n_target;
        int n_selected = cube_extract_selected(ws->prob, N, eps, col);
        if (n_selected != n_target) {
            PutRNGstate();
            UNPROTECT(1);
            error(
                "cube batch draw %d produced size %d, expected %d",
                s + 1,
                n_selected,
                n_target
            );
        }
    }
    PutRNGstate();

    UNPROTECT(1);
    return result;
}

/*
 * C_cube_stratified: Stratified cube method (Chauvet & Tille, 2006/2009)
 *
 * Three phases:
 *   1. Per-stratum flight: balance within each stratum
 *   2. Global flight: balance across strata with stratum size constraints
 *   3. Per-stratum landing: resolve remaining undecided units
 */
SEXP C_cube_stratified(SEXP prob_sexp, SEXP X_sexp, SEXP strata_sexp,
                       SEXP eps_sexp) {
    const int N = length(prob_sexp);
    int p_orig = 0;
    cube_check_X_dim(X_sexp, N, &p_orig);

    if (length(strata_sexp) != N) {
        error("length(strata) = %d does not match length(prob) = %d",
              length(strata_sexp), N);
    }

    const double eps = asReal(eps_sexp);
    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);
    const int *strata = INTEGER(strata_sexp);
    const int n_strata = cube_strata_nlevels(strata, N);
    int *selected = (int *)R_alloc(N, sizeof(int));

    GetRNGstate();
    int n_selected = cube_draw_stratified(
        prob, X, strata, N, p_orig, n_strata, eps, selected
    );
    PutRNGstate();

    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
    UNPROTECT(1);
    return result;
}

SEXP C_cube_stratified_batch(SEXP prob_sexp, SEXP X_sexp, SEXP strata_sexp,
                             SEXP eps_sexp, SEXP nrep_sexp) {
    const int N = length(prob_sexp);
    int p_orig = 0;
    cube_check_X_dim(X_sexp, N, &p_orig);

    if (length(strata_sexp) != N) {
        error("length(strata) = %d does not match length(prob) = %d",
              length(strata_sexp), N);
    }

    const int nrep = INTEGER(nrep_sexp)[0];
    const double eps = asReal(eps_sexp);
    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);
    const int *strata = INTEGER(strata_sexp);
    const int n_strata = cube_strata_nlevels(strata, N);

    double n_sum = 0.0;
    for (int i = 0; i < N; i++) n_sum += prob[i];
    const int n_target = (int)(n_sum + 0.5);

    SEXP result = PROTECT(allocMatrix(INTSXP, n_target, nrep));
    int *res = INTEGER(result);

    if (nrep <= 0) {
        UNPROTECT(1);
        return result;
    }

    /* Precompute base state (constant across replicates) */
    double *base_prob = (double *)R_alloc(N, sizeof(double));
    memcpy(base_prob, prob, N * sizeof(double));

    int *base_sizes = (int *)R_alloc(n_strata, sizeof(int));
    memset(base_sizes, 0, n_strata * sizeof(int));
    for (int i = 0; i < N; i++) {
        if (PROB_IS_ONE(base_prob[i], eps)) {
            base_prob[i] = 1.0;
        } else if (PROB_IS_ZERO(base_prob[i], eps)) {
            base_prob[i] = 0.0;
        } else {
            base_sizes[strata[i] - 1]++;
        }
    }

    int *offsets = (int *)R_alloc(n_strata, sizeof(int));
    offsets[0] = 0;
    for (int h = 1; h < n_strata; h++) {
        offsets[h] = offsets[h - 1] + base_sizes[h - 1];
    }
    int total_active = offsets[n_strata - 1] + base_sizes[n_strata - 1];

    int *base_all_units = (int *)R_alloc(total_active > 0 ? total_active : 1,
                                         sizeof(int));
    int **base_units = (int **)R_alloc(n_strata, sizeof(int *));
    for (int h = 0; h < n_strata; h++) {
        base_units[h] = base_all_units + offsets[h];
    }

    int *fill_pos = (int *)R_alloc(n_strata, sizeof(int));
    memset(fill_pos, 0, n_strata * sizeof(int));
    for (int i = 0; i < N; i++) {
        if (!PROB_IS_INT(base_prob[i], eps)) {
            int h = strata[i] - 1;
            base_units[h][fill_pos[h]++] = i;
        }
    }

    /* Working copies restored from base each replicate */
    int *stratum_sizes = (int *)R_alloc(n_strata, sizeof(int));
    int *work_all_units = (int *)R_alloc(total_active > 0 ? total_active : 1,
                                         sizeof(int));
    int **stratum_units = (int **)R_alloc(n_strata, sizeof(int *));
    for (int h = 0; h < n_strata; h++) {
        stratum_units[h] = work_all_units + offsets[h];
    }
    int *active_levels = (int *)R_alloc(n_strata, sizeof(int));

    /* Precompute X/pi once (avoids N*p_orig divisions per replicate) */
    double *Xpi = (double *)R_alloc((size_t)N * p_orig, sizeof(double));
    build_x_over_pi(Xpi, prob, X, N, p_orig, eps);

    /* Allocate workspaces once (reused across replicates) */
    int p_stratum = p_orig + 1;
    int max_cand_stratum = p_stratum + 1;
    CubeWorkspace *ws = cube_workspace_alloc(N, p_stratum, eps);
    /* ws_global allocated with max p = n_strata + p_orig; actual p set per rep */
    CubeWorkspace *ws_global = cube_workspace_alloc(N, n_strata + p_orig, eps);

    GetRNGstate();
    for (int s = 0; s < nrep; s++) {
        if (s % 32 == 0) R_CheckUserInterrupt();

        /* Restore base state */
        memcpy(ws->prob, base_prob, N * sizeof(double));
        memcpy(stratum_sizes, base_sizes, n_strata * sizeof(int));
        if (total_active > 0) {
            memcpy(work_all_units, base_all_units, total_active * sizeof(int));
        }

        /* PHASE 1: Flight per stratum */
        for (int h = 0; h < n_strata; h++) {
            if (stratum_sizes[h] < max_cand_stratum) continue;

            indexlist_clear(ws->idx);
            for (int i = 0; i < stratum_sizes[h]; i++) {
                indexlist_add(ws->idx, stratum_units[h][i]);
            }
            indexlist_shuffle(ws->idx);

            build_amat_stratum_xpi(ws, Xpi, stratum_units[h], stratum_sizes[h]);
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

        /* PHASE 2: Global flight across strata */
        int total_remaining = 0;
        int active_strata = 0;
        for (int h = 0; h < n_strata; h++) {
            total_remaining += stratum_sizes[h];
            if (stratum_sizes[h] > 0) active_strata++;
        }

        if (total_remaining > 0 && active_strata > 0) {
            int p_global = active_strata + p_orig;
            int max_cand_global = p_global + 1;

            if (total_remaining >= max_cand_global) {
                ws_global->p = p_global;
                memcpy(ws_global->prob, ws->prob, N * sizeof(double));

                int ai = 0;
                for (int h = 0; h < n_strata; h++) {
                    if (stratum_sizes[h] > 0) {
                        active_levels[ai++] = h + 1;
                    }
                }

                indexlist_clear(ws_global->idx);
                for (int h = 0; h < n_strata; h++) {
                    for (int i = 0; i < stratum_sizes[h]; i++) {
                        indexlist_add(ws_global->idx, stratum_units[h][i]);
                    }
                }
                indexlist_shuffle(ws_global->idx);

                build_amat_global_xpi(
                    ws_global, Xpi, p_orig, active_strata,
                    strata, active_levels
                );
                cube_flight(ws_global);

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

                memcpy(ws->prob, ws_global->prob, N * sizeof(double));
            }
        }

        /* PHASE 3: Landing per stratum */
        for (int h = 0; h < n_strata; h++) {
            if (stratum_sizes[h] == 0) continue;

            indexlist_clear(ws->idx);
            for (int i = 0; i < stratum_sizes[h]; i++) {
                indexlist_add(ws->idx, stratum_units[h][i]);
            }

            build_amat_stratum_xpi(ws, Xpi, stratum_units[h], stratum_sizes[h]);
            cube_landing(ws);
        }

        int *col = res + s * n_target;
        int n_selected = cube_extract_selected(ws->prob, N, eps, col);
        if (n_selected != n_target) {
            PutRNGstate();
            UNPROTECT(1);
            error(
                "stratified cube batch draw %d produced size %d, expected %d. "
                "Use list-based path for non-integer per-stratum sizes.",
                s + 1,
                n_selected,
                n_target
            );
        }
    }
    PutRNGstate();

    UNPROTECT(1);
    return result;
}
