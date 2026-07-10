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
 * - Tripet, A. and Tillé, Y. (2026). Balanced sampling with
 *   inequalities. JASA, 121(553), 796-806. (Flight phase with
 *   inequality constraints B's <= r: step sizes are capped so the
 *   walk stays inside the polytope, and a constraint whose slack
 *   reaches 0 becomes an equality from then on.)
 *
 * NOTE on auxiliary variable ordering:
 *   The landing phase relaxes constraints starting from the LAST variable.
 *   Users must order auxiliary variables by importance (most important first).
 *   For the stratified cube, the within-stratum sample size constraint is
 *   automatically placed in position 0 (always preserved).
 *   With inequality bounds, the landing priority is: size constraint,
 *   then bounds that became equalities during flight, then auxiliary
 *   variables (still relaxed from the last).
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

/* Inequality row states (cube with bounds) */
#define ROW_INACTIVE 0  /* strict inequality: slack > 0, caps the step */
#define ROW_ACTIVE   1  /* slack hit 0: enforced as an equality */
#define ROW_DROPPED  2  /* relaxed during landing (infeasible to keep) */

/* A row is "on its face" (activatable) when its slack is this small.
 * Scaled by the bound magnitude so large-count bounds behave like
 * small ones. */
#define ROW_ACT_TOL(r) (1e-9 * (1.0 + fabs(r)))

typedef struct {
    int *list;      /* Active indices: list[0..len-1] are valid */
    int *reverse;   /* reverse[id] = position in list, or capacity if removed */
    int len;        /* Current number of active indices */
    int capacity;   /* Maximum capacity (= N) */
} IndexList;

typedef struct {
    int N;              /* Population size */
    int p;              /* Number of balancing variables */
    int q;              /* Number of inequality rows (0 = classic cube) */
    double eps;         /* Numerical tolerance */
    double *prob;       /* Working probabilities [N] */
    double *amat;       /* A = X / pi, column-major [N x p] */
    double *bmat;       /* Working matrix, row-major, for RREF */
    double *uvec;       /* Null space vector [p+q+1] */
    int *candidates;    /* Current candidate indices [p+q+1] */
    IndexList *idx;     /* Active index list */

    /* Inequality constraints B's <= r (unused when q == 0).
     * Constraint ids in `cons`: id < p refers to amat column id,
     * id >= p refers to bcols column id - p. */
    const double *bcols;    /* B, column-major [N x q] (borrowed) */
    double *rvec;           /* Bounds r [q] */
    double *slack;          /* r_j - B_j' prob [q] */
    double *dvec;           /* Work: B_j' u over the candidate window [q] */
    int *row_state;         /* ROW_INACTIVE / ROW_ACTIVE / ROW_DROPPED [q] */
    int *in_kernel;         /* Work: row currently enforced via kernel [q] */
    int *cons;              /* Equality constraint ids [p + q] */
    int *landing_cons;      /* Landing priority order [p + q] */
    int n_eq;               /* Current number of equality constraints */
    int n_dropped;          /* Rows dropped during landing */
    int bind1;              /* Row capping lambda1 in last update, or -1 */
    int bind2;              /* Row capping lambda2 in last update, or -1 */
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

/* Returns the rank (number of pivots placed). */
static int cube_rref(double *mat, int nrows, int ncols, double tol) {
    int lead = 0;
    int rank = 0;

    for (int r = 0; r < nrows; r++) {
        if (lead >= ncols)
            return rank;

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
        rank++;
    }

    return rank;
}

static void cube_null_vector(double *uvec, const double *mat, int nrows,
                             int ncols, double tol) {
    if (nrows <= 0) {
        for (int k = 0; k < ncols; k++) {
            uvec[k] = 0.0;
        }
        uvec[0] = 1.0;
        return;
    }

    if (nrows == ncols - 1 &&
        fabs(mat[IDX_RM(nrows - 1, nrows - 1, ncols)] - 1.0) < tol) {
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

static CubeWorkspace *cube_workspace_alloc(int N, int p, int q, double eps) {
    CubeWorkspace *ws = (CubeWorkspace *)R_alloc(1, sizeof(CubeWorkspace));
    int pq = p + q;

    ws->N = N;
    ws->p = p;
    ws->q = q;
    ws->eps = eps;

    ws->prob = (double *)R_alloc(N, sizeof(double));
    ws->amat = (double *)R_alloc((size_t)N * p, sizeof(double));
    ws->bmat = (double *)R_alloc((size_t)(pq + 1) * (pq + 1), sizeof(double));
    ws->uvec = (double *)R_alloc(pq + 1, sizeof(double));
    ws->candidates = (int *)R_alloc(pq + 1, sizeof(int));

    ws->idx = indexlist_alloc(N);

    ws->bcols = NULL;
    ws->rvec = (double *)R_alloc(q > 0 ? q : 1, sizeof(double));
    ws->slack = (double *)R_alloc(q > 0 ? q : 1, sizeof(double));
    ws->dvec = (double *)R_alloc(q > 0 ? q : 1, sizeof(double));
    ws->row_state = (int *)R_alloc(q > 0 ? q : 1, sizeof(int));
    ws->in_kernel = (int *)R_alloc(q > 0 ? q : 1, sizeof(int));
    ws->cons = (int *)R_alloc(pq > 0 ? pq : 1, sizeof(int));
    ws->landing_cons = (int *)R_alloc(pq > 0 ? pq : 1, sizeof(int));
    /* The first p entries are the original equality constraints and are
     * never reordered; activations append after them. */
    for (int j = 0; j < pq; j++) {
        ws->cons[j] = j;
    }
    ws->n_eq = p;
    ws->n_dropped = 0;
    ws->bind1 = -1;
    ws->bind2 = -1;

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

/* --- Inequality constraints (Tripet & Tillé 2026) --------------------- */

/* Turn inequality row j into an equality constraint: the walk has
 * reached the face B_j'prob = r_j and must stay on it. */
static void cube_row_activate(CubeWorkspace *ws, int j) {
    ws->row_state[j] = ROW_ACTIVE;
    ws->slack[j] = 0.0;
    ws->cons[ws->n_eq++] = ws->p + j;
}

/* Recompute slacks from the current working probabilities and reset
 * the constraint bookkeeping for a fresh draw. Rows starting on their
 * face (slack ~ 0, e.g. equality bounds lower == upper) are activated
 * immediately. Errors if the starting probabilities are infeasible. */
static void cube_ineq_reset(CubeWorkspace *ws) {
    ws->n_eq = ws->p;
    ws->n_dropped = 0;
    ws->bind1 = -1;
    ws->bind2 = -1;

    for (int j = 0; j < ws->q; j++) {
        double dot = 0.0;
        for (int i = 0; i < ws->N; i++) {
            dot += ws->bcols[IDX_CM(i, j, ws->N)] * ws->prob[i];
        }
        double g = ws->rvec[j] - dot;
        if (g < -1e-7 * (1.0 + fabs(ws->rvec[j]))) {
            error(
                "bounds are infeasible at the initial inclusion "
                "probabilities (constraint %d: value %g exceeds bound %g)",
                j + 1, dot, ws->rvec[j]
            );
        }
        ws->slack[j] = g;
        ws->row_state[j] = ROW_INACTIVE;
        if (g <= ROW_ACT_TOL(ws->rvec[j])) {
            cube_row_activate(ws, j);
        }
    }
}

/* One basic step of the flight/landing random walk over the candidate
 * window. The kernel direction is computed from the ncons constraint
 * ids in cons_ids (id < p: amat column; id >= p: inequality row
 * id - p). ncons may exceed n_cand - 1 during landing: dependent
 * constraints (e.g. redundant margin systems) then cost nothing, and
 * only a genuinely full-rank system reports "no kernel".
 *
 * With inequality rows (q > 0), the step sizes are additionally capped
 * so B'prob <= r stays feasible, and any row whose slack reaches 0 is
 * converted to an equality (Tripet & Tillé 2026, Algorithm 2).
 *
 * Returns 1 when progress was made (a move and/or an activation),
 * 0 on a stall (degenerate direction, or a binding row that cannot be
 * activated — the landing phase then relaxes it), and -1 when the
 * constraint system has full column rank (no direction exists; the
 * caller must relax a constraint). No RNG is consumed unless a move
 * is made. */
static int cube_update(CubeWorkspace *ws, int n_cand,
                       const int *cons_ids, int ncons) {
    double eps = ws->eps;
    double tol = RREF_TOL(eps);

    int nrows = ncons;
    int ncols = n_cand;

    ws->bind1 = -1;
    ws->bind2 = -1;

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            int id = ws->candidates[j];
            int c = cons_ids[i];
            ws->bmat[IDX_RM(i, j, ncols)] = (c < ws->p)
                ? ws->amat[IDX_CM(id, c, ws->N)]
                : ws->bcols[IDX_CM(id, c - ws->p, ws->N)];
        }
    }

    int rank = cube_rref(ws->bmat, nrows, ncols, tol);
    if (rank >= ncols) {
        return -1;
    }
    cube_null_vector(ws->uvec, ws->bmat, nrows, ncols, tol);

    if (ws->q > 0) {
        /* Max-norm normalization so the inequality caps and the stall
         * test below are scale-free. Skipped for q == 0 to keep the
         * classic path bit-identical. */
        double umax = 0.0;
        for (int i = 0; i < n_cand; i++) {
            double a = fabs(ws->uvec[i]);
            if (a > umax) umax = a;
        }
        if (umax > 0.0) {
            for (int i = 0; i < n_cand; i++) {
                ws->uvec[i] /= umax;
            }
        }
    }

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

    /*
     * Two-tier guard:
     * (1) If no candidate had a meaningful null-vector entry (fabs(u) >=
     *     tol), lambda1 and lambda2 stayed at their DBL_MAX initializers.
     *     Summing would overflow to +Inf and defeat the < tol test, so
     *     check the sentinel explicitly.
     * (2) Otherwise, ensure the feasible movement range is non-degenerate.
     */
    if (lambda1 == DBL_MAX || lambda2 == DBL_MAX) return 0;

    if (ws->q > 0) {
        /* Which inequality rows are enforced through the kernel in this
         * step? Those move orthogonally (B_j'u = 0) and must not cap. */
        for (int j = 0; j < ws->q; j++) {
            ws->in_kernel[j] = 0;
        }
        for (int i = 0; i < ncons; i++) {
            if (cons_ids[i] >= ws->p) {
                ws->in_kernel[cons_ids[i] - ws->p] = 1;
            }
        }

        for (int j = 0; j < ws->q; j++) {
            ws->dvec[j] = 0.0;
            if (ws->row_state[j] == ROW_DROPPED || ws->in_kernel[j]) {
                continue;
            }
            double d = 0.0;
            for (int i = 0; i < n_cand; i++) {
                d += ws->bcols[IDX_CM(ws->candidates[i], j, ws->N)] *
                     ws->uvec[i];
            }
            ws->dvec[j] = d;
            if (d > tol) {
                double cap = ws->slack[j] / d;
                if (cap < lambda1) {
                    lambda1 = cap;
                    ws->bind1 = j;
                }
            } else if (d < -tol) {
                double cap = ws->slack[j] / (-d);
                if (cap < lambda2) {
                    lambda2 = cap;
                    ws->bind2 = j;
                }
            }
        }
        if (lambda1 < 0.0) lambda1 = 0.0;
        if (lambda2 < 0.0) lambda2 = 0.0;
    }

    /* Stall test. With bounds, either side capped at ~0 is a stall: the
     * martingale coin would pick the zero step with probability ~1, so
     * no unbiased progress is possible along this direction. (For
     * q == 0 the classic sum test is kept bit-identical.) */
    int stalled = (ws->q > 0)
        ? (lambda1 < tol || lambda2 < tol)
        : (lambda1 + lambda2 < tol);
    if (stalled) {
        /* If a binding row sits on its face, converting it to an
         * equality re-orients the next kernel direction — progress
         * without moving (and without a random draw, so E(s) = pik is
         * untouched). */
        int activated = 0;
        int binds[2] = {ws->bind1, ws->bind2};
        for (int b = 0; b < 2; b++) {
            int j = binds[b];
            if (j >= 0 && ws->row_state[j] == ROW_INACTIVE &&
                ws->slack[j] <= ROW_ACT_TOL(ws->rvec[j])) {
                cube_row_activate(ws, j);
                activated = 1;
            }
        }
        return activated;
    }

    double lambda = (unif_rand() * (lambda1 + lambda2) < lambda2)
                    ? lambda1 : -lambda2;

    for (int i = 0; i < n_cand; i++) {
        int id = ws->candidates[i];
        ws->prob[id] += lambda * ws->uvec[i];

        if (ws->prob[id] < eps) ws->prob[id] = 0.0;
        if (ws->prob[id] > 1.0 - eps) ws->prob[id] = 1.0;
    }

    if (ws->q > 0) {
        for (int j = 0; j < ws->q; j++) {
            double d = ws->dvec[j];
            if (d == 0.0) continue;  /* dropped, in kernel, or untouched */
            double g = ws->slack[j] - lambda * d;
            if (ws->row_state[j] == ROW_INACTIVE &&
                g <= ROW_ACT_TOL(ws->rvec[j])) {
                cube_row_activate(ws, j);
            } else {
                /* ROW_ACTIVE rows outside the kernel only occur during
                 * landing; their slack can grow back (the walk may leave
                 * the face inward — the bound itself still holds). */
                ws->slack[j] = g < 0.0 ? 0.0 : g;
            }
        }
    }

    return 1;
}

static void cube_flight(CubeWorkspace *ws) {
    double eps = ws->eps;
    int max_iter = (ws->N > 31622) ? 1000000000 : ws->N * ws->N;
    int iter = 0;

    if (ws->q == 0) {
        /* Legacy entry paths (stratified phases) manage prob/idx
         * directly; make sure the constraint count is in sync. */
        ws->n_eq = ws->p;
    }

    /* The window grows as inequality rows become equalities. */
    while (ws->idx->len >= ws->n_eq + 1) {
        R_CheckUserInterrupt();
        if (++iter > max_iter) {
            error("cube flight phase did not converge after %d iterations. "
                  "This usually indicates collinear balancing variables.",
                  max_iter);
        }

        int max_cand = ws->n_eq + 1;
        for (int i = 0; i < max_cand; i++) {
            ws->candidates[i] = ws->idx->list[i];
        }

        int progress = cube_update(ws, max_cand, ws->cons, ws->n_eq);

        for (int i = 0; i < max_cand; i++) {
            int id = ws->candidates[i];
            if (PROB_IS_INT(ws->prob[id], eps)) {
                indexlist_remove(ws->idx, id);
            }
        }

        /* With bounds, a no-progress step is deterministic (no RNG was
         * consumed), so retrying the same window cannot help. */
        if (ws->q > 0 && !progress) {
            error("cube flight phase stalled; this usually indicates "
                  "collinear balancing variables or badly scaled bounds");
        }
    }
}

static void cube_landing(CubeWorkspace *ws) {
    double eps = ws->eps;
    int max_iter = (ws->N > 31622) ? 1000000000 : ws->N * ws->N;
    int iter = 0;

    if (ws->q == 0) {
        ws->n_eq = ws->p;

        while (ws->idx->len > 1) {
            R_CheckUserInterrupt();
            if (++iter > max_iter) {
                error(
                    "cube landing phase did not converge after %d iterations.",
                    max_iter
                );
            }

            int n_cand = ws->idx->len;

            for (int i = 0; i < n_cand; i++) {
                ws->candidates[i] = ws->idx->list[i];
            }

            cube_update(ws, n_cand, ws->cons, n_cand - 1);

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
        return;
    }

    /* Landing with inequality bounds. The kernel keeps as many
     * constraints as the candidate count allows, in priority order:
     * the size constraint, then bounds that became equalities, then
     * auxiliary variables (relaxed from the last, as usual). All
     * remaining bounds keep capping the steps; a bound that provably
     * blocks completion is dropped (reported back for a warning). */
    while (ws->idx->len > 0) {
        R_CheckUserInterrupt();
        if (++iter > max_iter) {
            error("cube landing phase did not converge after %d iterations.",
                  max_iter);
        }

        int n_cand = ws->idx->len;
        for (int i = 0; i < n_cand; i++) {
            ws->candidates[i] = ws->idx->list[i];
        }

        /* Kernel priority: size constraint, then bound rows that are
         * tight NOW and act on the current candidates (rows whose walk
         * moved back inside, or with no support in the window, are
         * fully handled by the step caps), then auxiliary variables.
         * All of them are enforced at once — dependent rows (e.g.
         * redundant margin systems) cost nothing. Only when the system
         * has full column rank is the tail relaxed, one constraint at
         * a time: auxiliaries first (from the last, as in the classic
         * landing), then tight bounds, then the size constraint for
         * the final coin flip. */
        int m = 0;
        ws->landing_cons[m++] = 0;
        for (int i = ws->p; i < ws->n_eq; i++) {
            int c = ws->cons[i];
            int j = c - ws->p;
            if (ws->row_state[j] != ROW_ACTIVE) continue;
            if (ws->slack[j] > ROW_ACT_TOL(ws->rvec[j])) continue;
            int supported = 0;
            for (int k = 0; k < n_cand; k++) {
                if (ws->bcols[IDX_CM(ws->candidates[k], j, ws->N)] != 0.0) {
                    supported = 1;
                    break;
                }
            }
            if (supported) {
                ws->landing_cons[m++] = c;
            }
        }
        for (int c = 1; c < ws->p; c++) {
            ws->landing_cons[m++] = c;
        }

        int progress = -1;
        while (m >= 0) {
            progress = cube_update(ws, n_cand, ws->landing_cons, m);
            if (progress != -1) break;
            m--;
        }

        if (progress == 1) {
            for (int i = 0; i < n_cand; i++) {
                int id = ws->candidates[i];
                if (PROB_IS_INT(ws->prob[id], eps)) {
                    indexlist_remove(ws->idx, id);
                }
            }
            continue;
        }

        /* Stalled: relax the binding bound (there is no exact 0/1
         * completion satisfying it — e.g. an infeasible controlled
         * rounding structure). Unbiasedness is unaffected: no move was
         * made. If no bound binds, the stall is a degeneracy of the
         * balancing system itself; keep the classic behavior of
         * erroring out (via max_iter, message above). */
        int j = (ws->bind1 >= 0) ? ws->bind1 : ws->bind2;
        if (j >= 0) {
            ws->row_state[j] = ROW_DROPPED;
            ws->n_dropped++;
        }
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

/* Validate the inequality system (B, r) and return q = ncol(B).
 * The R layer always passes a double N x q matrix (q may be 0). */
static int cube_check_bounds(SEXP B_sexp, SEXP r_sexp, int N) {
    SEXP dim = getAttrib(B_sexp, R_DimSymbol);
    if (isNull(dim)) {
        error("'B' must be a numeric matrix");
    }
    int B_nrow = INTEGER(dim)[0];
    int q = INTEGER(dim)[1];
    if (q > 0 && B_nrow != N) {
        error("nrow(B) = %d does not match length(prob) = %d", B_nrow, N);
    }
    if (length(r_sexp) != q) {
        error("length(r) = %d does not match ncol(B) = %d",
              length(r_sexp), q);
    }
    return q;
}

/* --- Internal draw helpers (no RNG state management) -------------------- */

static int cube_draw_unstratified(const double *prob, const double *X,
                                  const double *B, const double *r,
                                  int N, int p, int q, double eps,
                                  int *out_idx, int *n_dropped) {
    CubeWorkspace *ws = cube_workspace_alloc(N, p, q, eps);
    if (q > 0) {
        ws->bcols = B;
        memcpy(ws->rvec, r, (size_t)q * sizeof(double));
    }
    cube_workspace_init(ws, prob, X);
    cube_ineq_reset(ws);
    cube_flight(ws);
    cube_landing(ws);
    *n_dropped = ws->n_dropped;
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

    CubeWorkspace *ws = cube_workspace_alloc(N, p_stratum, 0, eps);
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
                cube_workspace_alloc(N, p_global, 0, eps);
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

SEXP C_cube(SEXP prob_sexp, SEXP X_sexp, SEXP B_sexp, SEXP r_sexp,
            SEXP eps_sexp) {
    const int N = length(prob_sexp);
    int p = 0;
    cube_check_X_dim(X_sexp, N, &p);
    const int q = cube_check_bounds(B_sexp, r_sexp, N);

    const double eps = asReal(eps_sexp);
    const double *prob = REAL(prob_sexp);
    const double *X = REAL(X_sexp);
    const double *B = REAL(B_sexp);
    const double *r = REAL(r_sexp);
    int *selected = (int *)R_alloc(N, sizeof(int));
    int n_dropped = 0;

    GetRNGstate();
    int n_selected = cube_draw_unstratified(
        prob, X, B, r, N, p, q, eps, selected, &n_dropped
    );
    PutRNGstate();

    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
    if (n_dropped > 0) {
        setAttrib(result, install("relaxed"), ScalarInteger(n_dropped));
    }
    UNPROTECT(1);
    return result;
}

SEXP C_cube_batch(SEXP prob_sexp, SEXP X_sexp, SEXP B_sexp, SEXP r_sexp,
                  SEXP eps_sexp, SEXP nrep_sexp) {
    const int N = length(prob_sexp);
    int p = 0;
    cube_check_X_dim(X_sexp, N, &p);
    const int q = cube_check_bounds(B_sexp, r_sexp, N);

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
    CubeWorkspace *ws = cube_workspace_alloc(N, p, q, eps);
    if (q > 0) {
        ws->bcols = REAL(B_sexp);
        memcpy(ws->rvec, REAL(r_sexp), (size_t)q * sizeof(double));
    }
    cube_workspace_build_amat(ws, prob, X);

    int relaxed_reps = 0;

    GetRNGstate();
    for (int s = 0; s < nrep; s++) {
        if (s % 32 == 0) R_CheckUserInterrupt();

        cube_workspace_reset(ws, prob);
        cube_ineq_reset(ws);
        cube_flight(ws);
        cube_landing(ws);
        if (ws->n_dropped > 0) relaxed_reps++;

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

    if (relaxed_reps > 0) {
        setAttrib(result, install("relaxed_reps"), ScalarInteger(relaxed_reps));
    }
    UNPROTECT(1);
    return result;
}

/*
 * C_cube_stratified: Stratified cube method (Chauvet & Tillé, 2006/2009)
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
    CubeWorkspace *ws = cube_workspace_alloc(N, p_stratum, 0, eps);
    /* ws_global allocated with max p = n_strata + p_orig; actual p set per rep */
    CubeWorkspace *ws_global = cube_workspace_alloc(N, n_strata + p_orig, 0, eps);

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
