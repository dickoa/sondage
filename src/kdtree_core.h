/*
 * kdtree_core.h: exact nearest-neighbor search over a shrinking set of
 * units, for the local pivotal method.
 *
 * A static k-d tree (Friedman, Bentley & Finkel, 1977) is built once over
 * the units that start undecided. Each node splits its units at the
 * median of the dimension with the largest extent and keeps the bounding
 * box of its units and a count of live units below it. Removing a unit
 * marks it dead and decrements the counts on its path to the root, so
 * empty subtrees are skipped without rebuilding.
 *
 * A query is an exact branch-and-bound search that returns every live
 * unit at the minimum distance, leaving the tie break to the caller.
 * Distances come from kd_dist(), which the linear scan in lpm.c uses as
 * well, so both searches see the same values. A subtree is pruned only
 * when the lower bound from its box exceeds the best distance by a
 * relative margin of a few ulps. Each squared gap in the bound is no
 * larger than the matching squared coordinate difference of any unit in
 * the box, and rounding is monotone, so the bound never exceeds a
 * computed distance. The margin keeps that true if a compiler contracts
 * the two sums differently into fused multiply-adds.
 *
 * The tree is clean-room code written from the paper:
 *   Friedman, J.H., Bentley, J.L. and Finkel, R.A. (1977). An algorithm
 *     for finding best matches in logarithmic expected time. ACM
 *     Transactions on Mathematical Software, 3(3), 209-226.
 *
 * Building costs O(N log N) expected time and O(N) memory. No random
 * numbers are used.
 */

#ifndef SONDAGE_KDTREE_CORE_H
#define SONDAGE_KDTREE_CORE_H

#include <R.h>
#include <Rinternals.h>
#include <float.h>
#include <string.h>

#define KD_LEAF_SIZE 8
#define KD_STACK_SIZE 256

typedef struct {
    int lo, hi;      /* the node's units are perm[lo .. hi - 1] */
    int left, right; /* children, or -1 at a leaf */
    int parent;      /* -1 at the root */
    int live;        /* live units in the subtree */
} KdNode;

typedef struct {
    const double *z; /* column-major N x d coordinates (borrowed) */
    int N;
    int d;
    int n_units;     /* units in the tree */
    int *perm;       /* unit ids grouped by leaf */
    int *leaf_of;    /* unit id -> leaf node, or -1 if not in the tree */
    unsigned char *alive;  /* unit id -> 1 while live */
    unsigned char *alive0; /* alive at build time, for kd_reset() */
    KdNode *node;
    int n_node;
    int cap_node;
    int *live0;      /* node live counts at build time */
    double *bmin;    /* node bounding boxes, n_node x d, row-major */
    double *bmax;
    double margin;   /* relative pruning margin */
    int *ties;       /* units at the minimum distance after kd_nearest() */
    int n_ties;
    int *stack;
    double *stack_lb;
} KdTree;

/* Squared distance from unit b to unit a, summed over coordinates in
 * order. Stops early, returning the partial sum, once it exceeds
 * `bound`. */
static inline double kd_dist(const double *z, int N, int d, int a, int b,
                             double bound) {
    double dist = 0.0;
    for (int c = 0; c < d; c++) {
        double diff = z[a + (R_xlen_t)c * N] - z[b + (R_xlen_t)c * N];
        dist += diff * diff;
        if (dist > bound) break;
    }
    return dist;
}

/* Rearrange perm[lo .. hi - 1] so that perm[kth] holds the unit with the
 * kth smallest key, smaller-or-equal keys before it and greater-or-equal
 * after. Hoare partition with a median-of-three pivot, no RNG. */
static void kd_select(int *perm, const double *key, int lo, int hi,
                      int kth) {
    hi--;
    while (hi > lo) {
        int mid = lo + (hi - lo) / 2;
        double a = key[perm[lo]], b = key[perm[mid]], c = key[perm[hi]];
        double pivot;
        if (a < b) {
            pivot = (b < c) ? b : ((a < c) ? c : a);
        } else {
            pivot = (a < c) ? a : ((b < c) ? c : b);
        }
        int i = lo, j = hi;
        while (i <= j) {
            while (key[perm[i]] < pivot) i++;
            while (key[perm[j]] > pivot) j--;
            if (i <= j) {
                int tmp = perm[i];
                perm[i] = perm[j];
                perm[j] = tmp;
                i++;
                j--;
            }
        }
        if (kth <= j) {
            hi = j;
        } else if (kth >= i) {
            lo = i;
        } else {
            return;
        }
    }
}

static int kd_build_node(KdTree *t, int lo, int hi, int parent) {
    if (t->n_node >= t->cap_node) {
        error("internal error: k-d tree node capacity exceeded");
    }
    const int k = t->n_node++;
    const int d = t->d;
    const int N = t->N;
    double *bmin = t->bmin + (size_t)k * d;
    double *bmax = t->bmax + (size_t)k * d;

    for (int c = 0; c < d; c++) {
        const double *zc = t->z + (R_xlen_t)c * N;
        double lo_c = zc[t->perm[lo]], hi_c = lo_c;
        for (int i = lo + 1; i < hi; i++) {
            double v = zc[t->perm[i]];
            if (v < lo_c) lo_c = v;
            if (v > hi_c) hi_c = v;
        }
        bmin[c] = lo_c;
        bmax[c] = hi_c;
    }

    t->node[k].lo = lo;
    t->node[k].hi = hi;
    t->node[k].parent = parent;
    t->node[k].live = hi - lo;
    t->node[k].left = -1;
    t->node[k].right = -1;

    if (hi - lo <= KD_LEAF_SIZE) {
        for (int i = lo; i < hi; i++) t->leaf_of[t->perm[i]] = k;
        return k;
    }

    int split = 0;
    double widest = bmax[0] - bmin[0];
    for (int c = 1; c < d; c++) {
        if (bmax[c] - bmin[c] > widest) {
            widest = bmax[c] - bmin[c];
            split = c;
        }
    }
    const int mid = lo + (hi - lo) / 2;
    kd_select(t->perm, t->z + (R_xlen_t)split * N, lo, hi, mid);

    int left = kd_build_node(t, lo, mid, k);
    int right = kd_build_node(t, mid, hi, k);
    t->node[k].left = left;
    t->node[k].right = right;
    return k;
}

/* Build the tree over the units with 0 < prob < 1 (the units a pivotal
 * draw starts with undecided). All memory comes from R_alloc. */
static KdTree *kd_build(const double *z, int N, int d, const double *prob) {
    KdTree *t = (KdTree *)R_alloc(1, sizeof(KdTree));
    t->z = z;
    t->N = N;
    t->d = d;
    t->perm = (int *)R_alloc(N > 0 ? N : 1, sizeof(int));
    t->leaf_of = (int *)R_alloc(N > 0 ? N : 1, sizeof(int));
    t->alive = (unsigned char *)R_alloc(N > 0 ? N : 1, 1);
    t->alive0 = (unsigned char *)R_alloc(N > 0 ? N : 1, 1);
    t->ties = (int *)R_alloc(N > 0 ? N : 1, sizeof(int));
    t->stack = (int *)R_alloc(KD_STACK_SIZE, sizeof(int));
    t->stack_lb = (double *)R_alloc(KD_STACK_SIZE, sizeof(double));
    t->n_ties = 0;
    t->margin = 1.0 + 8.0 * (double)(d + 1) * DBL_EPSILON;

    int n = 0;
    for (int i = 0; i < N; i++) {
        t->leaf_of[i] = -1;
        if (prob[i] > 0.0 && prob[i] < 1.0) {
            t->perm[n++] = i;
            t->alive[i] = 1;
        } else {
            t->alive[i] = 0;
        }
    }
    t->n_units = n;

    /* Splits at the midpoint keep every leaf above KD_LEAF_SIZE / 2
     * units, so a tree has fewer than 4 * n / KD_LEAF_SIZE + 2 nodes. */
    t->cap_node = 4 * (n / KD_LEAF_SIZE) + 4;
    t->node = (KdNode *)R_alloc(t->cap_node, sizeof(KdNode));
    t->bmin = (double *)R_alloc((size_t)t->cap_node * d, sizeof(double));
    t->bmax = (double *)R_alloc((size_t)t->cap_node * d, sizeof(double));
    t->live0 = (int *)R_alloc(t->cap_node, sizeof(int));
    t->n_node = 0;
    if (n > 0) kd_build_node(t, 0, n, -1);

    for (int k = 0; k < t->n_node; k++) t->live0[k] = t->node[k].live;
    memcpy(t->alive0, t->alive, (size_t)N);
    return t;
}

/* Restore the build-time live set, for the next replicate. */
static void kd_reset(KdTree *t) {
    memcpy(t->alive, t->alive0, (size_t)t->N);
    for (int k = 0; k < t->n_node; k++) t->node[k].live = t->live0[k];
}

static void kd_remove(KdTree *t, int id) {
    if (!t->alive[id]) return;
    t->alive[id] = 0;
    for (int k = t->leaf_of[id]; k >= 0; k = t->node[k].parent) {
        t->node[k].live--;
    }
}

/* Squared distance from unit q to the box of node k, a lower bound on
 * kd_dist() from q to any unit in the box. */
static inline double kd_box_bound(const KdTree *t, int k, int q) {
    const int d = t->d;
    const double *bmin = t->bmin + (size_t)k * d;
    const double *bmax = t->bmax + (size_t)k * d;
    double bound = 0.0;
    for (int c = 0; c < d; c++) {
        double x = t->z[q + (R_xlen_t)c * t->N];
        double gap = 0.0;
        if (x < bmin[c]) {
            gap = bmin[c] - x;
        } else if (x > bmax[c]) {
            gap = x - bmax[c];
        }
        bound += gap * gap;
    }
    return bound;
}

/* Collect in t->ties every live unit other than q at the minimum
 * distance from q, in tree order. Returns the number of ties (0 when q
 * is the only live unit). */
static int kd_nearest(KdTree *t, int q) {
    const double *z = t->z;
    const int N = t->N;
    const int d = t->d;
    double best = R_PosInf;
    int top = 0;

    t->n_ties = 0;
    if (t->n_node == 0) return 0;
    t->stack[top] = 0;
    t->stack_lb[top] = 0.0;
    top++;

    while (top > 0) {
        top--;
        const int k = t->stack[top];
        const KdNode *nd = &t->node[k];
        if (nd->live == 0 || t->stack_lb[top] > best * t->margin) continue;

        if (nd->left < 0) {
            for (int i = nd->lo; i < nd->hi; i++) {
                const int id = t->perm[i];
                if (!t->alive[id] || id == q) continue;
                double dist = kd_dist(z, N, d, id, q, best);
                if (dist < best) {
                    best = dist;
                    t->ties[0] = id;
                    t->n_ties = 1;
                } else if (dist == best) {
                    t->ties[t->n_ties++] = id;
                }
            }
            continue;
        }

        const int l = nd->left, r = nd->right;
        double lb_l = kd_box_bound(t, l, q);
        double lb_r = kd_box_bound(t, r, q);
        if (top + 2 > KD_STACK_SIZE) {
            error("internal error: k-d tree search stack overflow");
        }
        /* Push the farther child first so the nearer one is searched
         * first and tightens `best` sooner. */
        if (lb_l <= lb_r) {
            t->stack[top] = r;
            t->stack_lb[top++] = lb_r;
            t->stack[top] = l;
            t->stack_lb[top++] = lb_l;
        } else {
            t->stack[top] = l;
            t->stack_lb[top++] = lb_l;
            t->stack[top] = r;
            t->stack_lb[top++] = lb_r;
        }
    }
    return t->n_ties;
}

#endif
