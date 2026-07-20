#ifndef SONDAGE_SPATIAL_CORE_H
#define SONDAGE_SPATIAL_CORE_H

#include <R.h>
#include <Rinternals.h>

/* Active unit ids with O(1) removal by id. */
typedef struct {
    int *list;
    int *reverse;
    int len;
    int capacity;
} SpatialPool;

static void spatial_pool_init(SpatialPool *pool, int capacity) {
    pool->list = (int *)R_alloc(capacity, sizeof(int));
    pool->reverse = (int *)R_alloc(capacity, sizeof(int));
    pool->len = 0;
    pool->capacity = capacity;
}

static SpatialPool *spatial_pool_alloc(int capacity) {
    SpatialPool *pool = (SpatialPool *)R_alloc(1, sizeof(SpatialPool));
    spatial_pool_init(pool, capacity);
    return pool;
}

static void spatial_pool_fill(SpatialPool *pool, const double *prob, int N) {
    pool->len = 0;
    for (int i = 0; i < N; i++) {
        if (prob[i] > 0.0 && prob[i] < 1.0) {
            pool->list[pool->len] = i;
            pool->reverse[i] = pool->len;
            pool->len++;
        } else {
            pool->reverse[i] = pool->capacity;
        }
    }
}

static void spatial_pool_remove(SpatialPool *pool, int id) {
    int position = pool->reverse[id];
    if (position >= pool->len) return;

    pool->len--;
    pool->reverse[id] = pool->capacity;
    if (position < pool->len) {
        int last = pool->list[pool->len];
        pool->list[position] = last;
        pool->reverse[last] = position;
    }
}

static int spatial_check_spread(SEXP spread, int N) {
    if (!isMatrix(spread) || !isReal(spread)) {
        error("'spread' must be a numeric matrix");
    }
    if (nrows(spread) != N) {
        error("nrow(spread) = %d does not match length(pik) = %d",
              nrows(spread), N);
    }
    int d = ncols(spread);
    if (d < 1) error("'spread' must have at least one column");
    return d;
}

#endif
