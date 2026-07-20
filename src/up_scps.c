/*
 * up_scps.c: Spatially correlated Poisson sampling (Grafstrom, 2012)
 *
 * This is an independent implementation of the maximal-weight strategy
 * described in:
 *   Grafstrom, A. (2012). Spatially correlated Poisson sampling.
 *   Journal of Statistical Planning and Inference, 142(1), 139-147.
 *
 * At each step an undecided unit i is chosen uniformly at random and its
 * outcome is fixed by a Bernoulli draw with its current conditional
 * probability p_i.  The resulting probability displacement is distributed
 * to nearby undecided units with non-negative weights that sum to one.
 * Each neighbour receives as much weight as its feasibility bound allows,
 * beginning with the closest distance group.  Equal-distance units share
 * the remaining weight as evenly as their bounds permit.
 *
 * Finding the distance at which the cumulative feasible weight reaches one
 * is a weighted-selection problem.  We solve it with randomized three-way
 * quickselect, avoiding a full sort at every step.  A draw therefore costs
 * expected O(N^2 * d) time and O(N) workspace, and stores no distance matrix.
 */

#include "sampling_core.h"
#include "spatial_core.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int id;
    double distance;
    double capacity;
    double weight;
} ScpsCandidate;

typedef struct {
    SpatialPool pool;
    ScpsCandidate *candidate;
} ScpsWorkspace;

static void scps_workspace_init(ScpsWorkspace *work, int N) {
    spatial_pool_init(&work->pool, N);
    work->candidate = (ScpsCandidate *)R_alloc(N, sizeof(ScpsCandidate));
}

static void scps_swap(ScpsCandidate *left, ScpsCandidate *right) {
    ScpsCandidate tmp = *left;
    *left = *right;
    *right = tmp;
}

/* Return the smallest distance whose cumulative candidate capacity reaches
 * target.  Equal distances are kept together, which is important on grids.
 * The randomized pivot gives expected linear work per selection. */
static double scps_weighted_cutoff(ScpsCandidate *candidate, int n,
                                   double target) {
    int lo = 0;
    int hi = n;

    while (hi - lo > 1) {
        int pivot_position = lo + (int)(unif_rand() * (hi - lo));
        if (pivot_position >= hi) pivot_position = hi - 1;
        double pivot = candidate[pivot_position].distance;

        int less = lo;
        int scan = lo;
        int greater = hi;
        double less_capacity = 0.0;
        double equal_capacity = 0.0;

        while (scan < greater) {
            double value = candidate[scan].distance;
            if (value < pivot) {
                less_capacity += candidate[scan].capacity;
                scps_swap(candidate + less, candidate + scan);
                less++;
                scan++;
            } else if (value > pivot) {
                greater--;
                scps_swap(candidate + scan, candidate + greater);
            } else {
                equal_capacity += candidate[scan].capacity;
                scan++;
            }
        }

        /* If the lower group alone can carry target, continue there.  The
         * tolerance only handles a target equal to the lower capacity up to
         * roundoff; choosing pivot then simply assigns zero to its tie group. */
        double tol = 64.0 * DBL_EPSILON * (1.0 + target);
        if (target < less_capacity - tol) {
            hi = less;
        } else if (target <= less_capacity + equal_capacity + tol) {
            return pivot;
        } else {
            target -= less_capacity + equal_capacity;
            lo = greater;
        }
    }

    return candidate[lo].distance;
}

static int scps_capacity_compare(const void *left, const void *right) {
    double a = ((const ScpsCandidate *)left)->capacity;
    double b = ((const ScpsCandidate *)right)->capacity;
    return (a > b) - (a < b);
}

/* Assign maximal weights once the cutoff distance is known. */
static void scps_assign_weights(ScpsCandidate *candidate, int n,
                                double cutoff) {
    int closer = 0;
    int scan = 0;
    int farther = n;

    /* Regroup so the cutoff-distance candidates are contiguous. */
    while (scan < farther) {
        double value = candidate[scan].distance;
        if (value < cutoff) {
            scps_swap(candidate + closer, candidate + scan);
            closer++;
            scan++;
        } else if (value > cutoff) {
            farther--;
            scps_swap(candidate + scan, candidate + farther);
        } else {
            scan++;
        }
    }

    double remaining = 1.0;
    for (int k = 0; k < closer; k++) {
        candidate[k].weight = candidate[k].capacity;
        remaining -= candidate[k].weight;
    }
    for (int k = farther; k < n; k++) candidate[k].weight = 0.0;

    if (remaining < 0.0 && remaining > -128.0 * DBL_EPSILON) {
        remaining = 0.0;
    }

    int tied = farther - closer;
    qsort(candidate + closer, (size_t)tied, sizeof(ScpsCandidate),
          scps_capacity_compare);

    /* Water filling: equal weights within the cutoff group, except where a
     * unit's feasibility capacity is too small to receive the equal share. */
    int first_equal = closer;
    while (first_equal < farther) {
        int left = farther - first_equal;
        double share = remaining / left;
        if (candidate[first_equal].capacity <= share) {
            candidate[first_equal].weight = candidate[first_equal].capacity;
            remaining -= candidate[first_equal].weight;
            first_equal++;
        } else {
            for (int k = first_equal; k < farther; k++) {
                candidate[k].weight = share;
            }
            remaining = 0.0;
            break;
        }
    }

    /* Absorb the final few ulps so the weights sum exactly to one in the
     * arithmetic subsequently used by the probability update. */
    if (remaining != 0.0 && farther > closer) {
        candidate[farther - 1].weight += remaining;
    }
}

static void scps_settle(SpatialPool *pool, double *prob, int id, double eps) {
    if (prob[id] <= eps) {
        prob[id] = 0.0;
        spatial_pool_remove(pool, id);
    } else if (prob[id] >= 1.0 - eps) {
        prob[id] = 1.0;
        spatial_pool_remove(pool, id);
    }
}

/* Restore the probability mass lost or gained through floating-point
 * arithmetic and boundary snapping in one update.  In exact arithmetic the
 * step-unit displacement and the weighted neighbour displacements cancel.
 * Keeping that invariant numerically prevents tiny errors from accumulating
 * until a later maximal-weight system appears infeasible. */
static int scps_correct_mass(SpatialPool *pool, double *prob, double correction,
                             double eps) {
    double tolerance = 4.0 * DBL_EPSILON * (1.0 + fabs(correction));

    while (fabs(correction) > tolerance && pool->len > 0) {
        int best_id = -1;
        double best_room = 0.0;
        for (int k = 0; k < pool->len; k++) {
            int id = pool->list[k];
            double room = (correction > 0.0) ? 1.0 - prob[id] : prob[id];
            if (room > best_room) {
                best_room = room;
                best_id = id;
            }
        }

        if (best_id < 0 || best_room <= 0.0) return 0;
        double amount = fmin(fabs(correction), best_room);
        double change = (correction > 0.0) ? amount : -amount;
        double before = prob[best_id];
        prob[best_id] += change;
        double actual_change = prob[best_id] - before;
        if (actual_change == 0.0) return fabs(correction) <= 1e-12;
        correction -= actual_change;
        scps_settle(pool, prob, best_id, eps);
    }

    return fabs(correction) <= 1e-12;
}

static int scps_build_candidates(ScpsWorkspace *work, const double *prob,
                                 const double *spread, int N, int d,
                                 int step_id, double step_prob) {
    int n = 0;
    for (int k = 0; k < work->pool.len; k++) {
        int id = work->pool.list[k];
        if (id == step_id) continue;

        double distance = 0.0;
        for (int column = 0; column < d; column++) {
            double delta = spread[id + column * N] -
                           spread[step_id + column * N];
            distance += delta * delta;
        }

        double p = prob[id];
        double capacity_selected = (1.0 - p) / step_prob;
        double capacity_rejected = p / (1.0 - step_prob);

        work->candidate[n].id = id;
        work->candidate[n].distance = distance;
        work->candidate[n].capacity = fmin(capacity_selected,
                                            capacity_rejected);
        work->candidate[n].weight = 0.0;
        n++;
    }
    return n;
}

static int scps_run(double *prob, const double *spread, int N, int d,
                    double eps, ScpsWorkspace *work) {
    spatial_pool_fill(&work->pool, prob, N);

    int steps = 0;
    while (work->pool.len > 1) {
        if ((++steps & 255) == 0) R_CheckUserInterrupt();

        int position = (int)(unif_rand() * work->pool.len);
        if (position >= work->pool.len) position = work->pool.len - 1;
        int step_id = work->pool.list[position];
        double step_prob = prob[step_id];

        int n_candidates = scps_build_candidates(
            work, prob, spread, N, d, step_id, step_prob
        );

        double total_capacity = 0.0;
        for (int k = 0; k < n_candidates; k++) {
            total_capacity += work->candidate[k].capacity;
        }

        /* With an integer probability sum this condition follows from the
         * CPS feasibility bounds.  Allow a small floating-point deficit. */
        double feasibility_tolerance =
            128.0 * work->pool.capacity * DBL_EPSILON;
        if (total_capacity < 1.0 - feasibility_tolerance) {
            return 0;
        }

        double cutoff_target = fmin(1.0, total_capacity);
        double cutoff = scps_weighted_cutoff(
            work->candidate, n_candidates, cutoff_target
        );
        scps_assign_weights(work->candidate, n_candidates, cutoff);

        double outcome = (unif_rand() < step_prob) ? 1.0 : 0.0;
        double displacement = outcome - step_prob;
        double mass_change = displacement;
        prob[step_id] = outcome;
        spatial_pool_remove(&work->pool, step_id);

        for (int k = 0; k < n_candidates; k++) {
            double weight = work->candidate[k].weight;
            if (weight == 0.0) continue;
            int id = work->candidate[k].id;
            double old_prob = prob[id];
            prob[id] -= displacement * weight;
            scps_settle(&work->pool, prob, id, eps);
            mass_change += prob[id] - old_prob;
        }

        if (!scps_correct_mass(
                &work->pool, prob, -mass_change, eps
            )) {
            return 0;
        }
    }

    /* A lone unresolved probability can only be a rounding residue because
     * the initial total is an integer and every maximal-weight update keeps
     * the total fixed. */
    if (work->pool.len == 1) {
        int id = work->pool.list[0];
        prob[id] = (prob[id] >= 0.5) ? 1.0 : 0.0;
        spatial_pool_remove(&work->pool, id);
    }
    return 1;
}

SEXP C_scps(SEXP pik_sexp, SEXP spread_sexp, SEXP eps_sexp) {
    const int N = length(pik_sexp);
    const int d = spatial_check_spread(spread_sexp, N);

    const double eps = asReal(eps_sexp);
    double *prob = (double *)R_alloc(N, sizeof(double));
    int *selected = (int *)R_alloc(N, sizeof(int));
    memcpy(prob, REAL(pik_sexp), (size_t)N * sizeof(double));

    ScpsWorkspace work;
    scps_workspace_init(&work, N);

    GetRNGstate();
    int feasible = scps_run(prob, REAL(spread_sexp), N, d, eps, &work);
    PutRNGstate();
    if (!feasible) error("SCPS maximal weights are numerically infeasible");

    int n_selected = sampling_extract_selected(prob, N, 0.5, selected, N);
    SEXP result = PROTECT(allocVector(INTSXP, n_selected));
    memcpy(INTEGER(result), selected, (size_t)n_selected * sizeof(int));
    UNPROTECT(1);
    return result;
}

SEXP C_scps_batch(SEXP pik_sexp, SEXP spread_sexp, SEXP eps_sexp,
                  SEXP nrep_sexp) {
    const int N = length(pik_sexp);
    const int d = spatial_check_spread(spread_sexp, N);

    const int nrep = INTEGER(nrep_sexp)[0];
    const double eps = asReal(eps_sexp);
    const double *pik = REAL(pik_sexp);
    const double *spread = REAL(spread_sexp);

    double total = 0.0;
    for (int i = 0; i < N; i++) total += pik[i];
    const int n_target = (int)(total + 0.5);

    SEXP result = PROTECT(allocMatrix(INTSXP, n_target, nrep));
    int *out = INTEGER(result);
    double *prob = (double *)R_alloc(N, sizeof(double));
    ScpsWorkspace work;
    scps_workspace_init(&work, N);

    GetRNGstate();
    for (int draw = 0; draw < nrep; draw++) {
        if ((draw & 31) == 0) R_CheckUserInterrupt();
        memcpy(prob, pik, (size_t)N * sizeof(double));
        int feasible = scps_run(prob, spread, N, d, eps, &work);
        if (!feasible) {
            PutRNGstate();
            UNPROTECT(1);
            error("SCPS maximal weights are numerically infeasible in draw %d",
                  draw + 1);
        }

        int *column = out + (R_xlen_t)draw * n_target;
        int n_selected = sampling_extract_selected(
            prob, N, 0.5, column, n_target
        );
        if (n_selected != n_target) {
            PutRNGstate();
            UNPROTECT(1);
            error("SCPS batch draw %d produced size %d, expected %d",
                  draw + 1, n_selected, n_target);
        }
    }
    PutRNGstate();

    UNPROTECT(1);
    return result;
}
