#' Balanced Sampling Without Replacement
#'
#' Draws a balanced sample using the cube method (Deville & \enc{Tillé}{Tille}, 2004),
#' or a spatially balanced (well-spread) sample using the local pivotal
#' method (\enc{Grafström}{Grafstrom}, \enc{Lundström}{Lundstrom} & Schelin,
#' 2012) or spatially correlated Poisson sampling
#' (\enc{Grafström}{Grafstrom}, 2012).
#' A balanced sample satisfies (approximately) the balancing equations
#' \eqn{\sum_{k \in S} x_k / \pi_k \approx \sum_{k \in U} x_k}{sum_S x_k / pi_k ~ sum_U x_k} for each
#' auxiliary variable \eqn{x}. A well-spread sample selects units that
#' are far apart in the space spanned by the spreading variables.
#'
#' @param pik A numeric vector of inclusion probabilities (length N).
#'   `sum(pik)` must be an integer to floating-point accuracy. See
#'   [unequal_prob_wor()] for the exact-0/1 handling of boundary
#'   values. The target sample size, `sum(pik)`, must be at least 1.
#' @param aux An optional numeric matrix (N x p) of auxiliary balancing
#'   variables. Each column defines a balancing constraint. The sample
#'   size constraint is always included automatically and `aux` specifies
#'   additional variables to balance on.
#'   When `NULL`, only the sample size is balanced (equivalent to an
#'   unbalanced fixed-size design).
#' @param strata An optional integer vector (length N) of stratum
#'   indicators (positive integers). Uses the stratified cube method
#'   (Chauvet & \enc{Tillé}{Tille}, 2006; Chauvet, 2009) to preserve within-stratum
#'   sample sizes while balancing on `aux`. Requires `sum(pik)` within
#'   each stratum to be close to an integer for exact sizes. If not, a
#'   warning is issued and `fixed_size` is set to `FALSE`.
#' @param spread An optional numeric matrix (N x d) of spatial
#'   coordinates (or other spreading variables) for well-spread,
#'   spatially balanced sampling. Required by the built-in `"lpm2"`
#'   and `"scps"` methods, and supported by methods registered via
#'   [register_method()] with `supports_spread = TRUE`. The built-in
#'   `"cube"` method does not use it and will error. A non-matrix
#'   vector is treated as a single spreading variable.
#' @param bounds An optional list with elements `B`, `lower`, and
#'   `upper` describing linear inequality constraints on the realized
#'   sample (Tripet & \enc{Tillé}{Tille}, 2026): `B` is a numeric
#'   matrix (N x q) whose columns are constraint variables, and the
#'   sample `s` is drawn so that
#'   \eqn{lower_j \le \sum_{k \in s} B_{kj} \le upper_j}{lower_j <= sum_s B_kj <= upper_j} for every
#'   constraint `j`. `-Inf` / `Inf` entries make a constraint
#'   one-sided, and `lower[j] == upper[j]` enforces an exact equality.
#'   The starting probabilities must be feasible:
#'   `lower <= colSums(B * pik) <= upper`. Only supported by the
#'   built-in `"cube"` method. See **Inequality constraints** below.
#' @param method The sampling method. `"cube"` (the default) balances
#'   on `aux`. `"lpm2"` (local pivotal method 2) and `"scps"`
#'   (spatially correlated Poisson sampling) spread on `spread`. The
#'   name of a balanced method added via [register_method()] is also
#'   accepted.
#' @param nrep Number of replicate samples (default 1). When `nrep > 1`,
#'   `$sample` holds a matrix (n x nrep) for fixed-size designs, or a
#'   list of integer vectors when within-stratum sizes are not exact.
#' @param prn Optional vector of permanent random numbers (length N,
#'   values in the open interval (0, 1)) for sample coordination.
#'   Supported by `"scps"` and by methods registered with
#'   `supports_prn = TRUE`. With `"scps"`, the sample is then a function
#'   of `pik`, `spread`, `prn` and the row order alone. See **Sample
#'   coordination** below. Cannot be combined with `nrep > 1`.
#' @param ... Additional arguments passed to methods:
#'   \describe{
#'     \item{`eps`}{Boundary tolerance (default `1e-10`): decides when
#'       an updated \emph{working} probability has numerically reached
#'       0 or 1 (during the cube flight phase, or during the pivotal
#'       steps of `"lpm2"` or `"scps"`). It never reclassifies the input.
#'       Supplied `pik` inside `(0, eps]` or `[1 - eps, 1)` are rejected. Only
#'       `pik` of exactly 0 or exactly 1 enter the algorithm as
#'       already-resolved units.}
#'     \item{`condition_aux`}{Logical. If `TRUE`, pre-conditions `aux` to
#'       improve numerical stability with ill-conditioned or collinear
#'       auxiliary variables (default `FALSE`). Each column is centered and
#'       scaled in the space the cube balances, `aux / pik`, so that the
#'       balancing equations are unchanged, and columns that depend on
#'       earlier ones or on the sample size constraint are dropped. The
#'       other columns keep their order.}
#'     \item{`qr_tol`}{Tolerance for QR rank detection when
#'       `condition_aux = TRUE`. It must be one finite, non-negative
#'       number (default `sqrt(.Machine$double.eps)`).}
#'   }
#'
#' @details
#' The cube method proceeds in two phases:
#' \describe{
#'   \item{Flight phase}{Probabilities are moved toward 0 or 1 while
#'     maintaining all balancing constraints. Each step resolves at
#'     least one unit. Terminates when fewer than p+1 undecided units
#'     remain.}
#'   \item{Landing phase}{Remaining undecided units are resolved by
#'     progressively relaxing balancing constraints, starting from the
#'     \strong{last} column of `aux`. Users should order auxiliary
#'     variables by importance (most important first).}
#' }
#'
#' The sample size constraint is always placed first (never relaxed
#' during landing). For stratified designs, within-stratum size
#' constraints are also placed first.
#'
#' Joint inclusion probabilities are approximated via the high-entropy
#' approximation (Brewer & Donadio, 2003). It uses only `pik`, so it ignores
#' `aux`, `strata`, and `bounds`, and it returns the same matrix as for an
#' unconstrained design. It can give a positive probability to a pair the
#' design excludes. For example, with strata of two units and one selection
#' in each, two units of the same stratum never occur together. The
#' marginal-defect warning of [joint_inclusion_prob()] does not detect this.
#' For a stratified or strongly constrained design, prefer Monte Carlo
#' estimates of the joint probabilities or a variance method suited to
#' balanced sampling.
#'
#' @section Inequality constraints:
#'
#' `bounds` implements the cube method with inequality constraints of
#' Tripet & \enc{Tillé}{Tille} (2026). During the flight phase, steps
#' are capped so every constraint stays feasible, and a constraint
#' whose slack reaches zero becomes an equality from then on. The
#' inclusion probabilities are respected exactly (`E(s) = pik`),
#' unlike rejective procedures.
#'
#' `B` applies to the realized sample directly (counts or raw sums),
#' so it is \strong{not} divided by `pik`. To bound a Horvitz-Thompson
#' estimator, pass `x / pik` as the constraint column.
#'
#' One application is controlled selection with bounding
#' category counts to the integers adjacent to their expectation.
#' With indicator columns `B` and `S <- colSums(B * pik)`, use `lower = floor(S)`,
#' `upper = ceiling(S)`. Categories may overlap (e.g. row and column
#' margins of a two-way control table, as in NAEP-style designs).
#'
#' Integer-valued bound systems on partitions or two-way margins are
#' satisfied exactly. For structures with no exact integer solution,
#' some three-way controlled rounding problems, or continuous-valued
#' constraints that end the flight phase tight against a boundary,
#' bounds that provably block the landing phase are relaxed one at a
#' time, with a warning, and `E(s) = pik` still holds. In other words,
#' the bounds are guaranteed whenever no relaxation warning is
#' raised.
#'
#' The high-entropy approximation used by [joint_inclusion_prob()] is
#' less accurate under tight bounds, which distort the design away
#' from maximum entropy. Tripet & \enc{Tillé}{Tille} (2026) recommend
#' Monte Carlo estimation of joint inclusion probabilities in that
#' case.
#'
#' @section Spatially balanced sampling (`lpm2` and `scps`):
#'
#' `method = "lpm2"` implements the local pivotal method 2 of
#' \enc{Grafström}{Grafstrom}, \enc{Lundström}{Lundstrom} & Schelin
#' (2012). Repeatedly, a randomly chosen undecided unit and its
#' nearest undecided neighbor in the `spread` space compete in a
#' pivotal step (Deville & \enc{Tillé}{Tille}, 1998) that resolves at
#' least one of them to 0 or 1 while preserving the inclusion
#' probabilities exactly (`E(s) = pik`). Nearby units thereby tend to
#' exclude each other, spreading the sample over the population.
#' Nearest-neighbor distance ties (common on gridded coordinates)
#' are broken uniformly at random.
#'
#' Spreading variables should be on comparable scales, since
#' nearness is plain Euclidean distance in the `spread` columns.
#' Rescale them (e.g. with [scale()]) when they are not. `"lpm2"` and
#' `"scps"` use `spread` only: they do not accept `aux`, `strata`, or
#' `bounds`. To exactly balance covariate totals \emph{and} spread
#' spatially, register a method that supports both `aux` and `spread`
#' (for example, a local-cube implementation). Running `"cube"` with
#' coordinates in `aux` balances their totals but does not itself enforce
#' spatial spread.
#'
#' `method = "scps"` implements spatially correlated Poisson sampling
#' with the maximal-weight strategy of \enc{Grafström}{Grafstrom} (2012).
#' At each step, a randomly chosen undecided unit is accepted or rejected
#' using its current conditional probability. Its probability displacement
#' is then distributed to the nearest undecided units, subject to feasibility
#' bounds that keep every working probability in \eqn{[0, 1]}.
#' Equal-distance
#' units share weight as evenly as their bounds allow. Random selection of
#' the step unit avoids dependence on input row order. With `prn`, units
#' are visited in row order instead (see **Sample coordination** below).
#'
#' Both spatial methods deliberately drive joint inclusion probabilities of
#' nearby units toward zero, so the design is \emph{not} high
#' entropy and no joint-probability approximation is provided:
#' [joint_inclusion_prob()] errors for these designs. Variance for
#' well-spread samples is usually estimated with local-neighborhood
#' estimators (e.g. \enc{Grafström}{Grafstrom} & Schelin, 2014).
#'
#' Both methods are \emph{spread-only}: neither exactly balances the totals
#' of `aux`. They are dispatched by `balanced_wor()` because "spatially
#' balanced sampling" is the standard name for well-spread fixed-size
#' designs, and because the same interface accommodates local-cube methods
#' that combine exact balancing with spread. Capability metadata keeps the
#' distinction explicit: both report `supports_aux = FALSE` and
#' `supports_spread = TRUE` through [method_spec()].
#'
#' LPM2 finds nearest neighbors with a k-d tree (Friedman et al., 1977) when `spread` has at most
#' eight columns, so a draw costs about O(N log N) time. Many units at
#' the same distance, as with duplicated points, make the tree visit all of
#' them, and a draw then approaches O(N^2). With more columns
#' a tree prunes too little and LPM2 scans the undecided units instead, at
#' O(N^2 * d) per draw. SCPS uses weighted quickselect to
#' find the distance at which its maximal weights sum to one, avoiding a
#' full sort of the remaining units at each step. Its expected cost is
#' O(N^2 * d). Both implementations use O(N) workspace and store no
#' distance matrix. Units at the cutoff distance share their weight
#' equally, up to their feasibility bounds, and SCPS sorts that group only
#' when one of its units cannot take an equal share.
#'
#' @section Sample coordination:
#'
#' With `prn`, `method = "scps"` is the list-sequential correlated Poisson
#' scheme of Bondesson & Thorburn (2008) with
#' \enc{Grafström}{Grafstrom}'s maximal weights, coordinated through
#' permanent random numbers as in \enc{Grafström}{Grafstrom} & Matei
#' (2018). Units are visited in row order. Unit \eqn{k} is selected when
#' \eqn{u_k < \pi_k^{(k-1)}}{u_k < pi_k^(k-1)}, its current conditional
#' probability, and its
#' probability displacement goes only to undecided units later in the
#' list. The inclusion probabilities are still `pik` exactly. As with the
#' `"poisson"`, `"sps"` and `"pareto"` methods of [unequal_prob_wor()], a
#' small permanent random number makes selection more likely.
#'
#' Drawing two samples with the same `prn` coordinates them positively
#' (large overlap), and drawing the second with `1 - prn` coordinates them
#' negatively (small overlap). The sample depends on the row order, so
#' coordinated draws must use the same ordering of the frame, for example
#' sorted by a permanent unit identifier.
#'
#' Matei, Smith, Smeets & Klingwort (2023) replace the coordinates in
#' `spread` by a measure of response burden, such as an indicator of
#' highly burdened units or the inclusion probabilities, and draw the
#' surveys with `prn` and `1 - prn`. Spreading on the burden measure keeps
#' burdened units out of the same sample, and negative coordination keeps
#' them out of both samples. This adapted SCP sampling is
#' `balanced_wor(pik, spread = burden, method = "scps", prn = u)`.
#'
#' @return An object of class
#'   `c("balanced", "unequal_prob", "wor", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector of selected unit
#'   indices. When `nrep > 1`, `$sample` is a matrix (n x nrep) for
#'   fixed-size designs, or a list of integer vectors when `fixed_size`
#'   is `FALSE` (e.g., stratified with non-integer per-stratum sizes).
#'
#' @references
#' Bondesson, L. and Thorburn, D. (2008). A list sequential sampling
#'   method suitable for real-time sampling. \emph{Scandinavian Journal of
#'   Statistics}, 35(3), 466-483.
#'   \doi{10.1111/j.1467-9469.2008.00596.x}
#'
#' Deville, J.C. and \enc{Tillé}{Tille}, Y. (1998). Unequal probability
#'   sampling without replacement through a splitting method.
#'   \emph{Biometrika}, 85(1), 89-101.
#'
#' Deville, J.C. and \enc{Tillé}{Tille}, Y. (2004). Efficient balanced sampling: the
#'   cube method. \emph{Biometrika}, 91(4), 893-912.
#'
#' Chauvet, G. and \enc{Tillé}{Tille}, Y. (2006). A fast algorithm for balanced sampling.
#'   \emph{Computational Statistics}, 21(1), 53-62.
#'
#' Chauvet, G. (2009). Stratified balanced sampling. \emph{Survey
#'   Methodology}, 35, 115-119.
#'
#' Friedman, J.H., Bentley, J.L. and Finkel, R.A. (1977). An algorithm for
#'   finding best matches in logarithmic expected time. \emph{ACM
#'   Transactions on Mathematical Software}, 3(3), 209-226.
#'   \doi{10.1145/355744.355745}
#'
#' \enc{Grafström}{Grafstrom}, A., \enc{Lundström}{Lundstrom}, N.L.P.
#'   and Schelin, L. (2012). Spatially balanced sampling through the
#'   pivotal method. \emph{Biometrics}, 68(2), 514-520.
#'   \doi{10.1111/j.1541-0420.2011.01699.x}
#'
#' \enc{Grafström}{Grafstrom}, A. (2012). Spatially correlated Poisson
#'   sampling. \emph{Journal of Statistical Planning and Inference},
#'   142(1), 139-147. \doi{10.1016/j.jspi.2011.07.003}
#'
#' \enc{Grafström}{Grafstrom}, A. and Matei, A. (2018). Coordination of
#'   spatially balanced samples. \emph{Survey Methodology}, 44(2),
#'   215-238.
#'
#' \enc{Grafström}{Grafstrom}, A. and Schelin, L. (2014). How to select
#'   representative samples. \emph{Scandinavian Journal of Statistics},
#'   41(2), 277-290. \doi{10.1111/sjos.12016}
#'
#' Matei, A., Smith, P.A., Smeets, M.J.E. and Klingwort, J. (2023).
#'   Targetted double control of burden in multiple surveys.
#'   \emph{Survey Methodology}, 49(2), 363-384.
#'
#' Tripet, A. and \enc{Tillé}{Tille}, Y. (2026). Balanced sampling with
#'   inequalities: application to category bounding, matrix rounding,
#'   and spread sampling. \emph{Journal of the American Statistical
#'   Association}, 121(553), 796-806.
#'   \doi{10.1080/01621459.2025.2550667}
#'
#' @seealso [unequal_prob_wor()] for unbalanced designs,
#'   [inclusion_prob()] to compute inclusion probabilities from size measures.
#'
#' @examples
#' # Unequal probability balanced sample
#' pik <- c(0.3, 0.6, 0.2, 0.4, 0.5)
#' x <- matrix(c(10, 20, 15, 25, 30))
#' set.seed(1)
#' s <- balanced_wor(pik, aux = x)
#' s$sample
#'
#' # Check balancing: HT estimate of aux totals vs population totals
#' colSums(x[s$sample, , drop = FALSE] / pik[s$sample]) - colSums(x)
#'
#' # Stratified balanced sample
#' N <- 20
#' pik <- rep(0.4, N)
#' x <- matrix(as.double(1:N), ncol = 1)
#' strata <- rep(1:4, each = 5)
#' set.seed(1)
#' s <- balanced_wor(pik, aux = x, strata = strata)
#' s$sample
#'
#' # Controlled selection: bound category counts to the integers
#' # adjacent to their expectation (floor/ceil)
#' pik <- rep(0.5, 12)  # n = 6
#' groups <- rep(c("a", "b", "c"), each = 4)
#' B <- sapply(unique(groups), function(g) as.double(groups == g))
#' S <- colSums(B * pik)  # 2, 2, 2 per category
#' set.seed(1)
#' s <- balanced_wor(
#'   pik,
#'   bounds = list(B = B, lower = floor(S), upper = ceiling(S))
#' )
#' table(groups[s$sample])  # exactly 2 per category
#'
#' # Spatially balanced (well-spread) sample: local pivotal method 2
#' N <- 100
#' coords <- cbind(runif(N), runif(N))
#' pik <- rep(0.1, N)
#' set.seed(1)
#' s <- balanced_wor(pik, spread = coords, method = "lpm2")
#' s$sample
#'
#' # Spatially correlated Poisson sampling uses the same spread contract
#' set.seed(1)
#' s_scps <- balanced_wor(pik, spread = coords, method = "scps")
#' s_scps$sample
#'
#' # Two surveys negatively coordinated with permanent random numbers,
#' # spread on an indicator of highly burdened units (adapted SCP sampling)
#' burdened <- as.double(seq_len(N) %in% c(4, 12, 21, 22, 32, 44))
#' pik1 <- inclusion_prob(runif(N, 1, 5), 10)
#' pik2 <- inclusion_prob(runif(N, 1, 5), 6)
#' u <- runif(N)
#' s1 <- balanced_wor(pik1, spread = burdened, method = "scps", prn = u)
#' s2 <- balanced_wor(pik2, spread = burdened, method = "scps", prn = 1 - u)
#' intersect(s1$sample, s2$sample)
#'
#' @export
balanced_wor <- function(
  pik,
  aux = NULL,
  strata = NULL,
  spread = NULL,
  bounds = NULL,
  method = c("cube", "lpm2", "scps"),
  nrep = 1L,
  prn = NULL,
  ...
) {
  .stop_if_n_arg(sys.call(), "balanced_wor", "sum(pik)")
  if (.is_method_name(method) && is_registered_method(method)) {
    if (!is.null(bounds)) {
      stop(
        "registered balanced methods do not support 'bounds'",
        call. = FALSE
      )
    }
    return(
      .dispatch_registered_balanced(
        pik,
        aux,
        strata,
        spread,
        method,
        nrep,
        prn,
        list(...)
      )
    )
  }
  method <- .match_choice(method, c("cube", "lpm2", "scps"), "method")
  if (...length()) {
    allowed_dots <- switch(
      method,
      cube = c("eps", "condition_aux", "qr_tol"),
      lpm2 = ,
      scps = "eps"
    )
    .check_dots(...length(), ...names(), allowed = allowed_dots)
  }
  nrep <- .check_nrep_prn(
    nrep,
    prn,
    method,
    supports_prn = .method_supports_prn(method, "balanced"),
    supported = "scps"
  )

  if (method %in% c("lpm2", "scps")) {
    if (!is.null(aux)) {
      stop(
        sprintf(
          "method '%s' does not use auxiliary balancing variables; ",
          method
        ),
        "use method = \"cube\" to balance on 'aux'",
        call. = FALSE
      )
    }
    if (!is.null(strata)) {
      stop(
        sprintf("method '%s' does not support 'strata'", method),
        call. = FALSE
      )
    }
    if (!is.null(bounds)) {
      stop(
        "'bounds' is only supported by the built-in \"cube\" method",
        call. = FALSE
      )
    }
    if (is.null(spread)) {
      stop(
        sprintf("method '%s' requires 'spread' ", method),
        "(a matrix of spatial coordinates or other spreading variables)",
        call. = FALSE
      )
    }
    .check_pik(pik, fixed_size = TRUE)
    return(.spatial_wor(pik, spread, method, nrep = nrep, prn = prn, ...))
  }

  if (!is.null(spread)) {
    stop(
      "method 'cube' does not support spatial spreading; ",
      "use method = \"lpm2\", method = \"scps\", or a method registered with ",
      "supports_spread = TRUE",
      call. = FALSE
    )
  }
  if (!is.null(bounds) && !is.null(strata)) {
    stop(
      "'bounds' cannot be combined with 'strata'; ",
      "express within-stratum size constraints through 'bounds' ",
      "(with lower = upper) instead",
      call. = FALSE
    )
  }

  .check_pik(pik, fixed_size = TRUE)

  if (nrep == 1L) {
    .cube_sample(pik, aux = aux, strata = strata, bounds = bounds, ...)
  } else {
    .batch_balanced_wor(pik, aux, strata, bounds, method, nrep, ...)
  }
}

#' Draw a spatial fixed-size WOR sample and wrap the design object.
#'
#' @noRd
.spatial_draw_fns <- list(
  lpm2 = list(
    single = function(pik, spread, eps, nrep) {
      .Call(C_lpm2, pik, spread, eps)
    },
    batch = function(pik, spread, eps, nrep) {
      .Call(C_lpm2_batch, pik, spread, eps, nrep)
    }
  ),
  scps = list(
    single = function(pik, spread, eps, nrep, prn = NULL) {
      .Call(C_scps, pik, spread, eps, prn)
    },
    batch = function(pik, spread, eps, nrep) {
      .Call(C_scps_batch, pik, spread, eps, nrep)
    }
  )
)

#' @noRd
.spatial_wor <- function(
  pik,
  spread,
  method,
  nrep = 1L,
  prn = NULL,
  eps = 1e-10,
  ...
) {
  N <- length(pik)
  eps <- .check_eps(eps)
  .check_cube_eps_classification(pik, eps)
  spread <- .check_cube_aux(spread, N, what = "spread")
  if (ncol(spread) == 0L) {
    stop("'spread' must have at least one column", call. = FALSE)
  }
  n <- as.integer(round(sum(pik)))
  draw <- .spatial_draw_fns[[method]][[
    if (nrep == 1L) "single" else "batch"
  ]]
  draw_args <- list(
    as.double(pik),
    spread,
    as.double(eps),
    as.integer(nrep)
  )
  # prn reaches only methods whose spec allows it (checked upstream), and
  # never a batch draw (prn with nrep > 1 is rejected upstream).
  if (!is.null(prn)) {
    .check_prn(prn, N)
    draw_args$prn <- as.double(prn)
  }
  sample_data <- do.call(draw, draw_args)

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = n,
    N = N,
    method = method,
    fixed_size = TRUE,
    prob_class = "unequal_prob",
    extra_class = "balanced"
  )
}

#' @noRd
.cube_sample <- function(
  pik,
  aux = NULL,
  strata = NULL,
  bounds = NULL,
  eps = 1e-10,
  ...
) {
  N <- length(pik)
  dots <- list(...)
  options <- .check_cube_options(
    pik,
    eps,
    dots[["condition_aux"]],
    dots[["qr_tol"]]
  )
  eps <- options$eps
  condition_aux <- options$condition_aux
  qr_tol <- options$qr_tol

  strata_fixed <- TRUE
  if (is.null(strata)) {
    X <- .build_cube_aux(
      pik,
      aux,
      N,
      prepend_pik = TRUE,
      condition_aux = condition_aux,
      qr_tol = qr_tol
    )
    ineq <- .check_cube_bounds(bounds, pik, N)
    idx <- .Call(
      C_cube,
      as.double(pik),
      X,
      ineq$B,
      ineq$r,
      as.double(eps)
    )
    idx <- .warn_relaxed_bounds(idx)
  } else {
    strata_int <- .check_strata(strata, N)
    strata_fixed <- .check_stratum_sizes(pik, strata_int)
    X <- .build_cube_aux(
      pik,
      aux,
      N,
      prepend_pik = FALSE,
      condition_aux = condition_aux,
      qr_tol = qr_tol,
      strata = strata_int
    )
    idx <- .Call(
      C_cube_stratified,
      as.double(pik),
      X,
      strata_int,
      as.double(eps)
    )
  }

  out <- .new_wor_sample(
    sample = idx,
    pik = pik,
    n = ifelse(strata_fixed, as.integer(round(sum(pik))), sum(pik)),
    N = N,
    method = "cube",
    fixed_size = strata_fixed,
    prob_class = "unequal_prob",
    extra_class = "balanced"
  )
  if (!is.null(bounds)) {
    out$bounds <- bounds
  }
  out
}

#' Reject pik values that the flight-phase tolerance would reclassify.
#'
#' `eps` is the cube flight-phase boundary tolerance: it decides when an
#' updated *working* probability has numerically reached 0 or 1. It must
#' not reinterpret the *input* design. A supplied pik inside (0, eps] or
#' [1 - eps, 1) would be treated as already resolved, silently changing
#' the design and breaking the fixed-size contract, so it is rejected.
#'
#' @noRd
.check_cube_eps_classification <- function(pik, eps) {
  bad <- (pik > 0 & pik <= eps) | (pik >= 1 - eps & pik < 1)
  if (any(bad)) {
    stop(
      sprintf(
        paste0(
          "%d 'pik' value(s) lie within eps = %.2g of 0 or 1 without being ",
          "exactly 0 or exactly 1. 'eps' is the flight-phase boundary ",
          "tolerance, not a trimming rule: units are excluded or selected ",
          "with certainty only when pik is exactly 0 or exactly 1. ",
          "Round those values yourself or use a smaller eps."
        ),
        sum(bad),
        eps
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate options shared by single and batch cube draws.
#'
#' @noRd
.check_cube_options <- function(
  pik,
  eps = NULL,
  condition_aux = NULL,
  qr_tol = NULL
) {
  if (is.null(eps)) {
    eps <- 1e-10
  }
  eps <- .check_eps(eps)
  .check_cube_eps_classification(pik, eps)
  condition_aux <- if (is.null(condition_aux)) {
    FALSE
  } else {
    .check_flag(condition_aux, "condition_aux")
  }
  if (is.null(qr_tol)) {
    qr_tol <- sqrt(.Machine$double.eps)
  } else {
    qr_tol <- .check_number(qr_tol, "qr_tol")
    if (qr_tol < 0) {
      stop("'qr_tol' must be non-negative", call. = FALSE)
    }
  }
  list(eps = eps, condition_aux = condition_aux, qr_tol = qr_tol)
}

#' Build the balancing matrix for the cube C code.
#'
#' For non-stratified: prepends pik as column 1 (sample size constraint).
#' For stratified: the C code auto-adds within-stratum size constraints,
#' so pik is not prepended; returns a 0-column matrix when aux is NULL.
#'
#' @noRd
.build_cube_aux <- function(
  pik,
  aux,
  N,
  prepend_pik,
  condition_aux = FALSE,
  qr_tol = sqrt(.Machine$double.eps),
  strata = NULL
) {
  pik_d <- as.double(pik)

  if (is.null(aux)) {
    if (prepend_pik) {
      return(matrix(pik_d, ncol = 1))
    }
    return(matrix(0, nrow = N, ncol = 0))
  }

  aux <- .check_cube_aux(aux, N)

  if (isTRUE(condition_aux) && ncol(aux) > 0L) {
    aux <- .condition_cube_aux(aux, pik_d, qr_tol = qr_tol, strata = strata)
  }

  if (prepend_pik) {
    cbind(pik_d, aux, deparse.level = 0)
  } else {
    aux
  }
}

#' Validate an auxiliary design matrix (balancing or spreading).
#'
#' Coerces to a double matrix and checks dimensions and finiteness.
#' Shared by the built-in cube path and registered balanced methods;
#' `what` names the argument in error messages ("aux" or "spread").
#'
#' @noRd
.check_cube_aux <- function(aux, N, what = "aux") {
  if (!is.matrix(aux)) {
    aux <- as.matrix(aux)
  }
  if (!is.numeric(aux)) {
    stop(
      sprintf("'%s' must be a numeric vector or matrix", what),
      call. = FALSE
    )
  }
  if (anyNA(aux) || any(!is.finite(aux))) {
    stop(
      sprintf("%s must not contain NA, NaN, or Inf values", what),
      call. = FALSE
    )
  }
  if (nrow(aux) != N) {
    stop(
      sprintf(
        "nrow(%s) = %d does not match length(pik) = %d",
        what,
        nrow(aux),
        N
      ),
      call. = FALSE
    )
  }
  storage.mode(aux) <- "double"
  aux
}

#' Validate inequality bounds and convert to the one-sided form B's <= r
#' used by the C code (Tripet & Tillé 2026).
#'
#' Each two-sided constraint `lower_j <= t(B[, j]) s <= upper_j` becomes
#' up to two one-sided rows: `(+B[, j], upper_j)` and `(-B[, j], -lower_j)`,
#' skipping infinite sides. Always returns a valid (possibly 0-column)
#' system so the C call site does not need to branch.
#'
#' @return `list(B = N x q' double matrix, r = double vector length q')`.
#' @noRd
.check_cube_bounds <- function(bounds, pik, N) {
  if (is.null(bounds)) {
    return(list(B = matrix(numeric(0), nrow = N, ncol = 0), r = double(0)))
  }
  if (!is.list(bounds)) {
    stop(
      "'bounds' must be a list with elements 'B', 'lower', and 'upper'",
      call. = FALSE
    )
  }
  required <- c("B", "lower", "upper")
  missing_el <- setdiff(required, names(bounds))
  if (length(missing_el) > 0L) {
    stop(
      "'bounds' is missing element(s): ",
      paste(missing_el, collapse = ", "),
      call. = FALSE
    )
  }
  extra <- setdiff(names(bounds), required)
  if (length(extra) > 0L) {
    stop(
      "unknown element(s) in 'bounds': ",
      paste(extra, collapse = ", "),
      call. = FALSE
    )
  }

  B <- .check_cube_aux(bounds$B, N, what = "bounds$B")
  q <- ncol(B)
  if (q == 0L) {
    stop("'bounds$B' must have at least one column", call. = FALSE)
  }

  lower <- bounds$lower
  upper <- bounds$upper
  for (side in c("lower", "upper")) {
    v <- if (side == "lower") lower else upper
    if (!is.numeric(v) || !is.null(dim(v)) || length(v) != q) {
      stop(
        sprintf(
          "'bounds$%s' must be a numeric vector of length ncol(bounds$B) = %d",
          side,
          q
        ),
        call. = FALSE
      )
    }
    if (anyNA(v)) {
      stop(
        sprintf("there are missing values in 'bounds$%s'", side),
        call. = FALSE
      )
    }
  }
  if (any(lower > upper)) {
    stop("'bounds$lower' must be <= 'bounds$upper'", call. = FALSE)
  }
  no_side <- !is.finite(lower) & !is.finite(upper)
  if (any(no_side)) {
    stop(
      "constraint(s) ",
      paste(which(no_side), collapse = ", "),
      " in 'bounds' have no finite bound",
      call. = FALSE
    )
  }

  # Feasibility: the walk starts at pik and must start inside the polytope
  v <- as.numeric(crossprod(B, as.double(pik)))
  infeasible <- (is.finite(upper) & v > upper + 1e-8 * (1 + abs(upper))) |
    (is.finite(lower) & v < lower - 1e-8 * (1 + abs(lower)))
  if (any(infeasible)) {
    stop(
      "initial inclusion probabilities violate 'bounds' for constraint(s) ",
      paste(which(infeasible), collapse = ", "),
      "; 'lower <= colSums(B * pik) <= upper' is required",
      call. = FALSE
    )
  }

  up <- which(is.finite(upper))
  lo <- which(is.finite(lower))
  list(
    B = cbind(B[, up, drop = FALSE], -B[, lo, drop = FALSE]),
    r = as.double(c(upper[up], -lower[lo]))
  )
}

#' Warn when the C code had to relax bounds during landing, and strip
#' the signalling attributes from the result.
#' @noRd
.warn_relaxed_bounds <- function(x, nrep = 1L) {
  relaxed <- attr(x, "relaxed", exact = TRUE)
  if (!is.null(relaxed)) {
    warning(
      relaxed,
      " bound(s) could not be met exactly and were relaxed during landing",
      call. = FALSE
    )
    attr(x, "relaxed") <- NULL
  }
  relaxed_reps <- attr(x, "relaxed_reps", exact = TRUE)
  if (!is.null(relaxed_reps)) {
    warning(
      "bounds were relaxed during landing in ",
      relaxed_reps,
      " of ",
      nrep,
      " replicate(s)",
      call. = FALSE
    )
    attr(x, "relaxed_reps") <- NULL
  }
  x
}

#' Condition auxiliary matrix for numerical stability in cube updates.
#'
#' The cube balances `aux / pik`, with `pik / pik = 1` as the size
#' constraint, or the stratum indicators with `strata`. Subtracting `c * pik`
#' from a raw column, within one stratum, subtracts `c` from its `aux / pik`
#' there, a multiple of a size constraint, so the balancing equations are
#' unchanged. The columns are centered that way, within each stratum when
#' `strata` is given, so a column constant within strata in `aux / pik`
#' (such as a region dummy under equal `pik`) is recognized as redundant.
#' They are then scaled and pruned of dependent columns in
#' `aux / pik` space on the units the flight moves (`0 < pik < 1`). The
#' kept columns stay in their original order, which the landing phase
#' relaxes from last to first. A column is dropped as constant when its
#' spread is at most `qr_tol` times its size in that space, so the result
#' does not depend on the scale of a column.
#'
#' Centering the raw columns instead added `-mu / pik` to each constraint,
#' which is outside the size constraint when `pik` varies, and dropped a
#' constant raw column that is not constant after division by `pik`.
#'
#' @param aux Numeric matrix (N x p).
#' @param pik Numeric vector of inclusion probabilities (length N).
#' @param qr_tol Tolerance for QR rank detection.
#' @param strata Dense integer stratum labels, or `NULL`.
#' @return Conditioned matrix, possibly with fewer columns.
#' @noRd
.condition_cube_aux <- function(aux, pik, qr_tol = sqrt(.Machine$double.eps),
                                strata = NULL) {
  N <- nrow(aux)
  if (N == 0L || ncol(aux) == 0L) {
    return(aux)
  }
  active <- pik > 0 & pik < 1
  if (!any(active)) {
    return(matrix(0, nrow = N, ncol = 0))
  }

  # The constant comes from the units the flight moves, so a column that
  # is constant there has no spread left after centering.
  flight <- function(x) x[active, , drop = FALSE] / pik[active]
  w <- pik[active] / sum(pik[active])
  rms <- sqrt(colSums(w * flight(aux)^2))
  if (is.null(strata)) {
    aux_c <- aux - outer(pik, colSums(aux[active, , drop = FALSE]) /
      sum(pik[active]))
  } else {
    # One multiple of pik per stratum, from its moving units. A stratum
    # with none keeps its values, which the flight never reads.
    g <- factor(strata)
    num <- rowsum(aux[active, , drop = FALSE], g[active])
    den <- as.numeric(rowsum(pik[active], g[active]))
    ratio <- matrix(0, nlevels(g), ncol(aux))
    ratio[match(rownames(num), levels(g)), ] <- num / den
    aux_c <- aux - pik * ratio[as.integer(g), , drop = FALSE]
  }
  sdw <- sqrt(pmax(colSums(w * flight(aux_c)^2), 0))
  # Relative to the column's own size, so rescaling a column cannot drop
  # it: an absolute cutoff dropped 1e-8 * aux, an independent constraint.
  keep_scale <- sdw > qr_tol * rms
  if (!any(keep_scale)) {
    return(matrix(0, nrow = N, ncol = 0))
  }
  aux_c <- sweep(
    aux_c[, keep_scale, drop = FALSE],
    2,
    sdw[keep_scale],
    "/",
    check.margin = FALSE
  )

  if (ncol(aux_c) <= 1L) {
    return(aux_c)
  }
  # LINPACK path: rank-revealing QR that honors `tol` and moves only
  # near-dependent columns to the end. The LAPACK path ignores `tol` and
  # always reports full column rank, so it cannot detect dependent
  # balancing constraints.
  q <- qr(flight(aux_c), tol = qr_tol)
  r <- q$rank
  if (r <= 0L) {
    return(matrix(0, nrow = N, ncol = 0))
  }
  aux_c[, sort(q$pivot[seq_len(r)]), drop = FALSE]
}

#' @noRd
.check_strata <- function(strata, N) {
  if (!is.numeric(strata) || !is.null(dim(strata))) {
    stop(
      "'strata' must be a numeric vector of positive integers",
      call. = FALSE
    )
  }
  if (length(strata) != N) {
    stop(
      sprintf("length(strata) = %d does not match N = %d", length(strata), N),
      call. = FALSE
    )
  }
  if (anyNA(strata)) {
    stop("there are missing values in 'strata'", call. = FALSE)
  }
  if (
    any(!is.finite(strata)) ||
      any(strata < 1) ||
      any(strata != floor(strata)) ||
      any(strata > .Machine$integer.max)
  ) {
    stop("'strata' values must be positive integers", call. = FALSE)
  }
  strata <- as.integer(strata)
  # Remap to dense 1:H so the C code doesn't over-allocate for sparse labels
  as.integer(factor(strata))
}

#' Check per-stratum sum(pik). Warns and returns FALSE if any is non-integer.
#'
#' A stratum sum counts as whole within the tolerance `.check_whole_sum()`
#' applies to the whole of `pik`, scaled by the stratum's own count and
#' sum. A fixed `1e-4` took a stratum summing to 5e-5 as whole, and the
#' stratified batch kernel then failed on a draw of the wrong size.
#'
#' @return TRUE if all per-stratum sums are close to an integer, FALSE otherwise.
#' @noRd
.check_stratum_sizes <- function(pik, strata, tol = NULL) {
  stratum_sums <- tapply(pik, strata, sum)
  if (is.null(tol)) {
    counts <- tapply(pik, strata, length)
    tol <- pmax(1e-10, 64 * .Machine$double.eps * (counts + abs(stratum_sums)))
  }
  not_int <- abs(stratum_sums - round(stratum_sums)) > tol
  if (any(not_int)) {
    bad <- names(stratum_sums)[not_int]
    warning(
      "per-stratum sum(pik) is not close to an integer for stratum ",
      paste(bad, collapse = ", "),
      "; within-stratum sample sizes will not be exact",
      call. = FALSE
    )
    return(FALSE)
  }
  TRUE
}

#' @noRd
.batch_balanced_wor <- function(pik, aux, strata, bounds, method, nrep, ...) {
  N <- length(pik)
  n <- sum(pik)
  dots <- list(...)
  options <- .check_cube_options(
    pik,
    dots[["eps"]],
    dots[["condition_aux"]],
    dots[["qr_tol"]]
  )
  eps <- options$eps
  condition_aux <- options$condition_aux
  qr_tol <- options$qr_tol
  pik_d <- as.double(pik)

  if (is.null(strata)) {
    X <- .build_cube_aux(
      pik,
      aux,
      N,
      prepend_pik = TRUE,
      condition_aux = condition_aux,
      qr_tol = qr_tol
    )
    ineq <- .check_cube_bounds(bounds, pik, N)
    sample_data <- .Call(
      C_cube_batch,
      pik_d,
      X,
      ineq$B,
      ineq$r,
      as.double(eps),
      as.integer(nrep)
    )
    sample_data <- .warn_relaxed_bounds(sample_data, nrep = nrep)
    fixed_size <- TRUE
  } else {
    strata_int <- .check_strata(strata, N)
    fixed_size <- .check_stratum_sizes(pik, strata_int)
    X <- .build_cube_aux(
      pik,
      aux,
      N,
      prepend_pik = FALSE,
      condition_aux = condition_aux,
      qr_tol = qr_tol,
      strata = strata_int
    )

    if (fixed_size) {
      sample_data <- .Call(
        C_cube_stratified_batch,
        pik_d,
        X,
        strata_int,
        as.double(eps),
        as.integer(nrep)
      )
    } else {
      sample_data <- vector("list", nrep)
      for (i in seq_len(nrep)) {
        sample_data[[i]] <- .Call(
          C_cube_stratified,
          pik_d,
          X,
          strata_int,
          as.double(eps)
        )
      }
    }
  }

  out <- .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = if (fixed_size) as.integer(round(n)) else n,
    N = N,
    method = "cube",
    fixed_size = fixed_size,
    prob_class = "unequal_prob",
    extra_class = "balanced"
  )
  if (!is.null(bounds)) {
    out$bounds <- bounds
  }
  out
}
