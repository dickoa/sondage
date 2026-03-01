#' Balanced Sampling Without Replacement
#'
#' Draws a balanced sample using the cube method (Deville & Tille, 2004).
#' A balanced sample satisfies (approximately) the balancing equations
#' \eqn{\sum_{k \in S} x_k / \pi_k \approx \sum_{k \in U} x_k} for each
#' auxiliary variable \eqn{x}.
#'
#' @param pik A numeric vector of inclusion probabilities (length N).
#'   `sum(pik)` must be close to an integer.
#' @param aux An optional numeric matrix (N x p) of auxiliary balancing
#'   variables. Each column defines a balancing constraint. The sample
#'   size constraint is always included automatically; `aux` specifies
#'   additional variables to balance on.
#'   When `NULL`, only the sample size is balanced (equivalent to an
#'   unbalanced fixed-size design).
#' @param strata An optional integer vector (length N) of stratum
#'   indicators (positive integers). When provided, uses the stratified
#'   cube method (Chauvet & Tille, 2006; Chauvet, 2009) which preserves
#'   within-stratum sample sizes while balancing on `aux`. Exact
#'   preservation requires `sum(pik)` within each stratum to be close
#'   to an integer; otherwise sizes will vary around the target and a
#'   warning is issued and the design is marked as random-size
#'   (`fixed_size = FALSE`), so batch replicates are returned as a list.
#' @param method The sampling method. Currently only `"cube"`.
#' @param nrep Number of replicate samples (default 1). When `nrep > 1`,
#'   `$sample` holds a matrix (n x nrep) for fixed-size designs, or a
#'   list of integer vectors when within-stratum sizes are not exact.
#' @param ... Additional arguments passed to methods:
#'   \describe{
#'     \item{`eps`}{Boundary tolerance for deciding 0/1 (default `1e-10`).}
#'     \item{`condition_aux`}{Logical; if `TRUE`, pre-conditions `aux` by
#'       weighted centering/scaling and QR-pivot rank pruning to improve
#'       numerical stability with ill-conditioned or collinear auxiliary
#'       variables (default `FALSE`).}
#'     \item{`qr_tol`}{Tolerance for QR rank detection when
#'       `condition_aux = TRUE` (default `sqrt(.Machine$double.eps)`).}
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
#' approximation (Brewer & Donadio, 2003), which is appropriate since
#' the cube produces a near-maximum-entropy design.
#'
#' @return An object of class `c("unequal_prob", "wor", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector of selected unit
#'   indices. When `nrep > 1`, `$sample` is a matrix (n x nrep) for
#'   fixed-size designs, or a list of integer vectors when `fixed_size`
#'   is `FALSE` (e.g., stratified with non-integer per-stratum sizes).
#'
#' @references
#' Deville, J.C. and Tille, Y. (2004). Efficient balanced sampling: the
#'   cube method. \emph{Biometrika}, 91(4), 893-912.
#'
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm for balanced sampling.
#'   \emph{Computational Statistics}, 21(1), 53-62.
#'
#' Chauvet, G. (2009). Stratified balanced sampling. \emph{Survey
#'   Methodology}, 35, 115-119.
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
#' @export
balanced_wor <- function(
  pik,
  aux = NULL,
  strata = NULL,
  method = "cube",
  nrep = 1L,
  ...
) {
  method <- match.arg(method)
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  check_pik(pik, fixed_size = TRUE)

  if (nrep == 1L) {
    .cube_sample(pik, aux = aux, strata = strata, ...)
  } else {
    .batch_balanced_wor(pik, aux, strata, method, nrep, ...)
  }
}

#' @noRd
.cube_sample <- function(pik, aux = NULL, strata = NULL, eps = 1e-10, ...) {
  N <- length(pik)
  dots <- list(...)
  condition_aux <- isTRUE(dots[["condition_aux"]])
  qr_tol <- dots[["qr_tol"]]
  if (is.null(qr_tol)) qr_tol <- sqrt(.Machine$double.eps)
  eps <- check_eps(eps)

  strata_fixed <- TRUE
  if (is.null(strata)) {
    X <- .build_cube_aux(pik, aux, N, prepend_pik = TRUE,
                         condition_aux = condition_aux, qr_tol = qr_tol)
    idx <- .Call(C_cube, as.double(pik), X, as.double(eps))
  } else {
    strata_int <- .check_strata(strata, N)
    strata_fixed <- .check_stratum_sizes(pik, strata_int)
    X <- .build_cube_aux(pik, aux, N, prepend_pik = FALSE,
                         condition_aux = condition_aux, qr_tol = qr_tol)
    idx <- .Call(
      C_cube_stratified,
      as.double(pik),
      X,
      strata_int,
      as.double(eps)
    )
  }

  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = if (strata_fixed) as.integer(round(sum(pik))) else sum(pik),
    N = N,
    method = "cube",
    fixed_size = strata_fixed,
    prob_class = "unequal_prob"
  )
}

#' Build the balancing matrix for the cube C code.
#'
#' For non-stratified: prepends pik as column 1 (sample size constraint).
#' For stratified: the C code auto-adds within-stratum size constraints,
#' so pik is not prepended; returns a 0-column matrix when aux is NULL.
#'
#' @noRd
.build_cube_aux <- function(pik, aux, N, prepend_pik,
                            condition_aux = FALSE,
                            qr_tol = sqrt(.Machine$double.eps)) {
  pik_d <- as.double(pik)

  if (is.null(aux)) {
    if (prepend_pik) {
      return(matrix(pik_d, ncol = 1))
    }
    return(matrix(0, nrow = N, ncol = 0))
  }

  if (!is.matrix(aux)) {
    aux <- as.matrix(aux)
  }
  if (anyNA(aux) || any(!is.finite(aux))) {
    stop("aux must not contain NA, NaN, or Inf values", call. = FALSE)
  }
  if (nrow(aux) != N) {
    stop(
      sprintf("nrow(aux) = %d does not match length(pik) = %d", nrow(aux), N),
      call. = FALSE
    )
  }
  storage.mode(aux) <- "double"

  if (isTRUE(condition_aux) && ncol(aux) > 0L) {
    aux <- .condition_cube_aux(aux, pik_d, qr_tol = qr_tol)
  }

  if (prepend_pik) {
    cbind(pik_d, aux, deparse.level = 0)
  } else {
    aux
  }
}

#' Condition auxiliary matrix for numerical stability in cube updates.
#'
#' Applies weighted centering/scaling (weights proportional to pik), then
#' QR with column pivoting to drop near-dependent columns while preserving
#' the effective constraint span.
#'
#' @param aux Numeric matrix (N x p).
#' @param pik Numeric vector of inclusion probabilities (length N).
#' @param qr_tol Tolerance for QR rank detection.
#' @return Conditioned matrix, possibly with fewer columns.
#' @noRd
.condition_cube_aux <- function(aux, pik, qr_tol = sqrt(.Machine$double.eps)) {
  N <- nrow(aux)
  if (N == 0L || ncol(aux) == 0L) return(aux)

  # Weighted centering
  w <- pik / sum(pik)
  mu <- as.numeric(crossprod(w, aux))
  aux_cs <- sweep(aux, 2, mu, "-", check.margin = FALSE)

  # Weighted scaling — drop zero-variance columns
  sdw <- sqrt(pmax(as.numeric(crossprod(w, aux_cs^2)), 0))
  keep_scale <- sdw > qr_tol
  if (!any(keep_scale)) {
    return(matrix(0, nrow = N, ncol = 0))
  }
  aux_cs <- aux_cs[, keep_scale, drop = FALSE]
  sdw <- sdw[keep_scale]
  aux_cs <- sweep(aux_cs, 2, sdw, "/", check.margin = FALSE)

  # QR with column pivoting to prune linearly dependent columns
  if (ncol(aux_cs) <= 1L) return(aux_cs)
  q <- qr(aux_cs, tol = qr_tol, LAPACK = TRUE)
  r <- q$rank
  if (r <= 0L) {
    return(matrix(0, nrow = N, ncol = 0))
  }
  aux_cs[, q$pivot[seq_len(r)], drop = FALSE]
}

#' @noRd
.check_strata <- function(strata, N) {
  if (length(strata) != N) {
    stop(
      sprintf("length(strata) = %d does not match N = %d", length(strata), N),
      call. = FALSE
    )
  }
  strata <- as.integer(strata)
  if (anyNA(strata)) {
    stop("there are missing values in 'strata'", call. = FALSE)
  }
  if (any(strata < 1L)) {
    stop("'strata' values must be positive integers", call. = FALSE)
  }
  # Remap to dense 1:H so the C code doesn't over-allocate for sparse labels
  as.integer(factor(strata))
}

#' Check per-stratum sum(pik). Warns and returns FALSE if any is non-integer.
#' @return TRUE if all per-stratum sums are close to an integer, FALSE otherwise.
#' @noRd
.check_stratum_sizes <- function(pik, strata, tol = 1e-4) {
  stratum_sums <- tapply(pik, strata, sum)
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
.batch_balanced_wor <- function(pik, aux, strata, method, nrep, ...) {
  N <- length(pik)
  n <- sum(pik)
  dots <- list(...)
  eps <- dots[["eps"]]
  if (is.null(eps)) eps <- 1e-10
  eps <- check_eps(eps)
  condition_aux <- isTRUE(dots[["condition_aux"]])
  qr_tol <- dots[["qr_tol"]]
  if (is.null(qr_tol)) qr_tol <- sqrt(.Machine$double.eps)
  pik_d <- as.double(pik)

  if (is.null(strata)) {
    X <- .build_cube_aux(pik, aux, N, prepend_pik = TRUE,
                         condition_aux = condition_aux, qr_tol = qr_tol)
    sample_data <- .Call(
      C_cube_batch,
      pik_d,
      X,
      as.double(eps),
      as.integer(nrep)
    )
    fixed_size <- TRUE
  } else {
    strata_int <- .check_strata(strata, N)
    fixed_size <- .check_stratum_sizes(pik, strata_int)
    X <- .build_cube_aux(pik, aux, N, prepend_pik = FALSE,
                         condition_aux = condition_aux, qr_tol = qr_tol)

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

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = if (fixed_size) as.integer(round(n)) else n,
    N = N,
    method = "cube",
    fixed_size = fixed_size,
    prob_class = "unequal_prob"
  )
}
