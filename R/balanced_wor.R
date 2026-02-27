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
#' @param ... Additional arguments passed to methods (e.g., `eps` for
#'   boundary tolerance).
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
#'   indices. When `nrep > 1`, `$sample` is a matrix (n x nrep).
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
  method = c("cube"),
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

  strata_fixed <- TRUE
  if (is.null(strata)) {
    X <- .build_cube_aux(pik, aux, N, prepend_pik = TRUE)
    idx <- .Call(C_cube, as.double(pik), X, as.double(eps))
  } else {
    strata_int <- .check_strata(strata, N)
    strata_fixed <- .check_stratum_sizes(pik, strata_int)
    X <- .build_cube_aux(pik, aux, N, prepend_pik = FALSE)
    idx <- .Call(C_cube_stratified, as.double(pik), X, strata_int, as.double(eps))
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
#' For non-stratified: prepends pik as column 0 (sample size constraint).
#' For stratified: the C code auto-adds within-stratum size constraints,
#' so pik is not prepended.
#'
#' @noRd
.build_cube_aux <- function(pik, aux, N, prepend_pik) {
  pik_d <- as.double(pik)

  if (is.null(aux)) {
    return(matrix(pik_d, ncol = 1))
  }

  if (!is.matrix(aux)) aux <- as.matrix(aux)
  if (nrow(aux) != N) {
    stop(
      sprintf("nrow(aux) = %d does not match length(pik) = %d", nrow(aux), N),
      call. = FALSE
    )
  }
  storage.mode(aux) <- "double"

  if (prepend_pik) {
    cbind(pik_d, aux, deparse.level = 0)
  } else {
    aux
  }
}

#' @noRd
.check_strata <- function(strata, N) {
  if (length(strata) != N) {
    stop(
      sprintf("length(strata) = %d does not match N = %d",
              length(strata), N),
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

  first <- .cube_sample(pik, aux = aux, strata = strata, ...)
  fixed_size <- first$fixed_size

  if (fixed_size) {
    n_int <- as.integer(round(n))
    mat <- matrix(0L, n_int, nrep)
    mat[, 1L] <- first$sample
    for (i in seq_len(nrep - 1L) + 1L) {
      mat[, i] <- .cube_sample(pik, aux = aux, strata = strata, ...)$sample
    }
    sample_data <- mat
  } else {
    sample_data <- vector("list", nrep)
    sample_data[[1L]] <- first$sample
    for (i in seq_len(nrep - 1L) + 1L) {
      sample_data[[i]] <- .cube_sample(pik, aux = aux, strata = strata, ...)$sample
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
