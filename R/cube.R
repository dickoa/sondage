#' Cube Method for Balanced Sampling
#'
#' Selects balanced samples with prescribed inclusion probabilities
#' using the fast flight Cube Method.
#'
#' @param prob Numeric vector of inclusion probabilities (length N).
#' @param x Matrix of balancing variables (N x p). If \code{prob} is included
#'   as a balancing variable, the sample size will be fixed.
#' @param strata Optional integer vector of stratum identifiers (length N).
#'   When provided, uses stratified cube method guaranteeing exact sample
#'   size per stratum.
#' @param nrep Number of samples to draw (default: 1).
#' @param eps Numerical tolerance (default: 1e-12).
#'
#' @return If \code{nrep = 1}, an integer vector of selected indices.
#'   If \code{nrep > 1}, an integer matrix with n rows and \code{nrep} columns,
#'   where each column contains the indices for one replicate.
#'
#' @details
#' The cube method (Deville and Tillé, 2004) selects balanced samples that
#' satisfy (approximately) the balancing equations:
#' \deqn{\sum_{k \in S} x_k / \pi_k = \sum_{k \in U} x_k}
#'
#' This implementation uses the fast flight phase algorithm of Chauvet and
#' Tillé (2006), which processes only p+1 units at a time, achieving O(N)
#' complexity.
#'
#' The algorithm has two phases:
#' \itemize{
#'   \item Flight phase: Iteratively rounds probabilities to 0 or 1 while
#'     maintaining exact balance, until at most p non-integer probabilities
#'     remain.
#'   \item Landing phase: Handles the remaining units, relaxing the balancing
#'     constraints minimally.
#' }
#'
#' When \code{strata} is provided, the stratified cube method is used:
#' \enumerate{
#'   \item Flight per stratum: balance within each stratum
#'   \item Global flight: balance across strata with stratum indicators
#'   \item Landing per stratum: separate landing ensures exact n_h
#' }
#' This guarantees exact sample sizes per stratum.
#'
#' For repeated sampling (simulations), use the \code{nrep} parameter instead
#' of a loop for better performance.
#'
#' @references
#' Deville, J.C. and Tillé, Y. (2004). Efficient balanced sampling: the cube
#' method. \emph{Biometrika}, 91(4), 893-912.
#'
#' Chauvet, G. and Tillé, Y. (2006). A fast algorithm for balanced sampling.
#' \emph{Computational Statistics}, 21(1), 53-62.
#'
#' @examples
#' # Basic balanced sampling
#' set.seed(123)
#' N <- 100
#' n <- 10
#' prob <- rep(n / N, N)
#' x <- cbind(prob, matrix(runif(N * 2), ncol = 2))
#' s <- cube(prob, x)
#' print(s)
#'
#' # Stratified balanced sampling
#' set.seed(123)
#' N <- 200
#' strata <- rep(1:4, each = 50)
#' prob <- rep(10 / 50, N)  # 10 per stratum
#' x <- cbind(prob, runif(N))
#' s <- cube(prob, x, strata = strata)
#' table(strata[s])  # Exact 10 per stratum
#'
#' #' # Batch sampling for simulations
#' set.seed(1)
#' samples <- cube(prob, x, nrep = 100)
#' emp_prob <- tabulate(as.vector(samples), nbins = N) / 100
#'
#' @export
cube <- function(prob, x, strata = NULL, nrep = 1L, eps = 1e-12) {
  if (!is.numeric(prob)) {
    stop("'prob' must be a numeric vector")
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  N <- length(prob)
  if (nrow(x) != N) {
    stop("'prob' and 'x' must have the same number of rows")
  }
  if (any(prob < 0) || any(prob > 1)) {
    stop("'prob' must be in [0, 1]")
  }
  if (eps <= 0) {
    stop("'eps' must be positive")
  }
  nrep <- as.integer(nrep)
  if (nrep < 1L) {
    stop("'nrep' must be at least 1")
  }

  if (!is.null(strata)) {
    if (length(strata) != N) {
      stop("'strata' must have the same length as 'prob'")
    }
    strata <- as.integer(as.factor(strata)) # Convert to 1..H integers

    if (nrep == 1L) {
      .Call(C_cube_stratified, as.double(prob), x, strata, as.double(eps))
    } else {
      n <- round(sum(prob))
      result <- matrix(NA_integer_, nrow = n, ncol = nrep)
      for (i in seq_len(nrep)) {
        s <- .Call(
          C_cube_stratified,
          as.double(prob),
          x,
          strata,
          as.double(eps)
        )
        result[seq_along(s), i] <- s
      }
      result
    }
  } else {
    if (nrep == 1L) {
      .Call(C_cube, as.double(prob), x, as.double(eps))
    } else {
      .Call(C_cube_batch, as.double(prob), x, as.double(eps), nrep)
    }
  }
}
