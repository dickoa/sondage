#' Maximum Entropy Sampling (Conditional Poisson Sampling)
#'
#' Draws samples using the maximum entropy design, also known as Conditional
#' Poisson Sampling (CPS). This is the unique design that maximizes entropy
#' subject to fixed inclusion probabilities.
#'
#' @param pik A numeric vector of inclusion probabilities. The sum should be
#'   an integer representing the desired sample size.
#' @param nrep Number of sample replicates to draw. Default is 1.
#' @param eps A small threshold value for boundary cases. Default is 1e-06.
#'
#' @return If \code{nrep = 1}, an integer vector of selected indices.
#'   If \code{nrep > 1}, an integer matrix with n rows and \code{nrep} columns,
#'   where each column contains the indices for one replicate.
#'
#' @details
#' Maximum entropy sampling has several desirable properties:
#' \itemize{
#'   \item Fixed sample size: exactly `round(sum(pik))` units selected
#'   \item Exact inclusion probabilities: \eqn{E(I_k) = \pi_k}
#'   \item All joint inclusion probabilities are positive: \eqn{\pi_{kl} > 0}
#'   \item Maximum entropy among all designs with fixed \eqn{\pi_k}
#' }
#'
#' The implementation uses the sequential algorithm of Chen, Dempster, and Liu
#' (1994) as described in Tillé (2006). The q-values are computed on-the-fly
#' to reduce memory usage from O(N*n) for the full q-matrix to O(N) working
#' arrays.
#'
#' For repeated sampling (simulations), use the `nrep` parameter instead of
#' a loop for much better performance. The design is computed once and reused
#' for all replicates.
#'
#' @references
#' Chen, S. X., Dempster, A. P., & Liu, J. S. (1994). Weighted finite population
#'   sampling to maximize entropy. \emph{Biometrika}, 81(3), 457-469.
#'
#' Tillé, Y. (2006). \emph{Sampling Algorithms}. Springer Series in Statistics.
#'   Chapter 6.
#'
#' @seealso [up_brewer()] for Brewer's method, [up_systematic()] for systematic PPS
#'
#' @examples
#' pik <- c(0.2, 0.4, 0.6, 0.8)  # sum = 2
#'
#' # Single sample
#' set.seed(123)
#' idx <- up_maxent(pik)
#' idx
#'
#' # Multiple replicates for simulation
#' samples <- up_maxent(pik, nrep = 1000)
#' dim(samples)  # 2 x 1000
#'
#' # Verify inclusion probabilities
#' rowMeans(apply(samples, 2, function(s) 1:4 %in% s))  # close to pik
#'
#' @export
up_maxent <- function(pik, nrep = 1L, eps = 1e-06) {
  if (any(is.na(pik))) {
    stop("there are missing values in the pik vector", call. = FALSE)
  }
  if (!is.numeric(pik)) {
    stop("pik must be a numeric vector", call. = FALSE)
  }
  if (length(pik) == 0) {
    stop("pik vector is empty", call. = FALSE)
  }
  if (!is.numeric(nrep) || length(nrep) != 1 || nrep < 1) {
    stop("nrep must be a positive integer", call. = FALSE)
  }

  if (any(pik < 0 | pik > 1)) {
    stop("inclusion probabilities must be between 0 and 1", call. = FALSE)
  }

  nrep <- as.integer(nrep)

  if (nrep == 1L) {
    .Call(C_maxent_single, as.double(pik), as.double(eps))
  } else {
    design <- .Call(C_maxent_design, as.double(pik), as.double(eps))
    .Call(C_maxent_draw_batch, design, nrep)
  }
}
