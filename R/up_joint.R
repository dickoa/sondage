#' Joint Inclusion Probabilities for Poisson Sampling
#'
#' Computes joint inclusion probabilities for Poisson sampling.
#' Since selections are independent, \eqn{\pi_{ij} = \pi_i \times \pi_j}.
#'
#' @param pik numeric vector of first-order inclusion probabilities.
#'
#' @return A symmetric N×N matrix of joint inclusion probabilities.
#'
#' @details
#' For Poisson sampling, units are selected independently, so:
#' \deqn{\pi_{ij} = \pi_i \pi_j \quad \text{for } i \neq j}
#'
#' @examples
#' pik <- c(0.2, 0.5, 0.8)
#' pikl <- up_poisson_jip(pik)
#'
#' @export
up_poisson_jip <- function(pik) {
  check_pik(pik)
  pikl <- outer(pik, pik)
  diag(pikl) <- pik
  pikl
}

#' Joint Inclusion Probabilities for Maximum Entropy (CPS) Sampling
#'
#' Computes exact joint inclusion probabilities for Conditional Poisson
#' Sampling (CPS), also known as Maximum Entropy sampling.
#'
#' @param pik numeric vector of first-order inclusion probabilities.
#' @param eps tolerance for boundary detection (default 1e-6).
#'
#' @return A symmetric N×N matrix of joint inclusion probabilities.
#'
#' @details
#' Uses Aires' formula from Tillé (2006) Expression 5.20:
#' \deqn{\pi_{k\ell} = \frac{r_\ell \pi_k - r_k \pi_\ell}{r_\ell - r_k}}
#' where \eqn{r_k = \exp(\lambda_k)}.
#'
#' @references
#' Aires, N. (1999). Algorithms to find exact inclusion probabilities for
#' conditional Poisson sampling. \emph{Methodology and Computing in Applied
#' Probability}, 4, 457-469.
#'
#' Tillé, Y. (2006). \emph{Sampling Algorithms}. Springer.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' pikl <- up_maxent_jip(pik)
#'
#' @export
up_maxent_jip <- function(pik, eps = 1e-6) {
  check_pik(pik)
  .Call(C_up_maxent_jip, as.double(pik), as.double(eps))
}

#' Joint Inclusion Probabilities for Systematic Sampling
#'
#' Computes exact joint inclusion probabilities for systematic sampling
#' with unequal probabilities.
#'
#' @param pik numeric vector of first-order inclusion probabilities.
#' @param eps tolerance for boundary detection (default 1e-6).
#'
#' @return A symmetric N×N matrix of joint inclusion probabilities.
#'
#' @note Systematic sampling is NOT a high-entropy design. Many joint
#' probabilities will be exactly zero.
#'
#' @references
#' Tillé, Y. (2006). \emph{Sampling Algorithms}. Springer.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' pikl <- up_systematic_jip(pik)
#'
#' @export
up_systematic_jip <- function(pik, eps = 1e-6) {
  check_pik(pik)
  .Call(C_up_systematic_jip, as.double(pik), as.double(eps))
}

#' Joint Inclusion Probabilities using Brewer's Approximation
#'
#' Computes approximate joint inclusion probabilities using equation (18)
#' from Brewer & Donadio (2003).
#'
#' @param pik numeric vector of first-order inclusion probabilities.
#'
#' @return A symmetric N×N matrix of approximate joint probabilities.
#'
#' @details
#' Uses the approximation:
#' \deqn{c_i = \frac{n-1}{n - \frac{2n-1}{n-1}\pi_i + \frac{\sum \pi_k^2}{n-1}}}
#' \deqn{\tilde{\pi}_{ij} = \pi_i \pi_j \frac{c_i + c_j}{2}}
#'
#' This has the best model-assisted properties among Brewer's approximations.
#'
#' @references
#' Brewer, K.R.W. and Donadio, M.E. (2003). The High Entropy Variance of the
#' Horvitz-Thompson Estimator. \emph{Survey Methodology}, 29(2), 189-196.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' pikl <- up_brewer_jip(pik)
#'
#' @export
up_brewer_jip <- function(pik) {
  check_pik(pik)
  .Call(C_up_brewer_jip, as.double(pik))
}

#' Pairwise Expectations for Chromy Sampling
#'
#' Estimates pairwise hit expectations \eqn{E(n_i n_j)} for Chromy's sequential
#' PPS sampling via Monte Carlo simulation.
#'
#' @param x Numeric vector of positive size measures.
#' @param n Sample size.
#' @param nsim Number of simulations (default 10000).
#'
#' @return A symmetric N x N matrix of pairwise expectations. Entry (i, j)
#'   is \eqn{E(n_i n_j)}, the expected product of hit counts. Diagonal
#'   entries are \eqn{E(n_k^2)}.
#'
#' @details
#' Chromy's method is a minimum replacement (PMR) design where units with
#' large size measures can be selected multiple times. The appropriate
#' variance estimator uses pairwise expectations rather than joint inclusion
#' probabilities.
#'
#' Chromy (2009) gives the generalized Yates-Grundy variance:
#' \deqn{V(\hat{T}) = \frac{1}{2} \sum_{i \neq j} \{E(n_i)E(n_j) - E(n_i n_j)\}
#'   \left(\frac{y_i}{E(n_i)} - \frac{y_j}{E(n_j)}\right)^2}
#'
#' where \eqn{E(n_k) = n x_k / \sum x} is exact.
#'
#' In the without-replacement case (all \eqn{E(n_k) < 1}), this reduces to
#' the standard Sen-Yates-Grundy formula with \eqn{E(n_i n_j) = \pi_{ij}}.
#'
#' @references
#' Chromy, J.R. (2009). Some Generalizations of the Horvitz-Thompson Estimator.
#' \emph{Proceedings of the Survey Research Methods Section, ASA}, 216-227.
#'
#' Chauvet, G. (2019). Properties of Chromy's sampling procedure.
#' \emph{arXiv:1912.10896}.
#'
#' @seealso [up_chromy()] for sampling, [up_brewer_jip()] for joint
#'   inclusion probabilities under high-entropy WOR designs.
#'
#' @examples
#' x <- c(10, 20, 15, 25, 30)
#' pairexp <- up_chromy_pairexp(x, n = 3, nsim = 5000)
#'
#' # Expected hits (exact)
#' En <- 3 * x / sum(x)
#'
#' # In WOR case: diagonal approximates E(n_k)
#' diag(pairexp)  # Should be close to En
#'
#' # Covariance structure for variance estimation
#' En_outer <- outer(En, En)
#' Cov_nn <- pairexp - En_outer  # E(n_i n_j) - E(n_i)E(n_j)
#'
#' @export
up_chromy_pairexp <- function(x, n, nsim = 10000L) {
  check_mos(x)

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("n must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    stop("nsim must be a positive integer", call. = FALSE)
  }

  .Call(C_up_chromy_pairexp, as.double(x), as.integer(n), as.integer(nsim))
}
