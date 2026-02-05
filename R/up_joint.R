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
#' pikl <- up_poisson_joint(pik)
#'
#' @export
up_poisson_joint <- function(pik) {
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
#' pikl <- up_maxent_joint(pik)
#'
#' @export
up_maxent_joint <- function(pik, eps = 1e-6) {
  check_pik(pik)
  .Call(C_up_maxent_joint, as.double(pik), as.double(eps))
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
#' pikl <- up_systematic_joint(pik)
#'
#' @export
up_systematic_joint <- function(pik, eps = 1e-6) {
  check_pik(pik)
  .Call(C_up_systematic_joint, as.double(pik), as.double(eps))
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
#' pikl <- up_brewer_joint(pik)
#'
#' @export
up_brewer_joint <- function(pik) {
  check_pik(pik)
  .Call(C_up_brewer_joint, as.double(pik))
}

#' Joint Inclusion Probabilities for Chromy Sampling
#'
#' Estimates joint inclusion probabilities for Chromy's sequential PPS
#' sampling via Monte Carlo simulation.
#'
#' @param x Numeric vector of positive size measures.
#' @param n Sample size.
#' @param nsim Number of simulations (default 10000).
#'
#' @return A symmetric N×N matrix of joint inclusion probabilities.
#'   Diagonal entries are first-order probabilities π_k.
#'
#' @details
#' Chromy's method with randomization has all π_kl > 0, enabling unbiased
#' variance estimation. Chauvet (2019) derived exact formulas requiring
#' O(N³) computation; this function uses simulation instead.
#'
#' Standard error of estimated π_kl ≈ sqrt(π_kl(1-π_kl)/nsim).
#' With nsim=10000, SE ≈ 0.005 for π_kl ≈ 0.5.
#'
#' @references
#' Chauvet, G. (2019). Properties of Chromy's sampling procedure.
#' \emph{arXiv:1912.10896}.
#'
#' @seealso [up_chromy()] for sampling, [up_brewer_joint()] for an
#'   analytical approximation suitable for high-entropy designs.
#'
#' @examples
#' x <- c(10, 20, 15, 25, 30)
#' joint <- up_chromy_joint(x, n = 3, nsim = 5000)
#'
#' # Diagonal ≈ first-order probabilities
#' pik <- 3 * x / sum(x)
#' diag(joint)  # Should be close to pik
#'
#' @export
up_chromy_joint <- function(x, n, nsim = 10000L) {
  check_mos(x)

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("n must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    stop("nsim must be a positive integer", call. = FALSE)
  }

  .Call(C_up_chromy_joint, as.double(x), as.integer(n), as.integer(nsim))
}
