#' High-Entropy Approximation for Joint Inclusion Probabilities
#'
#' Computes the joint inclusion probability matrix using the
#' high-entropy approximation of Brewer & Donadio (2003, eq. 18):
#' \deqn{\pi_{ij} \approx \pi_i \pi_j \frac{c_i + c_j}{2}}
#' where \eqn{c_k = \frac{n-1}{n - \frac{2n-1}{n-1}\pi_k +
#' \frac{\sum_\ell \pi_\ell^2}{n-1}}} and \eqn{n = \sum_k \pi_k}.
#'
#' @param pik Numeric vector of inclusion probabilities
#'   (\eqn{0 \le \pi_k \le 1}).
#' @param sample_idx Integer vector of 1-based indices for the
#'   sampled units, or `NULL` (default) for the full population
#'   matrix. When non-`NULL`, returns the submatrix for those units
#'   only, without allocating the full N x N matrix.
#' @param eps Boundary tolerance (default 1e-6). Units with
#'   \eqn{\pi_k \ge 1 - \varepsilon} are treated as certainty
#'   selections; units with \eqn{\pi_k \le \varepsilon} are
#'   treated as zero.
#' @param ... Additional arguments (ignored). Present so that the
#'   function matches the `joint_fn` signature required by
#'   [register_method()].
#'
#' @details
#' The high-entropy (HE) approximation is the recommended default
#' for designs that are close to maximum entropy, which includes
#' most common unequal-probability without-replacement designs:
#' Brewer, Sampford, Tille, SPS, Pareto, cube, and CPS itself.
#'
#' The approximation guarantees symmetry,
#' \eqn{0 \le \pi_{ij} \le \min(\pi_i, \pi_j)}, and correct
#' diagonal (\eqn{\pi_{ii} = \pi_i}), but does not exactly
#' satisfy the marginal identity
#' \eqn{\sum_{j \neq i} \pi_{ij} = (n-1)\pi_i}. The defect is
#' typically small but grows with skewed \eqn{\pi_k} and small
#' \eqn{n}.
#'
#' Internally, this calls the same C implementation used by
#' [joint_inclusion_prob()] for Brewer, SPS, Pareto, and cube
#' methods. Exported so that custom methods registered via
#' [register_method()] can use it directly as their `joint_fn`:
#'
#' ```
#' register_method("my_method", sample_fn = my_fn, joint_fn = he_jip)
#' ```
#'
#' For a lighter alternative based on conditional Poisson theory,
#' see [hajek_jip()].
#'
#' @return A symmetric matrix of joint inclusion probabilities:
#'   N x N when `sample_idx` is `NULL`, or
#'   `length(sample_idx)` x `length(sample_idx)` otherwise.
#'   Diagonal entries are \eqn{\pi_i}.
#'
#' @references
#' Brewer, K.R.W. and Donadio, M.E. (2003). The high entropy
#'   variance of the Horvitz-Thompson estimator.
#'   \emph{Survey Methodology}, 29(2), 189--196.
#'
#' @seealso [hajek_jip()] for the Hajek approximation,
#'   [joint_inclusion_prob()] for design-based dispatch,
#'   [register_method()] for custom method registration.
#'
#' @examples
#' pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
#'
#' # Full N x N matrix
#' pikl <- he_jip(pik)
#' round(pikl, 4)
#'
#' # Submatrix for specific units
#' he_jip(pik, sample_idx = c(1, 3, 5))
#'
#' # Use as joint_fn in register_method()
#' register_method("my_method", sample_fn = function(pik, n, prn, ...) {
#'   sample.int(length(pik), n, prob = pik)
#' }, joint_fn = he_jip)
#' unregister_method("my_method")
#'
#' @export
he_jip <- function(pik, sample_idx = NULL, eps = 1e-6, ...) {
  check_pik(pik)
  eps <- check_eps(eps)

  if (is.null(sample_idx)) {
    .Call(C_high_entropy_jip, as.double(pik), as.double(eps))
  } else {
    .he_jip_sampled(pik, as.integer(sample_idx), eps)
  }
}


#' Hajek Approximation for Joint Inclusion Probabilities
#'
#' Computes the joint inclusion probability matrix using the
#' Hajek (1964) approximation based on conditional Poisson
#' (rejective) sampling theory:
#' \deqn{\pi_{ij} \approx \pi_i \pi_j
#'   \left[1 - \frac{(1-\pi_i)(1-\pi_j)}{D}\right]}
#' where \eqn{D = \sum_k \pi_k (1 - \pi_k)}.
#'
#' @inheritParams he_jip
#'
#' @details
#' The Hajek approximation is simpler and computationally lighter
#' than the high-entropy approximation ([he_jip()]), but generally
#' slightly less accurate. It is derived from the asymptotic
#' theory of rejective (conditional Poisson) sampling, where the
#' design is obtained by conditioning independent Poisson trials
#' on the total sample size.
#'
#' The formula is valid for any high-entropy design, but is most
#' accurate when the design is close to rejective sampling. For
#' maximum-entropy designs (CPS, Sampford), [he_jip()] tends to
#' give tighter results. In practice, both approximations agree
#' closely for moderate to large populations with well-spread
#' inclusion probabilities.
#'
#' Like [he_jip()], this function matches the `joint_fn`
#' signature required by [register_method()]:
#'
#' ```
#' register_method("my_method", sample_fn = my_fn, joint_fn = hajek_jip)
#' ```
#'
#' @section Properties:
#' \itemize{
#'   \item Symmetric: \eqn{\pi_{ij} = \pi_{ji}}
#'   \item Diagonal: \eqn{\pi_{ii} = \pi_i} (set directly)
#'   \item Bounded: \eqn{0 \le \pi_{ij} \le \min(\pi_i, \pi_j)}
#'     (clamped)
#'   \item Marginal defect:
#'     \eqn{|\sum_{j \ne i} \pi_{ij} - (n-1)\pi_i| =
#'     O(1/N)} for well-spread \eqn{\pi_k}
#' }
#'
#' @return A symmetric matrix of joint inclusion probabilities:
#'   N x N when `sample_idx` is `NULL`, or
#'   `length(sample_idx)` x `length(sample_idx)` otherwise.
#'   Diagonal entries are \eqn{\pi_i}.
#'
#' @references
#' Hajek, J. (1964). Asymptotic theory of rejective sampling with
#'   varying probabilities from a finite population.
#'   \emph{Annals of Mathematical Statistics}, 35(4), 1491--1523.
#'
#' @seealso [he_jip()] for the high-entropy approximation,
#'   [joint_inclusion_prob()] for design-based dispatch,
#'   [register_method()] for custom method registration.
#'
#' @examples
#' pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
#'
#' # Full N x N matrix
#' pikl <- hajek_jip(pik)
#' round(pikl, 4)
#'
#' # Compare with high-entropy approximation
#' he <- he_jip(pik)
#' max(abs(pikl - he))
#'
#' @export
hajek_jip <- function(pik, sample_idx = NULL, eps = 1e-6, ...) {
  check_pik(pik)
  eps <- check_eps(eps)

  N <- length(pik)
  cert <- which(pik >= 1 - eps)
  valid <- which(pik > eps & pik < 1 - eps)

  if (!is.null(sample_idx)) {
    sample_idx <- as.integer(sample_idx)
    ns <- length(sample_idx)
    pik_out <- pik[sample_idx]
  } else {
    ns <- N
    pik_out <- pik
  }

  pikl <- matrix(0, ns, ns)
  diag(pikl) <- pik_out

  # Handle certainty units: pi_ij = pik[other]
  if (length(cert) > 0L) {
    if (!is.null(sample_idx)) {
      cert_s <- which(sample_idx %in% cert)
    } else {
      cert_s <- cert
    }
    for (ci in cert_s) {
      pikl[ci, ] <- pik_out
      pikl[, ci] <- pik_out
      pikl[ci, ci] <- pik_out[ci]
    }
  }

  # n <= 1: no pairs possible, off-diagonal stays zero
  n_sum <- sum(pik[valid])
  if (n_sum <= 1 + eps) {
    return(pikl)
  }

  # D = sum of pi_k * (1 - pi_k) over valid (non-certainty, non-zero) units
  D <- sum(pik[valid] * (1 - pik[valid]))

  if (D < 1e-15 || length(valid) < 2L) {
    return(pikl)
  }

  # Valid output units
  if (!is.null(sample_idx)) {
    valid_s <- which(pik_out > eps & pik_out < 1 - eps)
  } else {
    valid_s <- valid
  }

  if (length(valid_s) < 2L) {
    return(pikl)
  }

  pik_v <- pik_out[valid_s]
  q_v <- 1 - pik_v

  # pi_ij = pi_i * pi_j * [1 - (1 - pi_i)(1 - pi_j) / D]
  J <- outer(pik_v, pik_v) * (1 - outer(q_v, q_v) / D)

  # Clamp to [0, min(pi_i, pi_j)]
  J <- pmin(J, outer(pik_v, pik_v, pmin))
  J[J < 0] <- 0
  diag(J) <- pik_v

  pikl[valid_s, valid_s] <- J
  pikl
}
