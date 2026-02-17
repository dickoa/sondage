#' Inclusion Probabilities
#'
#' Compute inclusion probabilities from a size measure, or extract them
#' from a without-replacement design object.
#'
#' @param x A numeric vector of positive size measures, or a
#'   without-replacement design object (class `"wor"`).
#' @param n The desired sample size. Required when `x` is a numeric vector,
#'   ignored when `x` is a design object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of inclusion probabilities. When applied to a
#'   design object, returns the \strong{target} inclusion probabilities
#'   (i.e., the `pik` vector passed to [unequal_prob_wor()]). For most
#'   methods (`cps`, `brewer`, `systematic`, `poisson`), these are the
#'   exact first-order inclusion probabilities of the design. For order
#'   sampling methods (`sps`, `pareto`), the true first-order inclusion
#'   probabilities are approximately but not exactly equal to the target
#'   for finite populations; the discrepancy vanishes as N grows.
#'
#' @seealso [expected_hits()] for the with-replacement analogue,
#'   [unequal_prob_wor()] for sampling with these probabilities.
#'
#' @examples
#' # From size measures
#' size <- c(10, 20, 30, 40)
#' pik <- inclusion_prob(size, n = 2)
#' sum(pik)  # 2
#'
#' # From a design object
#' s <- unequal_prob_wor(pik, method = "cps")
#' inclusion_prob(s)
#'
#' @export
inclusion_prob <- function(x, ...) UseMethod("inclusion_prob")

#' @rdname inclusion_prob
#' @export
inclusion_prob.wor <- function(x, ...) x$pik

#' @rdname inclusion_prob
#' @export
inclusion_prob.wr <- function(x, ...) {
  stop(
    "'inclusion_prob()' applies to without-replacement (WOR) designs only. ",
    "For with-replacement designs, use 'expected_hits()' instead.",
    call. = FALSE
  )
}

#' Expected Hits
#'
#' Compute expected hits from a size measure, or extract them from a
#' with-replacement design object.
#'
#' @param x A numeric vector of positive size measures, or a
#'   with-replacement design object (class `"wr"`).
#' @param n The desired sample size. Required when `x` is a numeric vector,
#'   ignored when `x` is a design object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of expected hits (may exceed 1 for WR designs).
#'
#' @seealso [inclusion_prob()] for the without-replacement analogue,
#'   [unequal_prob_wr()] for sampling with expected hits.
#'
#' @examples
#' # From size measures
#' x <- c(40, 80, 50, 60, 70)
#' hits <- expected_hits(x, n = 3)
#' sum(hits)  # 3
#'
#' # From a design object
#' s <- unequal_prob_wr(hits, method = "chromy")
#' expected_hits(s)
#'
#' @export
expected_hits <- function(x, ...) UseMethod("expected_hits")

#' @rdname expected_hits
#' @export
expected_hits.default <- function(x, n, ...) {
  if (missing(n)) {
    stop("'n' is required when 'x' is not a design object", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value", call. = FALSE)
  }
  if (is.na(n) || n < 0) {
    stop("'n' must be non-negative and not NA", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("there are missing values in 'x'", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' values must be finite (no Inf or NaN)", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("'x' values must be non-negative", call. = FALSE)
  }
  if (sum(x) == 0) {
    stop("sum of 'x' must be positive", call. = FALSE)
  }
  n * x / sum(x)
}

#' @rdname expected_hits
#' @export
expected_hits.wr <- function(x, ...) x$n * x$prob

#' @rdname expected_hits
#' @export
expected_hits.wor <- function(x, ...) {
  stop(
    "'expected_hits()' applies to with-replacement (WR) designs only. ",
    "For without-replacement designs, use 'inclusion_prob()' instead.",
    call. = FALSE
  )
}

#' Joint Inclusion Probabilities
#'
#' Computes the matrix of joint inclusion probabilities
#' \eqn{\pi_{ij} = P(i \in S \text{ and } j \in S)} for a
#' without-replacement sampling design.
#'
#' @param x A without-replacement design object (class `"wor"`).
#' @param ... Additional arguments passed to methods (e.g., `eps`
#'   for boundary tolerance).
#'
#' @details
#' The computation depends on the sampling method stored in the design
#' object. Not all methods yield exact joint probabilities:
#'
#' \describe{
#'   \item{Exact}{`cps` (from the CPS design matrix), `systematic`
#'     (combinatorial enumeration via C code), `poisson`
#'     (\eqn{\pi_{ij} = \pi_i \pi_j}, independent selections),
#'     `srs`, and `bernoulli`.}
#'   \item{Approximation}{`brewer`, `sps`, and `pareto` use the
#'     \strong{high-entropy approximation}:
#'     \eqn{\pi_{ij} \approx \pi_i \pi_j
#'       (1 - (1-\pi_i)(1-\pi_j) / \sum_k \pi_k(1-\pi_k))}.
#'     This is very accurate for these high-entropy designs
#'     (Hajek, 1964; Brewer & Donadio, 2003).}
#' }
#'
#' For \strong{systematic PPS} sampling, some off-diagonal entries may be
#' exactly zero (pairs of units that can never co-occur in the same
#' systematic sample). This has consequences for variance estimation;
#' see [sampling_cov()].
#'
#' @return A symmetric N x N matrix of joint inclusion probabilities.
#'   Diagonal entries are the first-order inclusion probabilities
#'   \eqn{\pi_i}.
#'
#' @seealso [joint_expected_hits()] for the with-replacement analogue,
#'   [sampling_cov()] for the covariance matrix.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' s <- unequal_prob_wor(pik, method = "cps")
#' joint_inclusion_prob(s)
#'
#' @export
joint_inclusion_prob <- function(x, ...) UseMethod("joint_inclusion_prob")

#' @rdname joint_inclusion_prob
#' @param eps Boundary tolerance for CPS and systematic methods
#'   (default 1e-6).
#' @export
joint_inclusion_prob.wor <- function(x, eps = 1e-6, ...) {
  pik <- x$pik
  N <- x$N
  n <- x$n

  switch(
    x$method,
    cps = .Call(C_cps_jip, as.double(pik), as.double(eps)),
    brewer = .Call(C_high_entropy_jip, as.double(pik)),
    systematic = {
      if (inherits(x, "equal_prob")) {
        # Equal prob systematic: same as SRS
        .jip_srs(n, N)
      } else {
        .Call(C_up_systematic_jip, as.double(pik), as.double(eps))
      }
    },
    sps = .Call(C_high_entropy_jip, as.double(pik)),
    pareto = .Call(C_high_entropy_jip, as.double(pik)),
    poisson = {
      pikl <- outer(pik, pik)
      diag(pikl) <- pik
      pikl
    },
    srs = .jip_srs(n, N),
    bernoulli = {
      p <- pik[1]
      pikl <- matrix(p * p, N, N)
      diag(pikl) <- p
      pikl
    },
    stop(
      sprintf("joint_inclusion_prob not implemented for method '%s'", x$method),
      call. = FALSE
    )
  )
}

#' @rdname joint_inclusion_prob
#' @export
joint_inclusion_prob.wr <- function(x, ...) {
  stop(
    "'joint_inclusion_prob()' applies to without-replacement (WOR) designs only. ",
    "For with-replacement designs, use 'joint_expected_hits()' instead.",
    call. = FALSE
  )
}

#' @rdname joint_inclusion_prob
#' @export
joint_inclusion_prob.default <- function(x, ...) {
  stop(
    "'joint_inclusion_prob()' requires a WOR design object created by ",
    "'equal_prob_wor()' or 'unequal_prob_wor()'.",
    call. = FALSE
  )
}

#' Joint Expected Hits
#'
#' Computes the matrix of pairwise expectations \eqn{E(n_i n_j)} for a
#' with-replacement sampling design, where \eqn{n_k} is the number of
#' times unit \eqn{k} is selected.
#'
#' @param x A with-replacement design object (class `"wr"`).
#' @param ... Additional arguments passed to methods (e.g., `nsim`
#'   for simulation-based methods).
#'
#' @details
#' The computation depends on the sampling method:
#'
#' \describe{
#'   \item{Exact analytic}{`multinomial`
#'     (\eqn{E(n_i n_j) = n(n-1) p_i p_j}) and `srs`
#'     (\eqn{E(n_i n_j) = n(n-1)/N^2}).}
#'   \item{Simulation-based}{`chromy`: pairwise expectations are
#'     estimated by Monte Carlo simulation (controlled by the `nsim`
#'     parameter, default 10 000). Increase `nsim` for more precise
#'     estimates at the cost of computation time.}
#' }
#'
#' @return A symmetric N x N matrix. Diagonal entries are
#'   \eqn{E(n_i^2)} and off-diagonal entries are \eqn{E(n_i n_j)}.
#'
#' @seealso [joint_inclusion_prob()] for the without-replacement analogue,
#'   [sampling_cov()] for the covariance matrix.
#'
#' @examples
#' x <- c(40, 80, 50, 60, 70)
#' hits <- expected_hits(x, n = 3)
#' s <- unequal_prob_wr(hits, method = "chromy")
#' joint_expected_hits(s)
#'
#' @export
joint_expected_hits <- function(x, ...) UseMethod("joint_expected_hits")

#' @rdname joint_expected_hits
#'
#' @param nsim Number of simulations for Chromy's pairwise expectations
#'   (default 10000).
#'
#' @export
joint_expected_hits.wr <- function(x, nsim = 10000L, ...) {
  prob <- x$prob
  n <- x$n
  N <- x$N

  switch(
    x$method,
    chromy = {
      if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
        stop("'nsim' must be a positive integer", call. = FALSE)
      }
      nsim <- check_integer(nsim, "nsim")
      .Call(C_chromy_joint_exp, as.double(prob), as.integer(n), nsim)
    },
    multinomial = {
      pikl <- n * (n - 1) * outer(prob, prob)
      diag(pikl) <- n * prob
      pikl
    },
    srs = {
      p <- 1 / N
      pikl <- matrix(n * (n - 1) * p * p, N, N)
      diag(pikl) <- n * p
      pikl
    },
    stop(
      sprintf("joint_expected_hits not implemented for method '%s'", x$method),
      call. = FALSE
    )
  )
}

#' @rdname joint_expected_hits
#' @export
joint_expected_hits.wor <- function(x, ...) {
  stop(
    "'joint_expected_hits()' applies to with-replacement (WR) designs only. ",
    "For without-replacement designs, use 'joint_inclusion_prob()' instead.",
    call. = FALSE
  )
}

#' @rdname joint_expected_hits
#' @export
joint_expected_hits.default <- function(x, ...) {
  stop(
    "'joint_expected_hits()' requires a WR design object created by ",
    "'equal_prob_wr()' or 'unequal_prob_wr()'.",
    call. = FALSE
  )
}

#' Sampling Covariance Matrix
#'
#' Computes the sampling covariance matrix used in variance estimation.
#'
#' For without-replacement designs:
#' \eqn{\Delta_{ij} = \pi_{ij} - \pi_i \pi_j}.
#'
#' For with-replacement designs:
#' \eqn{E(n_i n_j) - E(n_i) E(n_j)}.
#'
#' When `scaled = TRUE`, returns the Sen-Yates-Grundy check quantities:
#' \eqn{1 - \pi_i \pi_j / \pi_{ij}} for WOR,
#' \eqn{1 - E(n_i) E(n_j) / E(n_i n_j)} for WR.
#'
#' @param x A sampling design object (class `"sondage_sample"`).
#' @param scaled If `FALSE` (default), returns the raw covariance matrix.
#'   If `TRUE`, returns the scaled check quantities used in the
#'   Sen-Yates-Grundy variance estimator.
#' @param ... Additional arguments passed to [joint_inclusion_prob()] or
#'   [joint_expected_hits()].
#'
#' @details
#' \strong{Exact vs. approximate joint probabilities.}
#' The accuracy of the covariance matrix depends on the accuracy of
#' the underlying joint probabilities. For `cps`, `systematic`, `poisson`,
#' `srs`, and `bernoulli`, the joint probabilities are exact and so is
#' the covariance matrix. For `brewer`, `sps`, and `pareto`, the joint
#' probabilities are based on the high-entropy approximation, so the
#' covariance matrix is also approximate. For `chromy`, the pairwise
#' expectations are simulation-based (controlled by `nsim`).
#' See [joint_inclusion_prob()] and [joint_expected_hits()] for details.
#'
#' \strong{Zero joint inclusion probabilities.}
#' Some designs (notably systematic PPS) can produce \eqn{\pi_{ij} = 0}
#' for pairs of units that never co-occur in the same sample. When
#' `scaled = TRUE`, the quantity \eqn{1 - \pi_i \pi_j / \pi_{ij}} is
#' undefined for such pairs. These entries are set to `NA` and a
#' warning is issued. The raw covariance (`scaled = FALSE`) is
#' unaffected, since \eqn{\Delta_{ij} = 0 - \pi_i \pi_j} is finite.
#'
#' \strong{Implications for variance estimation.}
#' The Sen-Yates-Grundy variance estimator requires all pairwise
#' \eqn{\pi_{ij} > 0} in the observed sample. It is not applicable for
#' designs with zero joint probabilities (a well-known limitation;
#' see Tille, 2006, Ch. 5). Consider alternative variance estimators
#' for such designs, e.g. successive-differences or Hartley-Rao
#' approximations.
#'
#' @return A symmetric N x N matrix. For WOR designs with `scaled = FALSE`,
#'   off-diagonal entries are typically negative for well-behaved designs.
#'   With `scaled = TRUE`, off-diagonal entries are typically non-positive
#'   (entries where \eqn{\pi_{ij} = 0} are set to `NA`).
#'
#' @references
#' Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
#'   estimator. \emph{Memorial JSM}.
#'
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' s <- unequal_prob_wor(pik, method = "cps")
#'
#' # Raw covariance
#' sampling_cov(s)
#'
#' # SYG check quantities
#' sampling_cov(s, scaled = TRUE)
#'
#' @export
sampling_cov <- function(x, ...) UseMethod("sampling_cov")

#' @rdname sampling_cov
#' @export
sampling_cov.wor <- function(x, scaled = FALSE, ...) {
  pikl <- joint_inclusion_prob(x, ...)
  pik <- x$pik
  if (scaled) {
    pip <- outer(pik, pik)
    zero <- pikl == 0 & pip > 0
    if (any(zero[lower.tri(zero)])) {
      warning(
        "Some joint inclusion probabilities are zero while marginal ",
        "probabilities are positive. The Sen-Yates-Grundy variance ",
        "estimator is not applicable for this design (e.g. systematic ",
        "sampling). Affected entries are set to NA. Consider using an ",
        "approximation-based variance estimator instead.",
        call. = FALSE
      )
    }
    m <- 1 - pip / pikl
    m[zero] <- NA_real_
    diag(m) <- 1 - pik
    m
  } else {
    pikl - outer(pik, pik)
  }
}

#' @rdname sampling_cov
#' @export
sampling_cov.wr <- function(x, scaled = FALSE, ...) {
  pikl <- joint_expected_hits(x, ...)
  ehits <- expected_hits(x)
  if (scaled) {
    m <- 1 - outer(ehits, ehits) / pikl
    diag(m) <- 1 - ehits / diag(pikl)
    m
  } else {
    pikl - outer(ehits, ehits)
  }
}

#' @rdname sampling_cov
#' @export
sampling_cov.default <- function(x, ...) {
  stop(
    "'sampling_cov()' requires a design object created by one of the ",
    "sampling functions ('equal_prob_wor()', 'equal_prob_wr()', ",
    "'unequal_prob_wor()', or 'unequal_prob_wr()').",
    call. = FALSE
  )
}

#' @noRd
.jip_srs <- function(n, N) {
  if (n < 2 || N < 2) {
    pikl <- matrix(0, N, N)
    diag(pikl) <- n / N
    return(pikl)
  }
  pi_ij <- n * (n - 1) / (N * (N - 1))
  pikl <- matrix(pi_ij, N, N)
  diag(pikl) <- n / N
  pikl
}
