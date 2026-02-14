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
#' @return A numeric vector of inclusion probabilities.
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
        stop("nsim must be a positive integer", call. = FALSE)
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
#' @return A symmetric N x N matrix. For WOR designs with `scaled = FALSE`,
#'   off-diagonal entries are typically negative for well-behaved designs.
#'   With `scaled = TRUE`, off-diagonal entries are typically non-positive.
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
    m <- 1 - outer(pik, pik) / pikl
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
