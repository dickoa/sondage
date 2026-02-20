#' Unequal Probability Sampling Without Replacement
#'
#' Draws a sample with unequal inclusion probabilities, without replacement.
#'
#' @param pik A numeric vector of inclusion probabilities.
#'   For fixed-size methods, `sum(pik)` must be close to an integer.
#' @param method The sampling method:
#'   \describe{
#'     \item{`"cps"`}{Conditional Poisson Sampling (maximum entropy design;
#'       Chen, Dempster & Liu, 1994). Fixed sample size, exact first-order
#'       inclusion probabilities. Joint inclusion probabilities are exact
#'       (computed from the CPS design matrix); all \eqn{\pi_{ij} > 0}.
#'       Complexity: O(N^2).}
#'     \item{`"brewer"`}{Brewer's (1963) draw-by-draw method. Fixed sample
#'       size, exact first-order inclusion probabilities. Joint inclusion
#'       probabilities are \emph{approximated} via the high-entropy
#'       approximation (see [joint_inclusion_prob()]). Complexity: O(Nn).}
#'     \item{`"systematic"`}{Systematic PPS sampling. Fixed sample size,
#'       exact first-order inclusion probabilities. Joint inclusion
#'       probabilities are exact but \strong{some may be zero} (pairs of
#'       units that can never co-occur in the same systematic sample).
#'       This makes the Sen-Yates-Grundy variance estimator inapplicable;
#'       see [sampling_cov()]. Complexity: O(N).}
#'     \item{`"poisson"`}{Poisson sampling. \strong{Random} sample size
#'       (expected size \eqn{n = \sum \pi_k}). Each unit is selected
#'       independently with probability \eqn{\pi_k}, so joint inclusion
#'       probabilities are exact: \eqn{\pi_{ij} = \pi_i \pi_j}.
#'       Supports PRN for sample coordination. Complexity: O(N).
#'       Note: the realized sample size varies across draws and may
#'       occasionally be zero, particularly when inclusion probabilities
#'       are small.}
#'     \item{`"sps"`}{Sequential Poisson Sampling (Ohlsson, 1998).
#'       Implemented as order sampling (Rosen, 1997) with ranking key
#'       \eqn{\xi_k = u_k / \pi_k}: the \eqn{n} units with the smallest
#'       \eqn{\xi_k} are selected. This is equivalent to Ohlsson's
#'       sequential threshold adjustment. Fixed sample size, high-entropy
#'       design. Supports PRN for sample coordination. First-order
#'       inclusion probabilities are approximately (not exactly) equal
#'       to the target `pik` for finite populations; see
#'       [inclusion_prob()]. Joint inclusion probabilities are
#'       \emph{approximated} via the high-entropy approximation.
#'       Complexity: O(N log N).}
#'     \item{`"pareto"`}{Pareto sampling (Rosen, 1997). Order sampling
#'       with odds-ratio ranking key
#'       \eqn{\xi_k = [u_k/(1-u_k)] / [\pi_k/(1-\pi_k)]}. Fixed sample
#'       size, high-entropy design. Supports PRN for sample coordination.
#'       First-order inclusion probabilities are approximately (not
#'       exactly) equal to the target `pik` for finite populations; see
#'       [inclusion_prob()]. Joint inclusion probabilities are
#'       \emph{approximated} via the high-entropy approximation.
#'       Complexity: O(N log N).}
#'   }
#' @param nrep Number of replicate samples (default 1). When `nrep > 1`,
#'   `$sample` holds a matrix (fixed-size) or list (random-size) of all
#'   replicates. The design object and all generics remain usable.
#' @param prn Optional vector of permanent random numbers (length N,
#'   values in the open interval (0, 1)) for sample coordination.
#'   Supported by methods `"sps"`, `"pareto"`, and `"poisson"`.
#'   When `NULL`, random numbers are generated internally.
#' @param ... Additional arguments passed to methods (e.g., `eps` for
#'   boundary tolerance).
#'
#' @return An object of class `c("unequal_prob", "wor", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector of selected unit indices.
#'   When `nrep > 1`, `$sample` is a matrix (n x nrep) for fixed-size methods,
#'   or a list of integer vectors of varying lengths for random-size methods (`"poisson"`).
#'
#' @references
#' Chen, S. X., Dempster, A. P., & Liu, J. S. (1994). Weighted finite
#'   population sampling to maximize entropy. \emph{Biometrika}, 81(3),
#'   457-469.
#'
#' Brewer, K.R.W. (1963). A model of systematic sampling with unequal
#'   probabilities. \emph{Australian Journal of Statistics}, 5, 5-13.
#'
#' Ohlsson, E. (1998). Sequential Poisson sampling. \emph{Journal of
#'   Official Statistics}, 14(2), 149-162.
#'
#' Rosen, B. (1997). On sampling with probability proportional to size.
#'   \emph{Journal of Statistical Planning and Inference}, 62(2), 159-191.
#'
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer.
#'
#' @seealso [unequal_prob_wr()] for with-replacement designs,
#'   [equal_prob_wor()] for equal probability designs,
#'   [inclusion_prob()] to compute inclusion probabilities from size measures.
#'
#' @examples
#' pik <- c(0.2, 0.4, 0.6, 0.8)
#'
#' # Conditional Poisson Sampling
#' set.seed(123)
#' s <- unequal_prob_wor(pik, method = "cps")
#' s$sample
#'
#' # Brewer's method
#' s <- unequal_prob_wor(pik, method = "brewer")
#' s$sample
#'
#' # Sequential Poisson Sampling with PRN coordination
#' prn <- runif(4)
#' s <- unequal_prob_wor(pik, method = "sps", prn = prn)
#' s$sample
#'
#' # Pareto sampling
#' s <- unequal_prob_wor(pik, method = "pareto", prn = prn)
#' s$sample
#'
#' # Batch mode for simulations
#' sim <- unequal_prob_wor(pik, method = "cps", nrep = 1000)
#' dim(sim$sample)  # 2 x 1000
#'
#' @export
unequal_prob_wor <- function(
  pik,
  method = c("cps", "brewer", "systematic", "poisson", "sps", "pareto"),
  nrep = 1L,
  prn = NULL,
  ...
) {
  method <- match.arg(method)
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn) && !method %in% c("sps", "pareto", "poisson")) {
    warning(
      sprintf("prn is not used by method '%s' and will be ignored", method),
      call. = FALSE
    )
  }

  if (nrep == 1L) {
    switch(
      method,
      cps = .cps_sample(pik, ...),
      brewer = .brewer_sample(pik, ...),
      systematic = .systematic_pps_sample(pik, ...),
      poisson = .poisson_pps_sample(pik, prn = prn, ...),
      sps = .sps_sample(pik, prn = prn, ...),
      pareto = .pareto_sample(pik, prn = prn, ...)
    )
  } else {
    .batch_wor(pik, method, nrep, prn, ...)
  }
}

#' Unequal Probability Sampling With Replacement
#'
#' Draws a sample with unequal selection probabilities, with replacement
#' or minimum replacement.
#'
#' @param hits A numeric vector of expected hits (expected number of
#'   selections per unit). Typically computed via [expected_hits()].
#'   `sum(hits)` must be close to a positive integer.
#' @param method The sampling method:
#'   \describe{
#'     \item{`"chromy"`}{Chromy's (1979) sequential PPS with minimum
#'       replacement. Default method in SAS SURVEYSELECT. Pairwise
#'       expectations \eqn{E(n_i n_j)} are estimated by simulation;
#'       see [joint_expected_hits()]. Complexity: O(N + n).}
#'     \item{`"multinomial"`}{Multinomial PPS (independent draws).
#'       Units can be selected any number of times. Pairwise expectations
#'       are exact: \eqn{E(n_i n_j) = n(n-1) p_i p_j}.
#'       Complexity: O(n).}
#'   }
#' @param nrep Number of replicate samples (default 1).
#' @param prn Optional vector of permanent random numbers for sample
#'   coordination. Not currently used by any WR method; a warning is
#'   issued if provided.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `c("unequal_prob", "wr", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector and `$hits` is an
#'   integer vector. When `nrep > 1`, `$sample` is a matrix (n x nrep) and
#'   `$hits` is a matrix (N x nrep).
#'
#' @references
#' Chromy, J.R. (1979). Sequential sample selection methods.
#'   \emph{Proceedings of the Survey Research Methods Section, ASA}, 401-406.
#'
#' Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
#'   estimator. \emph{Memorial JSM}.
#'
#' @seealso [unequal_prob_wor()] for without-replacement designs,
#'   [expected_hits()] to compute expected hits from size measures.
#'
#' @examples
#' x <- c(40, 80, 50, 60, 70)
#' hits <- expected_hits(x, n = 3)
#'
#' set.seed(42)
#' s <- unequal_prob_wr(hits, method = "chromy")
#' s$sample
#' s$hits
#'
#' @export
unequal_prob_wr <- function(
  hits,
  method = c("chromy", "multinomial"),
  nrep = 1L,
  prn = NULL,
  ...
) {
  method <- match.arg(method)
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn)) {
    warning(
      sprintf("prn is not used by method '%s' and will be ignored", method),
      call. = FALSE
    )
  }

  if (nrep == 1L) {
    switch(
      method,
      chromy = .chromy_sample(hits, ...),
      multinomial = .multinomial_sample(hits, ...)
    )
  } else {
    .batch_wr(hits, method, nrep, prn, ...)
  }
}

#' @noRd
.cps_sample <- function(pik, eps = 1e-06, ...) {
  check_pik(pik, fixed_size = TRUE)
  idx <- .Call(C_cps_single, as.double(pik), as.double(eps))
  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = as.integer(round(sum(pik))),
    N = length(pik),
    method = "cps",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.brewer_sample <- function(pik, eps = 1e-06, ...) {
  check_pik(pik, fixed_size = TRUE)
  idx <- .Call(C_up_brewer, as.double(pik), as.double(eps))
  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = as.integer(round(sum(pik))),
    N = length(pik),
    method = "brewer",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.systematic_pps_sample <- function(pik, eps = 1e-06, ...) {
  check_pik(pik, fixed_size = TRUE)
  certain <- which(pik >= 1 - eps)
  valid <- which(pik > eps & pik < 1 - eps)
  pik_valid <- pik[valid]
  n_valid <- length(pik_valid)

  if (n_valid == 0) {
    idx <- certain
  } else {
    x <- (c(0, cumsum(pik_valid)) - runif(1)) %% 1
    selected_valid <- valid[x[1:n_valid] > x[2:(n_valid + 1)]]
    idx <- sort(c(certain, selected_valid))
  }

  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = as.integer(round(sum(pik))),
    N = length(pik),
    method = "systematic",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.poisson_pps_sample <- function(pik, prn = NULL, ...) {
  check_pik(pik)
  N <- length(pik)
  if (is.null(prn)) {
    sel <- as.logical(rbinom(N, 1, pik))
  } else {
    check_prn(prn, N)
    sel <- prn < pik
  }
  idx <- which(sel)
  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = sum(pik),
    N = N,
    method = "poisson",
    fixed_size = FALSE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.sps_sample <- function(pik, prn = NULL, eps = 1e-06, ...) {
  check_pik(pik, fixed_size = TRUE)
  N <- length(pik)
  n <- as.integer(round(sum(pik)))
  if (is.null(prn)) prn <- stats::runif(N)
  check_prn(prn, N)

  ta <- which(pik >= 1 - eps)
  ts <- which(pik > eps & pik < 1 - eps)
  n_ts <- n - length(ta)

  if (n_ts > 0L && length(ts) > 0L) {
    xi <- prn[ts] / pik[ts]
    ord <- order(xi)[seq_len(n_ts)]
    selected <- sort(c(ta, ts[ord]))
  } else {
    selected <- sort(ta)
  }

  .new_wor_sample(
    sample = selected,
    pik = pik,
    n = n,
    N = N,
    method = "sps",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.pareto_sample <- function(pik, prn = NULL, eps = 1e-06, ...) {
  check_pik(pik, fixed_size = TRUE)
  N <- length(pik)
  n <- as.integer(round(sum(pik)))
  if (is.null(prn)) prn <- stats::runif(N)
  check_prn(prn, N)

  ta <- which(pik >= 1 - eps)
  ts <- which(pik > eps & pik < 1 - eps)
  n_ts <- n - length(ta)

  if (n_ts > 0L && length(ts) > 0L) {
    p <- pik[ts]
    u <- prn[ts]
    xi <- (u / (1 - u)) / (p / (1 - p))
    ord <- order(xi)[seq_len(n_ts)]
    selected <- sort(c(ta, ts[ord]))
  } else {
    selected <- sort(ta)
  }

  .new_wor_sample(
    sample = selected,
    pik = pik,
    n = n,
    N = N,
    method = "pareto",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.chromy_sample <- function(hits, ...) {
  check_hits(hits)
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  prob <- hits / sum(hits)

  idx <- .Call(C_up_chromy, as.double(hits), n)
  realized_hits <- tabulate(idx, nbins = N)

  .new_wr_sample(
    sample = idx,
    prob = prob,
    hits = realized_hits,
    n = n,
    N = N,
    method = "chromy",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.multinomial_sample <- function(hits, ...) {
  check_hits(hits)
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  prob <- hits / sum(hits)

  if (n == 0L) {
    idx <- integer(0L)
  } else {
    idx <- sample.int(N, n, replace = TRUE, prob = prob)
  }
  realized_hits <- tabulate(idx, nbins = N)

  .new_wr_sample(
    sample = idx,
    prob = prob,
    hits = realized_hits,
    n = n,
    N = N,
    method = "multinomial",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.batch_wor <- function(pik, method, nrep, prn, ...) {
  N <- length(pik)
  n <- sum(pik)
  fixed_size <- method != "poisson"

  if (method == "cps") {
    eps <- list(...)$eps
    if (is.null(eps)) eps <- 1e-06
    design <- .Call(C_cps_design, as.double(pik), as.double(eps))
    sample_data <- .Call(C_cps_draw_batch, design, as.integer(nrep))
  } else if (method == "poisson") {
    sample_data <- lapply(seq_len(nrep), function(i) {
      .poisson_pps_sample(pik, prn = prn, ...)$sample
    })
  } else {
    n_int <- as.integer(round(n))
    mat <- matrix(0L, n_int, nrep)
    draw_fn <- switch(
      method,
      brewer = .brewer_sample,
      systematic = .systematic_pps_sample,
      sps = .sps_sample,
      pareto = .pareto_sample
    )
    for (i in seq_len(nrep)) {
      mat[, i] <- draw_fn(pik, prn = prn, ...)$sample
    }
    sample_data <- mat
  }

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = n,
    N = N,
    method = method,
    fixed_size = fixed_size,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.batch_wr <- function(hits, method, nrep, prn, ...) {
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  prob <- hits / sum(hits)

  sample_mat <- matrix(0L, n, nrep)
  hits_mat <- matrix(0L, N, nrep)
  draw_fn <- switch(
    method,
    chromy = .chromy_sample,
    multinomial = .multinomial_sample
  )
  for (i in seq_len(nrep)) {
    d <- draw_fn(hits, ...)
    sample_mat[, i] <- d$sample
    hits_mat[, i] <- d$hits
  }

  .new_wr_sample(
    sample = sample_mat,
    prob = prob,
    hits = hits_mat,
    n = n,
    N = N,
    method = method,
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}
