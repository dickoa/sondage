#' Unequal Probability Sampling Without Replacement
#'
#' Draws a sample with unequal inclusion probabilities, without replacement.
#'
#' @param pik A numeric vector of inclusion probabilities.
#'   For fixed-size methods, `sum(pik)` must be an integer to
#'   floating-point accuracy: an exact fixed-size design cannot have a
#'   non-integer sum, so looser sums are rejected rather than silently
#'   rounded. Units with `pik` of exactly 0 are never selected and
#'   units with exactly 1 are always selected; values in between --
#'   however close to the boundary -- are sampled as given.
#' @param method The sampling method:
#'   \describe{
#'     \item{`"cps"`}{Conditional Poisson Sampling (maximum entropy;
#'       Chen et al., 1994). Fixed size, exact joint probabilities
#'       with all \eqn{\pi_{ij} > 0}. Calibration and the conditional
#'       probability table are O(Nn) (probability-domain
#'       Poisson-binomial recurrence); each draw is O(N) thereafter.
#'       Equal `pik` are drawn directly as SRS.}
#'     \item{`"brewer"`}{Brewer's (1975) draw-by-draw method. Fixed
#'       size, approximate joint probabilities (high-entropy
#'       approximation; see [joint_inclusion_prob()]). O(Nn).}
#'     \item{`"systematic"`}{Systematic PPS. Fixed size, exact joint
#'       probabilities but \strong{some may be zero} (pairs that never
#'       co-occur), making the SYG estimator inapplicable; see
#'       [sampling_cov()]. O(N).}
#'     \item{`"poisson"`}{Poisson sampling. \strong{Random} sample
#'       size (expected \eqn{n = \sum \pi_k}). Units selected
#'       independently, so \eqn{\pi_{ij} = \pi_i \pi_j}. Supports
#'       PRN. O(N).}
#'     \item{`"sps"`}{Sequential Poisson Sampling (Ohlsson, 1998).
#'       Order sampling with key \eqn{\xi_k = u_k / \pi_k}; the
#'       \eqn{n} smallest are selected. Fixed size, high-entropy.
#'       Supports PRN. Approximate joint probabilities. The true
#'       first-order inclusion probabilities are approximately equal
#'       to the supplied `pik`; see [inclusion_prob()].
#'       Expected O(N). Tied keys (possible with duplicated `prn`
#'       values) are broken toward the smallest population index.}
#'     \item{`"pareto"`}{Pareto sampling (Rosen, 1997). Order
#'       sampling with odds-ratio key
#'       \eqn{\xi_k = [u_k/(1-u_k)] / [\pi_k/(1-\pi_k)]}. Same
#'       properties as `"sps"`. Expected O(N).}
#'   }
#' @param nrep Number of replicate samples (default 1). When `nrep > 1`,
#'   `$sample` holds a matrix (fixed-size) or list (random-size) of all
#'   replicates. The design object and all generics remain usable.
#' @param prn Optional vector of permanent random numbers (length N,
#'   values in the open interval (0, 1)) for sample coordination.
#'   Supported by methods `"sps"`, `"pareto"`, and `"poisson"`.
#'   When `NULL`, random numbers are generated internally.
#'   Cannot be used with `nrep > 1` (identical PRN would produce
#'   identical replicates). Use a loop with different PRN vectors
#'   for coordinated repeated sampling.
#' @param ... Additional arguments passed to methods registered via
#'   [register_method()]. Built-in methods take no additional
#'   arguments; the former `eps` boundary-trimming argument was
#'   removed because it silently changed the design.
#'
#' @return An object of class `c("unequal_prob", "wor", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector of selected unit indices.
#'   When `nrep > 1`, `$sample` is a matrix (n x nrep) for fixed-size methods,
#'   or a list of integer vectors of varying lengths for random-size methods (`"poisson"`).
#'   `$n` is an integer for fixed-size methods (realized size) and a
#'   double for `"poisson"` (expected size, `sum(pik)`); see
#'   [sondage_sample].
#'
#' @references
#' Chen, X. H., Dempster, A. P., & Liu, J. S. (1994). Weighted finite
#'   population sampling to maximize entropy. \emph{Biometrika}, 81(3),
#'   457-469.
#'
#' Brewer, K.R.W. (1975). A simple procedure for sampling pi-ps wor.
#'   \emph{Australian Journal of Statistics}, 17(3), 166-172.
#'
#' Ohlsson, E. (1998). Sequential Poisson sampling. \emph{Journal of
#'   Official Statistics}, 14(2), 149-162.
#'
#' Rosen, B. (1997). On sampling with probability proportional to size.
#'   \emph{Journal of Statistical Planning and Inference}, 62(2), 159-191.
#'
#' \enc{Tillé}{Tille}, Y. (2006). \emph{Sampling Algorithms}. Springer.
#'
#' @details
#' **Near-certainty inclusion probabilities (CPS).** The CPS
#' fixed-point calibration converges geometrically for well-spread
#' `pik`, but
#' asymptotes at a non-zero defect when some `pik` are within a few
#' decimal digits of 0 or 1 (e.g. 0.9999). When this happens the
#' function emits a "CPS calibration did not reach tolerance" warning
#' reporting the achieved `max_diff`. The realized first-order
#' inclusion probabilities differ from the target by up to `max_diff` --
#' typically 1e-5 or smaller for inputs in the 0.999-range, well within
#' Monte Carlo error for most estimators. If the warning is unwanted,
#' clip `pik` away from 0/1 before calling.
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
#' \donttest{
#' # Batch mode for simulations
#' sim <- unequal_prob_wor(pik, method = "cps", nrep = 1000)
#' dim(sim$sample)  # 2 x 1000
#' }
#'
#' @export
unequal_prob_wor <- function(
  pik,
  method = c("cps", "brewer", "systematic", "poisson", "sps", "pareto"),
  nrep = 1L,
  prn = NULL,
  ...
) {
  if (is.character(method) && length(method) == 1L && is_registered_method(method)) {
    return(.dispatch_registered_wor(pik, method, nrep, prn, ...))
  }
  method <- match.arg(method)
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn) && !.method_supports_prn(method, "wor")) {
    warning(
      sprintf("prn is not used by method '%s' and will be ignored", method),
      call. = FALSE
    )
  }

  if (!is.null(prn) && nrep > 1L) {
    stop(
      "prn and nrep > 1 cannot be used together. ",
      "Permanent random numbers produce identical samples across replicates. ",
      "Use a loop with different prn vectors for coordinated repeated sampling.",
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
      pareto = .pareto_sample(pik, prn = prn, ...),
      .stop_unknown_method(method) # nocov
    )
  } else {
    .batch_wor(pik, method, nrep, ...)
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
#'       expectations \eqn{E(n_i n_j)} are estimated by simulation.
#'       See [joint_expected_hits()]. Complexity: O(N + n).}
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
#'   estimator. \emph{Proceedings of the Survey Research Methods Section,
#'   American Statistical Association}.
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
  if (is.character(method) && length(method) == 1L && is_registered_method(method)) {
    return(.dispatch_registered_wr(hits, method, nrep, prn, ...))
  }
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
      multinomial = .multinomial_sample(hits, ...),
      .stop_unknown_method(method) # nocov
    )
  } else {
    .batch_wr(hits, method, nrep, ...)
  }
}

#' CPS with equal interior probabilities is exactly SRS on those units:
#' the maximum-entropy fixed-size design with equal marginals. Draw it
#' directly instead of building the O(N * n) conditional table.
#'
#' Returns an index vector (nrep = 1), an n x nrep matrix (nrep > 1), or
#' NULL when the fast path does not apply.
#'
#' @noRd
.cps_equal_fast_path <- function(pik, nrep = 1L) {
  valid <- which(pik > 0 & pik < 1)
  if (length(valid) < 2L) {
    return(NULL)
  }
  pv <- pik[valid]
  if (max(pv) - min(pv) > 16 * .Machine$double.eps * max(1, pv[[1L]])) {
    return(NULL)
  }
  certain <- which(pik >= 1)
  n_total <- as.integer(round(sum(pik)))
  n_draw <- n_total - length(certain)
  draw_one <- function(i) {
    sort(c(certain, valid[sample.int(length(valid), n_draw)]))
  }
  if (nrep == 1L) {
    draw_one(1L)
  } else {
    matrix(vapply(seq_len(nrep), draw_one, integer(n_total)), nrow = n_total)
  }
}

#' @noRd
.cps_sample <- function(pik, eps = NULL, ...) {
  if (!is.null(eps)) .stop_eps_removed()
  check_pik(pik, fixed_size = TRUE)
  idx <- .cps_equal_fast_path(pik)
  if (is.null(idx)) {
    idx <- .Call(C_cps_single, as.double(pik))
  }
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
.brewer_sample <- function(pik, eps = NULL, ...) {
  if (!is.null(eps)) .stop_eps_removed()
  check_pik(pik, fixed_size = TRUE)
  idx <- .Call(C_up_brewer, as.double(pik))
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
.systematic_pps_sample <- function(pik, eps = NULL, ...) {
  if (!is.null(eps)) .stop_eps_removed()
  check_pik(pik, fixed_size = TRUE)
  idx <- .Call(C_up_systematic, as.double(pik))
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
  if (!is.null(prn)) check_prn(prn, N)
  idx <- .Call(C_up_poisson, as.double(pik), if (is.null(prn)) NULL else as.double(prn))
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
.sps_sample <- function(pik, prn = NULL, eps = NULL, ...) {
  if (!is.null(eps)) .stop_eps_removed()
  check_pik(pik, fixed_size = TRUE)
  N <- length(pik)
  if (!is.null(prn)) check_prn(prn, N)
  idx <- .Call(C_up_sps, as.double(pik), if (is.null(prn)) NULL else as.double(prn))
  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = as.integer(round(sum(pik))),
    N = N,
    method = "sps",
    fixed_size = TRUE,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.pareto_sample <- function(pik, prn = NULL, eps = NULL, ...) {
  if (!is.null(eps)) .stop_eps_removed()
  check_pik(pik, fixed_size = TRUE)
  N <- length(pik)
  if (!is.null(prn)) check_prn(prn, N)
  idx <- .Call(C_up_pareto, as.double(pik), if (is.null(prn)) NULL else as.double(prn))
  .new_wor_sample(
    sample = idx,
    pik = pik,
    n = as.integer(round(sum(pik))),
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
.batch_wor <- function(pik, method, nrep, ...) {
  # prn is rejected upstream in unequal_prob_wor when nrep > 1, so batch
  # callers never carry a PRN vector, so no prn forwarding here.
  if (!is.null(list(...)[["eps"]])) .stop_eps_removed()
  fixed_size <- .method_is_fixed_size(method, "wor")
  check_pik(pik, fixed_size = fixed_size)

  N <- length(pik)
  n <- sum(pik)
  n_out <- if (fixed_size) as.integer(round(n)) else n

  if (method %in% .batch_optimised_methods) {
    sample_data <- .cps_equal_fast_path(pik, nrep)
    if (is.null(sample_data)) {
      design <- .Call(C_cps_design, as.double(pik))
      sample_data <- .Call(C_cps_draw_batch, design, as.integer(nrep))
    }
  } else if (!fixed_size) {
    sample_data <- lapply(seq_len(nrep), function(i) {
      .poisson_pps_sample(pik, ...)$sample
    })
  } else {
    # Validation already ran above; call the C sampler directly per
    # replicate instead of building a full design object each time.
    n_int <- as.integer(round(n))
    mat <- matrix(0L, n_int, nrep)
    pik_d <- as.double(pik)
    draw_fn <- switch(
      method,
      brewer = function() .Call(C_up_brewer, pik_d),
      systematic = function() .Call(C_up_systematic, pik_d),
      sps = function() .Call(C_up_sps, pik_d, NULL),
      pareto = function() .Call(C_up_pareto, pik_d, NULL),
      .stop_unknown_method(method) # nocov
    )
    for (i in seq_len(nrep)) {
      mat[, i] <- draw_fn()
    }
    sample_data <- mat
  }

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = n_out,
    N = N,
    method = method,
    fixed_size = fixed_size,
    prob_class = "unequal_prob"
  )
}

#' @noRd
.batch_wr <- function(hits, method, nrep, ...) {
  # prn is rejected upstream in unequal_prob_wr (WR methods do not support
  # PRN); no prn forwarding here.
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  prob <- hits / sum(hits)

  sample_mat <- matrix(0L, n, nrep)
  hits_mat <- matrix(0L, N, nrep)
  draw_fn <- switch(
    method,
    chromy = .chromy_sample,
    multinomial = .multinomial_sample,
    .stop_unknown_method(method) # nocov
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
