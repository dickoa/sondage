#' Equal Probability Sampling Without Replacement
#'
#' Draws a sample with equal inclusion probabilities, without replacement.
#'
#' @param N Population size (positive integer).
#' @param n Expected sample size. For `"srs"` and `"systematic"`, must be
#'   a non-negative integer not exceeding N. For `"bernoulli"`, this is the
#'   expected sample size and `p = n/N` is used as the selection probability.
#' @param method The sampling method:
#'   \describe{
#'     \item{`"srs"`}{Simple Random Sampling. Each possible sample of size n
#'       has equal probability. Fixed sample size.}
#'     \item{`"systematic"`}{Systematic sampling with interval `k = N/n`.
#'       A random start is drawn from `(0, k]`. Fixed sample size.
#'       Implicit stratification based on unit ordering.}
#'     \item{`"bernoulli"`}{Bernoulli sampling. Each unit selected
#'       independently with probability `p = n/N`. Random sample size.
#'       Note: the realized sample size varies across draws.}
#'   }
#' @param nrep Number of replicate samples (default 1). When `nrep > 1`,
#'   `$sample` holds a matrix (fixed-size) or list (random-size) of all
#'   replicates. The design object and all generics remain usable.
#' @param prn Optional vector of permanent random numbers (length N,
#'   values in the open interval (0, 1)) for sample coordination.
#'   Only supported by `"bernoulli"` method. Cannot be used with
#'   `nrep > 1` (identical PRN would produce identical replicates).
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `c("equal_prob", "wor", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector. When `nrep > 1`,
#'   `$sample` is a matrix (n x nrep) for fixed-size methods, or a list
#'   of integer vectors of varying lengths for `"bernoulli"`.
#'
#' @seealso [equal_prob_wr()] for with-replacement designs,
#'   [unequal_prob_wor()] for unequal probability designs.
#'
#' @examples
#' set.seed(1)
#' s <- equal_prob_wor(10, 3)
#' s$sample
#'
#' # Systematic sampling
#' s <- equal_prob_wor(12, 3, method = "systematic")
#' s$sample
#'
#' # Bernoulli sampling (random size, expected n = 30)
#' s <- equal_prob_wor(100, 30, method = "bernoulli")
#' length(s$sample)
#'
#' @export
equal_prob_wor <- function(
  N,
  n,
  method = c("srs", "systematic", "bernoulli"),
  nrep = 1L,
  prn = NULL,
  ...
) {
  method <- match.arg(method)
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn) && method != "bernoulli") {
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
      srs = .srs_wor_sample(N, n, ...),
      systematic = .systematic_ep_sample(N, n, ...),
      bernoulli = .bernoulli_sample(N, n, prn = prn, ...)
    )
  } else {
    .batch_ep_wor(N, n, method, nrep, prn, ...)
  }
}

#' Equal Probability Sampling With Replacement
#'
#' Draws a simple random sample with replacement.
#'
#' @param N Population size (positive integer).
#' @param n Sample size (non-negative integer).
#' @param method The sampling method. Currently only `"srs"`.
#' @param nrep Number of replicate samples (default 1).
#' @param prn Optional vector of permanent random numbers for sample
#'   coordination. Not currently used by any equal-probability WR method;
#'   a warning is issued if provided.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `c("equal_prob", "wr", "sondage_sample")`.
#'   When `nrep = 1`, `$sample` is an integer vector and `$hits` is an
#'   integer vector. When `nrep > 1`, `$sample` is a matrix (n x nrep)
#'   and `$hits` is a matrix (N x nrep).
#'
#' @seealso [equal_prob_wor()] for without-replacement designs,
#'   [unequal_prob_wr()] for unequal probability designs.
#'
#' @examples
#' set.seed(1)
#' s <- equal_prob_wr(10, 3)
#' s$sample
#' s$hits
#'
#' @export
equal_prob_wr <- function(N, n, method = "srs", nrep = 1L, prn = NULL, ...) {
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
    .srs_wr_sample(N, n, ...)
  } else {
    .batch_ep_wr(N, n, method, nrep, prn, ...)
  }
}

#' @noRd
.srs_wor_sample <- function(N, n, ...) {
  .check_ep_args(N, n, replace = FALSE)
  N <- check_integer(N, "N")
  n <- check_integer(n, "n")

  if (n == 0L) {
    idx <- integer(0)
  } else {
    idx <- sample.int(N, n)
  }

  .new_wor_sample(
    sample = idx,
    pik = rep(n / N, N),
    n = n,
    N = N,
    method = "srs",
    fixed_size = TRUE,
    prob_class = "equal_prob"
  )
}

#' @noRd
.systematic_ep_sample <- function(N, n, ...) {
  .check_ep_args(N, n, replace = FALSE)
  N <- check_integer(N, "N")
  n <- check_integer(n, "n")

  if (n == 0L) {
    idx <- integer(0)
  } else {
    k <- N / n
    u <- runif(1, 0, k)
    idx <- as.integer(ceiling(u + k * (0:(n - 1))))
  }

  .new_wor_sample(
    sample = idx,
    pik = rep(n / N, N),
    n = n,
    N = N,
    method = "systematic",
    fixed_size = TRUE,
    prob_class = "equal_prob"
  )
}

#' @noRd
.bernoulli_sample <- function(N, n, prn = NULL, ...) {
  # Reuse shared N/n validation; replace = TRUE since n = N is valid (p = 1)
  .check_ep_args(N, n, replace = TRUE)
  N <- check_integer(N, "N")
  # n is the expected sample size (not necessarily integer), but must not exceed N
  if (n > N) {
    stop("'n' cannot exceed 'N'", call. = FALSE)
  }

  p <- n / N
  if (is.null(prn)) {
    sel <- as.logical(rbinom(N, 1, p))
  } else {
    check_prn(prn, N)
    sel <- prn < p
  }
  idx <- which(sel)

  .new_wor_sample(
    sample = idx,
    pik = rep(p, N),
    n = n,
    N = N,
    method = "bernoulli",
    fixed_size = FALSE,
    prob_class = "equal_prob"
  )
}

#' @noRd
.srs_wr_sample <- function(N, n, ...) {
  .check_ep_args(N, n, replace = TRUE)
  N <- check_integer(N, "N")
  n <- check_integer(n, "n")

  prob <- rep(1 / N, N)

  if (n == 0L) {
    idx <- integer(0)
  } else {
    idx <- sample.int(N, n, replace = TRUE)
  }
  realized_hits <- tabulate(idx, nbins = N)

  .new_wr_sample(
    sample = idx,
    prob = prob,
    hits = realized_hits,
    n = n,
    N = N,
    method = "srs",
    fixed_size = TRUE,
    prob_class = "equal_prob"
  )
}

#' @noRd
.check_ep_args <- function(N, n, replace = FALSE) {
  if (!is.numeric(N) || length(N) != 1 || is.na(N) || N < 1) {
    stop("'N' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0) {
    stop("'n' must be a non-negative integer", call. = FALSE)
  }
  if (!replace && n > N) {
    stop(
      "'n' cannot exceed 'N' when sampling without replacement",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
.batch_ep_wor <- function(N, n, method, nrep, prn, ...) {
  N_int <- check_integer(N, "N")
  fixed_size <- method != "bernoulli"

  if (method == "bernoulli") {
    p <- n / N_int
    sample_data <- lapply(seq_len(nrep), function(i) {
      .bernoulli_sample(N, n, prn = prn, ...)$sample
    })
    pik <- rep(p, N_int)
  } else {
    n_int <- check_integer(n, "n")
    mat <- matrix(0L, n_int, nrep)
    draw_fn <- switch(
      method,
      srs = .srs_wor_sample,
      systematic = .systematic_ep_sample
    )
    for (i in seq_len(nrep)) {
      mat[, i] <- draw_fn(N, n, ...)$sample
    }
    sample_data <- mat
    pik <- rep(n_int / N_int, N_int)
  }

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = n,
    N = N_int,
    method = method,
    fixed_size = fixed_size,
    prob_class = "equal_prob"
  )
}

#' @noRd
.batch_ep_wr <- function(N, n, method, nrep, prn, ...) {
  N_int <- check_integer(N, "N")
  n_int <- check_integer(n, "n")
  prob <- rep(1 / N_int, N_int)

  sample_mat <- matrix(0L, n_int, nrep)
  hits_mat <- matrix(0L, N_int, nrep)
  for (i in seq_len(nrep)) {
    d <- .srs_wr_sample(N, n, ...)
    sample_mat[, i] <- d$sample
    hits_mat[, i] <- d$hits
  }

  .new_wr_sample(
    sample = sample_mat,
    prob = prob,
    hits = hits_mat,
    n = n_int,
    N = N_int,
    method = "srs",
    fixed_size = TRUE,
    prob_class = "equal_prob"
  )
}
