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
inclusion_prob <- function(x, ...) {
  UseMethod("inclusion_prob")
}

#' @rdname inclusion_prob
#' @export
inclusion_prob.wor <- function(x, ...) {
  x$pik
}

#' @rdname inclusion_prob
#' @export
inclusion_prob.wr <- function(x, ...) {
  stop(
    "'inclusion_prob()' applies to without-replacement (WOR) designs only. ",
    "For with-replacement designs, use 'expected_hits()' instead.",
    call. = FALSE
  )
}

#' @rdname inclusion_prob
#'
#' @details
#' When `x` is a numeric vector of size measures and `n` is provided,
#' computes \strong{exact} inclusion probabilities via an iterative
#' algorithm: initial probabilities \eqn{\pi_k = n x_k / \sum x_k} are
#' computed, then any unit with \eqn{\pi_k \ge 1} is set to 1 (certainty
#' selection) and the remaining probabilities are recomputed with a reduced
#' \eqn{n}. This process repeats until all probabilities are in \eqn{[0, 1]}.
#' The result always sums to exactly \eqn{n}.
#'
#' This differs from [expected_hits()], which uses simple proportional
#' allocation \eqn{n x_k / \sum x_k} without capping -- values can exceed 1.
#'
#' Negative values in `x` are treated as zero (with a warning).
#'
#' @examples
#' # With certainty selections (large units)
#' size <- c(1, 1, 1, 100)
#' pik <- inclusion_prob(size, n = 2)
#' pik  # Unit 4 gets probability 1
#'
#' @export
inclusion_prob.default <- function(x, n, ...) {
  if (missing(n)) {
    stop("'n' is required when 'x' is not a design object", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value", call. = FALSE)
  }
  if (is.na(n)) {
    stop("'n' must not be NA", call. = FALSE)
  }
  if (!is.finite(n)) {
    stop("'n' must be finite", call. = FALSE)
  }
  if (n < 0) {
    stop("'n' must be non-negative", call. = FALSE)
  }
  n <- check_integer(n, "n")

  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector", call. = FALSE)
  }

  if (n > length(x)) {
    stop("'n' cannot exceed length of 'x'", call. = FALSE)
  }

  if (anyNA(x)) {
    stop("there are missing values in 'x'", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' values must be finite (no Inf or NaN)", call. = FALSE)
  }

  storage.mode(x) <- "double"
  neg <- x < 0
  if (any(neg)) {
    warning(
      "there are ", sum(neg), " negative value(s) shifted to zero",
      call. = FALSE
    )
  }
  .Call(C_inclusion_prob, x, as.double(n))
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
expected_hits <- function(x, ...) {
  UseMethod("expected_hits")
}

#' @rdname expected_hits
#' @export
expected_hits.default <- function(x, n, ...) {
  if (missing(n)) {
    stop("'n' is required when 'x' is not a design object", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value", call. = FALSE)
  }
  if (is.na(n)) {
    stop("'n' must not be NA", call. = FALSE)
  }
  if (!is.finite(n)) {
    stop("'n' must be finite", call. = FALSE)
  }
  if (n < 0) {
    stop("'n' must be non-negative", call. = FALSE)
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
expected_hits.wr <- function(x, ...) {
  x$n * x$prob
}

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
#' @param sampled_only If `TRUE`, return only the n x n submatrix for
#'   the sampled units (requires `nrep = 1`). Useful when N is large
#'   but n is manageable. Default `FALSE`.
#' @param ... Additional arguments passed to methods (e.g., `eps`
#'   for boundary tolerance).
#'
#' @details
#' The computation depends on the sampling method stored in the design
#' object. Not all methods yield exact joint probabilities:
#'
#' \describe{
#'   \item{Exact (up to numerical precision)}{`cps` (from the CPS design
#'     matrix via Aires' formula and elementary symmetric polynomials),
#'     `systematic` (combinatorial enumeration via C code), `poisson`
#'     (\eqn{\pi_{ij} = \pi_i \pi_j}, independent selections),
#'     `srs`, and `bernoulli`.}
#'   \item{Approximation}{`brewer`, `sps`, `pareto`, and `cube` use the
#'     \strong{high-entropy approximation} (Brewer & Donadio, 2003,
#'     eq. 18):
#'     \eqn{\pi_{ij} \approx \pi_i \pi_j (c_i + c_j) / 2}
#'     where \eqn{c_k = (n-1) / (n - (2n-1)\pi_k/(n-1) +
#'     \sum_l \pi_l^2/(n-1))}.
#'     This approximation guarantees symmetry,
#'     \eqn{0 \leq \pi_{ij} \leq \min(\pi_i, \pi_j)}, and correct
#'     diagonal (\eqn{\pi_{ii} = \pi_i}), but does \strong{not}
#'     guarantee the fixed-size marginal identity
#'     \eqn{\sum_{j \neq i} \pi_{ij} = (n-1)\pi_i}. The marginal
#'     defect is typically small for well-spread inclusion
#'     probabilities but can become non-trivial for highly skewed
#'     \eqn{\pi_k} vectors, especially when \eqn{n} is small
#'     (e.g. \eqn{n = 2}). A warning is issued when the defect
#'     exceeds 5\% of \eqn{n}. Use `method = "cps"` when exact (up to
#'     numerical precision) second-order inclusion probabilities are
#'     required.}
#' }
#'
#' For \strong{systematic PPS} sampling, some off-diagonal entries may be
#' exactly zero (pairs of units that can never co-occur in the same
#' systematic sample). This has consequences for variance estimation;
#' see [sampling_cov()].
#'
#' When `sampled_only = TRUE`, only the n x n submatrix for sampled
#' units is returned. For methods using the high-entropy approximation
#' (`brewer`, `sps`, `pareto`, `cube`), the \eqn{c_k} coefficients are
#' computed from the full population \eqn{\pi_k} vector (O(N)), but only
#' the sampled pairs are assembled into a matrix (O(n^2)), avoiding the
#' O(N^2) allocation entirely. Similarly, `poisson`, `bernoulli`, and
#' `srs` compute the n x n matrix directly. These methods can therefore
#' handle arbitrarily large N (e.g. N = 50 000 with n = 200).
#'
#' For `cps` and `systematic`, the underlying C code requires the full
#' N x N matrix, so the N > 10 000 guard still applies even with
#' `sampled_only = TRUE`. For large populations with these methods,
#' consider switching to a high-entropy method (e.g. `brewer`).
#'
#' The marginal defect diagnostic is skipped when `sampled_only = TRUE`
#' because the row-sum identity only holds for the full matrix.
#'
#' @return A symmetric N x N matrix (or n x n if `sampled_only = TRUE`)
#'   of joint inclusion probabilities. Diagonal entries are the
#'   first-order inclusion probabilities \eqn{\pi_i}.
#'
#' @seealso [joint_expected_hits()] for the with-replacement analogue,
#'   [sampling_cov()] for the covariance matrix.
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5)
#' s <- unequal_prob_wor(pik, method = "cps")
#' joint_inclusion_prob(s)
#'
#' # Only the n x n submatrix for sampled units
#' joint_inclusion_prob(s, sampled_only = TRUE)
#'
#' @export
joint_inclusion_prob <- function(x, ...) {
  UseMethod("joint_inclusion_prob")
}

#' @rdname joint_inclusion_prob
#' @param eps Boundary tolerance for CPS and systematic methods
#'   (default 1e-6).
#' @export
joint_inclusion_prob.wor <- function(x, sampled_only = FALSE, eps = 1e-6, ...) {
  pik <- x$pik
  N <- x$N
  n <- x$n
  sample_idx <- x$sample
  eps <- check_eps(eps)

  if (sampled_only) {
    if (is.matrix(sample_idx) || is.list(sample_idx)) {
      stop(
        "sampled_only = TRUE is not supported for batch designs (nrep > 1). ",
        "Use sampled_only = FALSE and subset manually:\n",
        "  pikl[s$sample[, i], s$sample[, i]] for fixed-size designs, or\n",
        "  pikl[s$sample[[i]], s$sample[[i]]] for random-size designs.",
        call. = FALSE
      )
    }
    n_sampled <- length(sample_idx)
    if (n_sampled > 10000L) {
      stop(
        sprintf(
          "n = %d is too large for dense joint probability computation (n x n = %.1f GB).",
          n_sampled,
          n_sampled^2 * 8 / 1e9
        ),
        call. = FALSE
      )
    }
    # CPS and systematic still require N x N internally
    if (x$method %in% c("cps", "systematic") && N > 10000L) {
      stop(
        sprintf(
          "N = %d is too large for method '%s' even with sampled_only = TRUE, ",
          N,
          x$method
        ),
        "because this method requires an N x N intermediate matrix. ",
        "Consider using method = \"brewer\" or another high-entropy method, ",
        "which can compute the n x n submatrix directly.",
        call. = FALSE
      )
    }
  } else {
    if (N > 10000L) {
      stop(
        sprintf(
          "N = %d is too large for dense joint probability computation (N x N = %.1f GB). ",
          N,
          N^2 * 8 / 1e9
        ),
        "Consider using sampled_only = TRUE to compute the n x n submatrix, or ",
        "approximation-based variance estimators for large populations.",
        call. = FALSE
      )
    }
  }

  if (sampled_only) {
    pikl <- switch(
      x$method,
      cps = .Call(C_cps_jip, as.double(pik), as.double(eps))[
        sample_idx,
        sample_idx,
        drop = FALSE
      ],
      systematic = .Call(C_up_systematic_jip, as.double(pik), as.double(eps))[
        sample_idx,
        sample_idx,
        drop = FALSE
      ],
      brewer = ,
      sps = ,
      pareto = ,
      cube = .he_jip_sampled(pik, sample_idx, eps),
      poisson = {
        pik_s <- pik[sample_idx]
        J <- outer(pik_s, pik_s)
        diag(J) <- pik_s
        J
      },
      srs = {
        ns <- length(sample_idx)
        .jip_srs(n, N, ns)
      },
      bernoulli = {
        p <- pik[1]
        ns <- length(sample_idx)
        J <- matrix(p * p, ns, ns)
        diag(J) <- p
        J
      },
      stop(
        sprintf(
          "joint_inclusion_prob not implemented for method '%s'",
          x$method
        ),
        call. = FALSE
      )
    )
  } else {
    pikl <- switch(
      x$method,
      cps = .Call(C_cps_jip, as.double(pik), as.double(eps)),
      brewer = .Call(C_high_entropy_jip, as.double(pik), as.double(eps)),
      systematic = .Call(C_up_systematic_jip, as.double(pik), as.double(eps)),
      sps = .Call(C_high_entropy_jip, as.double(pik), as.double(eps)),
      pareto = .Call(C_high_entropy_jip, as.double(pik), as.double(eps)),
      cube = .Call(C_high_entropy_jip, as.double(pik), as.double(eps)),
      poisson = {
        J <- outer(pik, pik)
        diag(J) <- pik
        J
      },
      srs = .jip_srs(n, N),
      bernoulli = {
        p <- pik[1]
        J <- matrix(p * p, N, N)
        diag(J) <- p
        J
      },
      stop(
        sprintf(
          "joint_inclusion_prob not implemented for method '%s'",
          x$method
        ),
        call. = FALSE
      )
    )

    # Marginal defect diagnostic for high-entropy approximation
    if (x$method %in% c("brewer", "sps", "pareto", "cube")) {
      n_round <- round(n)
      if (n_round >= 2L) {
        defect <- max(abs(rowSums(pikl) - n * pik))
        if (defect / n > 0.05) {
          warning(
            sprintf(
              "High-entropy approximation: marginal defect = %.4f (%.1f%% of n). ",
              defect,
              100 * defect / n
            ),
            "The marginal identity sum(pi_ij, j!=i) = (n-1)*pi_i is not well ",
            "satisfied for this pik vector. Consider using method = \"cps\" for ",
            "exact joint inclusion probabilities.",
            call. = FALSE
          )
        }
      }
    }
  }

  pikl
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
#' @param sampled_only If `TRUE`, return only the submatrix for units
#'   selected at least once (requires `nrep = 1`). Useful when N is
#'   large but the number of distinct selected units is manageable.
#'   Default `FALSE`.
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
#' When `sampled_only = TRUE`, the submatrix is indexed by population
#' units that were selected at least once (i.e., units with
#' `hits > 0`). For `multinomial` and `srs`, the submatrix is computed
#' directly at O(n_s^2). For `chromy`, the full N x N matrix is
#' computed internally and then subset.
#'
#' @return A symmetric N x N matrix (or n_s x n_s if
#'   `sampled_only = TRUE`, where n_s is the number of distinct
#'   selected units). Diagonal entries are \eqn{E(n_i^2)} and
#'   off-diagonal entries are \eqn{E(n_i n_j)}.
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
#' # Only the submatrix for selected units
#' joint_expected_hits(s, sampled_only = TRUE)
#'
#' @export
joint_expected_hits <- function(x, ...) UseMethod("joint_expected_hits")

#' @rdname joint_expected_hits
#'
#' @param nsim Number of simulations for Chromy's pairwise expectations
#'   (default 10000).
#'
#' @export
joint_expected_hits.wr <- function(
  x,
  sampled_only = FALSE,
  nsim = 10000L,
  ...
) {
  prob <- x$prob
  n <- x$n
  N <- x$N

  if (sampled_only) {
    if (is.matrix(x$sample) || is.list(x$sample)) {
      stop(
        "sampled_only = TRUE is not supported for batch designs (nrep > 1). ",
        "Use sampled_only = FALSE and subset manually:\n",
        "  pikl[s$sample[, i], s$sample[, i]] for fixed-size designs, or\n",
        "  pikl[s$sample[[i]], s$sample[[i]]] for random-size designs.",
        call. = FALSE
      )
    }
    sample_idx <- which(x$hits > 0)
    n_sampled <- length(sample_idx)
    if (n_sampled > 10000L) {
      stop(
        sprintf(
          "n = %d is too large for dense joint probability computation (n x n = %.1f GB).",
          n_sampled,
          n_sampled^2 * 8 / 1e9
        ),
        call. = FALSE
      )
    }
    # Chromy still requires N x N internally
    if (x$method == "chromy" && N > 10000L) {
      stop(
        sprintf(
          "N = %d is too large for method 'chromy' even with sampled_only = TRUE, ",
          N
        ),
        "because this method requires an N x N intermediate matrix. ",
        "Consider using method = \"multinomial\" which can compute the ",
        "n x n submatrix directly.",
        call. = FALSE
      )
    }
  } else {
    if (N > 10000L) {
      stop(
        sprintf(
          "N = %d is too large for dense joint probability computation (N x N = %.1f GB). ",
          N,
          N^2 * 8 / 1e9
        ),
        "Consider using sampled_only = TRUE to compute the n x n submatrix, or ",
        "approximation-based variance estimators for large populations.",
        call. = FALSE
      )
    }
  }

  if (sampled_only) {
    prob_s <- prob[sample_idx]
    switch(
      x$method,
      chromy = {
        if (
          !is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1
        ) {
          stop("'nsim' must be a positive integer", call. = FALSE)
        }
        nsim <- check_integer(nsim, "nsim")
        full <- .Call(C_chromy_joint_exp, as.double(prob), as.integer(n), nsim)
        full[sample_idx, sample_idx, drop = FALSE]
      },
      multinomial = {
        pikl <- n * (n - 1) * outer(prob_s, prob_s)
        diag(pikl) <- n * prob_s * (1 - prob_s) + (n * prob_s)^2
        pikl
      },
      srs = {
        p <- 1 / N
        ns <- n_sampled
        pikl <- matrix(n * (n - 1) * p * p, ns, ns)
        diag(pikl) <- n * p * (1 - p) + (n * p)^2
        pikl
      },
      stop(
        sprintf(
          "joint_expected_hits not implemented for method '%s'",
          x$method
        ),
        call. = FALSE
      )
    )
  } else {
    switch(
      x$method,
      chromy = {
        if (
          !is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1
        ) {
          stop("'nsim' must be a positive integer", call. = FALSE)
        }
        nsim <- check_integer(nsim, "nsim")
        .Call(C_chromy_joint_exp, as.double(prob), as.integer(n), nsim)
      },
      multinomial = {
        pikl <- n * (n - 1) * outer(prob, prob)
        diag(pikl) <- n * prob * (1 - prob) + (n * prob)^2
        pikl
      },
      srs = {
        p <- 1 / N
        pikl <- matrix(n * (n - 1) * p * p, N, N)
        diag(pikl) <- n * p * (1 - p) + (n * p)^2
        pikl
      },
      stop(
        sprintf(
          "joint_expected_hits not implemented for method '%s'",
          x$method
        ),
        call. = FALSE
      )
    )
  }
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
#' When `weighted = TRUE`, returns the Sen-Yates-Grundy check quantities:
#' \eqn{1 - \pi_i \pi_j / \pi_{ij}} for WOR,
#' \eqn{1 - E(n_i) E(n_j) / E(n_i n_j)} for WR.
#'
#' @param x A sampling design object (class `"sondage_sample"`).
#' @param weighted If `FALSE` (default), returns the raw covariance matrix.
#'   If `TRUE`, returns the weighted check quantities used in the
#'   Sen-Yates-Grundy variance estimator.
#' @param sampled_only If `TRUE`, compute only the submatrix for sampled
#'   units. Passed to [joint_inclusion_prob()] or [joint_expected_hits()].
#'   Default `FALSE`.
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
#' `weighted = TRUE`, the quantity \eqn{1 - \pi_i \pi_j / \pi_{ij}} is
#' undefined for such pairs. These entries are set to `NA` and a
#' warning is issued. The raw covariance (`weighted = FALSE`) is
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
#' @return A symmetric N x N matrix (or n x n if `sampled_only = TRUE`).
#'   For WOR designs with `weighted = FALSE`, off-diagonal entries are
#'   typically negative for well-behaved designs. With `weighted = TRUE`,
#'   off-diagonal entries are typically non-positive (entries where
#'   \eqn{\pi_{ij} = 0} are set to `NA`).
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
#' sampling_cov(s, weighted = TRUE)
#'
#' # Covariance for sampled units only
#' sampling_cov(s, sampled_only = TRUE)
#'
#' @export
sampling_cov <- function(x, ...) {
  UseMethod("sampling_cov")
}

#' @rdname sampling_cov
#' @export
sampling_cov.wor <- function(x, weighted = FALSE, sampled_only = FALSE, ...) {
  pikl <- joint_inclusion_prob(x, sampled_only = sampled_only, ...)
  if (sampled_only) {
    pik <- x$pik[x$sample]
  } else {
    pik <- x$pik
  }
  if (weighted) {
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
    # Handle 0/0 case (both pikl and pip are zero)
    both_zero <- pikl == 0 & pip == 0
    m <- 1 - pip / pikl
    m[zero] <- NA_real_
    m[both_zero] <- NA_real_
    diag(m) <- 1 - pik
    m
  } else {
    pikl - outer(pik, pik)
  }
}

#' @rdname sampling_cov
#' @export
sampling_cov.wr <- function(x, weighted = FALSE, sampled_only = FALSE, ...) {
  pikl <- joint_expected_hits(x, sampled_only = sampled_only, ...)
  if (sampled_only) {
    sample_idx <- which(x$hits > 0)
    ehits <- expected_hits(x)[sample_idx]
  } else {
    ehits <- expected_hits(x)
  }
  if (weighted) {
    eep <- outer(ehits, ehits)
    zero <- pikl == 0 & eep > 0
    if (any(zero[lower.tri(zero)])) {
      warning(
        "Some pairwise expected hits are zero while marginal expected ",
        "hits are positive. The weighted covariance is undefined for ",
        "these entries. Affected entries are set to NA.",
        call. = FALSE
      )
    }
    both_zero <- pikl == 0 & eep == 0
    m <- 1 - eep / pikl
    m[zero] <- NA_real_
    m[both_zero] <- NA_real_
    d <- diag(pikl)
    diag(m) <- ifelse(d == 0, NA_real_, 1 - ehits^2 / d)
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
.jip_srs <- function(n, N, ns = N) {
  if (n < 2 || N < 2) {
    pikl <- matrix(0, ns, ns)
    diag(pikl) <- n / N
    return(pikl)
  }
  pi_ij <- n * (n - 1) / (N * (N - 1))
  pikl <- matrix(pi_ij, ns, ns)
  diag(pikl) <- n / N
  pikl
}

#' @noRd
.he_jip_sampled <- function(pik, sample_idx, eps = 1e-6) {
  pik_s <- pik[sample_idx]
  n_s <- length(pik_s)

  pikl <- matrix(0, n_s, n_s)
  diag(pikl) <- pik_s

  # Handle certainty units in sample
  cert_s <- which(pik_s >= 1 - eps)
  if (length(cert_s) > 0) {
    for (ci in cert_s) {
      pikl[ci, ] <- pik_s
      pikl[, ci] <- pik_s
      pikl[ci, ci] <- pik_s[ci]
    }
  }

  # Valid (non-certainty, non-zero) units from FULL population
  valid <- pik > eps & pik < 1 - eps
  n_sum <- sum(pik[valid])
  sum_pik2 <- sum(pik[valid]^2)

  if (n_sum <= 1 + eps) {
    return(pikl)
  }

  nm1 <- n_sum - 1
  coef1 <- (2 * n_sum - 1) / nm1
  coef2 <- sum_pik2 / nm1
  ck <- numeric(length(pik))
  ck[valid] <- nm1 / (n_sum - coef1 * pik[valid] + coef2)

  # Valid units in sample
  valid_s <- which(pik_s > eps & pik_s < 1 - eps)
  if (length(valid_s) < 2) {
    return(pikl)
  }

  pik_v <- pik_s[valid_s]
  ck_v <- ck[sample_idx[valid_s]]
  J <- outer(pik_v, pik_v) * (outer(ck_v, ck_v, "+") / 2)
  J <- pmin(J, outer(pik_v, pik_v, pmin))
  J[J < 0] <- 0
  diag(J) <- pik_v
  pikl[valid_s, valid_s] <- J

  pikl
}
