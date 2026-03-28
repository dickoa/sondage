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
#'   design object, returns the stored `pik` vector. For most methods
#'   this equals the true first-order inclusion probabilities. For
#'   order-sampling methods (`sps`, `pareto`), it is the target used
#'   to define the design; the true probabilities are approximately
#'   equal and converge as N grows.
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
#' When `x` is a numeric vector and `n` is provided, computes inclusion
#' probabilities via iterative capping: units with \eqn{\pi_k \ge 1}
#' are set to 1 (certainty selections) and the remaining probabilities
#' are recomputed with reduced \eqn{n}. The result sums to exactly
#' \eqn{n}. This differs from [expected_hits()], which does simple
#' proportional allocation without capping. Negative values in `x`
#' are treated as zero (with a warning).
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
      "there are ",
      sum(neg),
      " negative value(s) shifted to zero",
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
#' @return A numeric vector of expected hits. Values can exceed 1 for
#'   with-replacement designs.
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
#' The computation depends on the method stored in the design object:
#'
#' \describe{
#'   \item{Exact}{`cps` (elementary symmetric polynomials),
#'     `systematic` (combinatorial enumeration), `poisson`
#'     (\eqn{\pi_{ij} = \pi_i \pi_j}), `srs`, and `bernoulli`.}
#'   \item{Approximate}{`brewer`, `sps`, `pareto`, and `cube` use the
#'     high-entropy approximation (Brewer & Donadio, 2003, eq. 18):
#'     \eqn{\pi_{ij} \approx \pi_i \pi_j (c_i + c_j) / 2}.
#'     This guarantees symmetry,
#'     \eqn{0 \leq \pi_{ij} \leq \min(\pi_i, \pi_j)}, and correct
#'     diagonal, but \strong{not} the marginal identity
#'     \eqn{\sum_{j \neq i} \pi_{ij} = (n-1)\pi_i}. The defect is
#'     typically small but grows with skewed \eqn{\pi_k} and small
#'     \eqn{n}. A warning is issued when it exceeds 5\% of \eqn{n}.
#'     Use `method = "cps"` when exact second-order probabilities are
#'     needed. For `sps` and `pareto`, the approximation uses the
#'     stored target `pik`, not the exact finite-population
#'     probabilities.}
#' }
#'
#' For \strong{systematic PPS}, some off-diagonal entries may be
#' exactly zero (pairs that never co-occur). See [sampling_cov()].
#'
#' When `sampled_only = TRUE`, only the n x n submatrix for sampled
#' units is returned. All methods compute this directly without
#' allocating the full N x N matrix, so large N is feasible
#' (e.g. N = 50 000 with n = 200). The marginal defect diagnostic
#' is skipped because the row-sum identity only holds for the full
#' matrix.
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
      cps = .Call(
        C_cps_jip_sub,
        as.double(pik),
        as.double(eps),
        as.integer(sample_idx)
      ),
      systematic = .Call(
        C_up_systematic_jip_sub,
        as.double(pik),
        as.double(eps),
        as.integer(sample_idx)
      ),
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
      .registered_joint_or_stop(
        x$method, pik, sample_idx, "joint_inclusion_prob", ...
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
      .registered_joint_or_stop(
        x$method, pik, NULL, "joint_inclusion_prob", ...
      )
    )

    # Marginal defect diagnostic for high-entropy approximation
    if (x$method %in% .he_jip_methods) {
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
#' The computation depends on the method:
#'
#' \describe{
#'   \item{Exact}{`multinomial`
#'     (\eqn{E(n_i n_j) = n(n-1) p_i p_j}) and `srs`
#'     (\eqn{E(n_i n_j) = n(n-1)/N^2}).}
#'   \item{Simulation}{`chromy`: estimated by Monte Carlo (controlled
#'     by `nsim`, default 10 000).}
#' }
#'
#' When `sampled_only = TRUE`, only the submatrix for units with
#' `hits > 0` is returned. All methods compute it directly without
#' allocating the full N x N matrix.
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
        .Call(
          C_chromy_joint_exp_sub,
          as.double(prob),
          as.integer(n),
          nsim,
          as.integer(sample_idx)
        )
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
      .registered_joint_or_stop(
        x$method, prob, sample_idx, "joint_expected_hits", ...
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
      .registered_joint_or_stop(
        x$method, prob, NULL, "joint_expected_hits", ...
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
#' Accuracy depends on the underlying joint probabilities. For `cps`,
#' `systematic`, `poisson`, `srs`, and `bernoulli`, joint probabilities
#' are exact and so is the covariance. For `brewer`, `sps`, `pareto`,
#' and `cube`, they use the high-entropy approximation. For `chromy`,
#' they are simulation-based (see `nsim`). See [joint_inclusion_prob()]
#' and [joint_expected_hits()].
#'
#' Some designs (notably systematic PPS) produce \eqn{\pi_{ij} = 0} for
#' pairs that never co-occur. When `weighted = TRUE`, the SYG quantity
#' \eqn{1 - \pi_i \pi_j / \pi_{ij}} is undefined for such pairs and
#' set to `NA` with a warning. The raw covariance (`weighted = FALSE`)
#' is unaffected. The Sen-Yates-Grundy estimator is not applicable
#' for these designs (Tille, 2006, Ch. 5).
#'
#' @return A symmetric N x N matrix (or n x n if `sampled_only = TRUE`).
#'   For WOR designs with `weighted = FALSE`, off-diagonal entries are
#'   typically negative for well-behaved designs. With `weighted = TRUE`,
#'   off-diagonal entries are typically non-positive (entries where
#'   \eqn{\pi_{ij} = 0} are set to `NA`).
#'
#' @references
#' Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
#'   estimator. \emph{Proceedings of the Survey Research Methods Section,
#'   American Statistical Association}.
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
