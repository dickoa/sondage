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
  if (is.na(n) || n < 0) {
    stop("'n' must be non-negative and not NA", call. = FALSE)
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
    warning("there are ", sum(neg), " negative value(s) shifted to zero")
  }
  .Call(C_inclusion_prob, x, as.double(n))
}
