#' Random Samples and Permutations
#'
#' `sample_int` takes a sample of the specified size from `1:n` using either with or without replacement.
#'
#' @importFrom stats runif
#'
#' @param n A positive number, the number of items to choose from. See "Details."
#' @param size A non-negative integer giving the number of items to
#'   choose. The default is `n` except for unequal-probability sampling
#'   without replacement, where it is `sum(prob)`.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the elements
#'   of the vector being sampled.
#' @param useHash Logical indicating if the hash-version of
#'   the algorithm should be used. Can only be used for `replace = FALSE`,
#'   `prob = NULL`, and `size <= n/2`, and really should be used for large `n`,
#'   as `useHash=FALSE` will use memory proportional to `n`.
#' @param method When `replace=FALSE` and `prob` is not `NULL`, the
#'   method to use for sampling, see Details below.
#'
#' @details
#' It is allowed to ask for `size = 0` samples with `n = 0` or
#' a length-zero `x`, but otherwise `n > 0` or positive
#' `length(x)` is required.
#'
#' Non-integer positive numerical values of `n` will be
#' truncated to the next smallest integer, which has to be no larger than
#' `.Machine$integer.max`.
#'
#' The optional `prob` argument can be used to give a vector of
#' weights for obtaining the elements of the vector being sampled.
#' They need not sum to one, but they should be non-negative and not all zero.
#' If `replace` is true, Walker's alias method (Ripley, 1987) is
#' used when there are more than 200 reasonably probable values: this
#' gives results incompatible with those from R < 2.2.0.
#'
#' If `replace` is false and `prob` is supplied there are three
#' options, controlled by `method`. The number of nonzero weights
#' must be at least `size` in this case and the weights should sum
#' to the desired sample size.
#'
#' All three of the methods have disadvantages. The default,
#' compatible with R < 4.4.0 is sequential sampling also know as the Yat, that
#' is,  the probability of choosing the next item is proportional to the
#' weights amongst the remaining items. Sequential sampling is fast,
#' *but the probability of being sampled is not equal or
#' proportional to `prob`*. Using `"marginal"`
#' draws a sample so that the inclusion probabilities are the values of
#' `prob` and by default `size=sum(prob)`. It uses the Brewer procedures, Algorithm
#' 6.10 from Tille (2006). This is much slower for large `n`
#' than the other methods.  Finally, `"poisson"` does Poisson
#' sampling, where the inclusion probabilities are the values of
#' `prob` but the *sample size is random* with mean `size`
#' rather than being fixed at `size`; it is most useful when
#' `size` is large and the variability is thus relatively small.
#'
#' @return an integer vector of length `size` with
#' elements from `1:n`, or a double vector if
#' \eqn{n \ge 2^{31}}{n >= 2^31}.
#'
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
#' *The New S Language*.
#' Wadsworth & Brooks/Cole.
#'
#' Ripley, B. D. (1987) *Stochastic Simulation*. Wiley.
#'
#' Tille, Y (2006) *Sampling Algorithms* Springer.
#'
#' @seealso
#' `RNGkind(sample.kind = ..)` about random number generation,
#' notably the change of `sample()` results with R version 3.6.0.
#'
#' CRAN package `sampling` for other methods of weighted sampling
#' without replacement.
#'
#' @examples
#' ## R 3.0.0 and later
#' sample_int(1e10, 12, replace = TRUE)
#' sample_int(1e10, 12) # not that there is much chance of duplicates
#'
#' ## R 4.4.0 and later
#' ### Sequential sampling does not give the specified inclusion probability
#' @export
sample_int <- function(n, size = NULL, replace = FALSE, prob = NULL,
                       useHash = (n > 1e+07 && !replace && is.null(prob) && (!is.null(size)) && size <= n/2),
                       method = c("sequential", "marginal", "poisson")) {
  if (replace || is.null(prob)) {
    if (is.null(size)) {
      size <- n
    }
  } else {
    if (is.null(size))
      size <- sum(prob)
  }
  stopifnot(length(n) == 1L)
  if (useHash) {
    stopifnot(is.null(prob), !replace)
    .Internal(sample2(n, size))
  }
  else if (!is.null(prob) && !replace) {
    if (length(prob) != n)
      stop("incorrect number of probabilities")
    method <- match.arg(method)
    switch(method,
           sequential = sample.int(n, size = size, replace = replace, prob = prob, useHash = FALSE),
           marginal = sample_pps(n, size, prob),
           poisson = sample(seq.int(1, n)[runif(n) <= prob * size/sum(prob)]))
  } else {
    ## Will be the .Internal from sample.int
    sample.int(n, size = size, replace = replace, prob = prob, useHash = FALSE)
  }
}
