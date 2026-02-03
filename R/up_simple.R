#' Poisson Sampling
#'
#' Each unit is selected independently with its own inclusion probability.
#' Sample size is random.
#'
#' @param pik A numeric vector of inclusion probabilities between 0 and 1.
#'
#' @return An integer vector of selected indices.
#'   Length is random with expected value `sum(pik)`.
#'
#' @details
#' Poisson sampling is the simplest unequal probability method.
#' Each unit k is selected independently with probability \eqn{\pi_k}.
#'
#' Properties:
#' \itemize{
#'   \item Random sample size with expectation \eqn{\sum \pi_k}
#'   \item Exact inclusion probabilities
#'   \item Independent selections
#'   \item Joint probabilities: \eqn{\pi_{kl} = \pi_k \times \pi_l}
#' }
#'
#' @seealso [up_maxent()] for fixed sample size maximum entropy sampling,
#'   [bernoulli()] for equal probability Bernoulli sampling
#'
#' @examples
#' pik <- c(0.2, 0.5, 0.8, 0.3)
#'
#' # Sample size varies
#' set.seed(1)
#' replicate(5, length(up_poisson(pik)))
#'
#' # Expected size = sum(pik) = 1.8
#' mean(replicate(1000, length(up_poisson(pik))))
#'
#' @export
up_poisson <- function(pik) {
  check_pik(pik)
  selected <- as.logical(rbinom(length(pik), 1, pik))
  which(selected)
}

#' Multinomial Sampling (PPS with Replacement)
#'
#' Draws n units with replacement, with selection probability proportional
#' to a size measure. A unit can be selected multiple times.
#'
#' @param x A numeric vector of positive size measures (e.g., population,
#'   revenue, area). Must be non-negative.
#' @param n The number of draws (sample size).
#'
#' @return An integer vector of n selected indices (1 to `length(x)`).
#'   May contain repeated values.
#'
#' @details
#' Multinomial sampling (PPS with replacement) makes n independent draws,
#' each with probability proportional to pik. This is equivalent to
#' `sampling::UPmultinomial()` but returns indices instead of counts.
#'
#' Properties:
#' \itemize{
#'   \item Fixed number of draws n = `round(sum(pik))`
#'   \item Units can be selected multiple times
#'   \item Selection probability per draw: \eqn{p_k = \pi_k / \sum \pi_k}
#' }
#'
#' Useful for Hansen-Hurwitz estimation and bootstrap procedures.
#'
#' @seealso [up_brewer()], [up_maxent()] for without-replacement methods,
#'   [srs()] with `replace = TRUE` for equal probability
#'
#' @examples
#' # Size measures
#' x <- c(10, 20, 30, 40)
#'
#' set.seed(1234)
#' idx <- up_multinomial(x, n = 5)
#' idx  # 5 indices, may repeat
#'
#' table(idx)  # Unit 4 likely appears most often
#'
#' # Use for selection from data frame
#' df <- data.frame(id = 1:4, size = x)
#' df[idx, ]  # 5 rows, some repeated
#'
#' @export
up_multinomial <- function(x, n) {
  if (!is.numeric(x)) {
    stop("x must be a numeric vector", call. = FALSE)
  }
  if (length(x) == 0L) {
    stop("x vector is empty", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("there are missing values in x", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("x values must be non-negative", call. = FALSE)
  }
  if (sum(x) == 0) {
    stop("sum of x must be positive", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L) {
    stop("n must be a single numeric value", call. = FALSE)
  }
  if (is.na(n) || n < 0) {
    stop("n must be non-negative and not NA", call. = FALSE)
  }

  n <- as.integer(n)
  if (n == 0L) {
    return(integer(0L))
  }

  sample.int(length(x), n, replace = TRUE, prob = x)
}

#' Systematic Sampling with Unequal Probabilities
#'
#' Fast systematic sampling using a single random start.
#' Sample size is fixed at `round(sum(pik))`.
#'
#' @param pik A numeric vector of inclusion probabilities. Should sum to
#'   an integer n (the sample size).
#' @param eps Threshold for boundary cases. Default is 1e-06.
#'
#' @return An integer vector of selected indices.
#'
#' @details
#' Systematic sampling is one of the fastest methods with O(N) time.
#' A single random number determines the entire sample.
#'
#' Properties:
#' \itemize{
#'   \item Fixed sample size n = `round(sum(pik))`
#'   \item Exact inclusion probabilities (on average)
#'   \item Very fast: O(N) time
#'   \item Some joint probabilities may be 0 (units far apart are always
#'     selected together or never)
#' }
#'
#' The order of units matters. For best results, sort by an auxiliary
#' variable first (creates implicit stratification).
#'
#' @references
#' Madow, W.G. (1949). On the theory of systematic sampling, II.
#' \emph{Annals of Mathematical Statistics}, 20, 333-354.
#'
#' @seealso [up_brewer()], [up_maxent()] for methods where all joint probs > 0,
#'   [systematic()] for equal probability systematic sampling
#'
#' @examples
#' pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)  # sum = 2
#'
#' set.seed(2)
#' idx <- up_systematic(pik)
#' idx
#'
#' set.seed(123)
#' n_sim <- 10000
#' counts <- integer(5)
#' for (i in 1:n_sim) {
#'  idx <- up_systematic(pik)
#'  counts[idx] <- counts[idx] + 1
#' }
#'
#' counts / n_sim
#'
#'
#' @export
up_systematic <- function(pik, eps = 1e-06) {
  check_pik(pik)

  certain <- which(pik >= 1 - eps)
  zero <- which(pik <= eps)
  valid <- which(pik > eps & pik < 1 - eps)

  pik_valid <- pik[valid]
  n_valid <- length(pik_valid)

  if (n_valid == 0) {
    return(certain)
  }

  x <- (c(0, cumsum(pik_valid)) - runif(1)) %% 1
  selected_valid <- valid[x[1:n_valid] > x[2:(n_valid + 1)]]
  sort(c(certain, selected_valid))
}
