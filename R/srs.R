#' Simple Random Sampling
#'
#' Draws a simple random sample of n units from a population of N units.
#'
#' @param n Sample size (number of units to select).
#' @param N Population size.
#' @param replace Logical. Sample with replacement? Default is `FALSE`.
#'
#' @return An integer vector of selected indices (1 to N).
#'   Length is `n`. With `replace = TRUE`, indices may repeat.
#'
#' @details
#' Without replacement (`replace = FALSE`):
#' \itemize{
#'   \item Each possible sample of size n has equal probability
#'   \item Inclusion probability: \eqn{\pi_k = n/N} for all units
#'   \item n must not exceed N
#' }
#'
#' With replacement (`replace = TRUE`):
#' \itemize{
#'   \item Each draw is independent with probability 1/N per unit
#'   \item Same unit can be selected multiple times
#'   \item n can exceed N
#' }
#'
#' @seealso [systematic()] for systematic sampling, [bernoulli()] for Bernoulli sampling
#'
#' @examples
#' # Without replacement
#' set.seed(42)
#' idx <- srs(3, 10)
#' idx
#'
#' # Select from data frame
#' df <- data.frame(id = 1:10, x = rnorm(10))
#' df[idx, ]
#'
#' # With replacement (can have repeats)
#' set.seed(42)
#' idx <- srs(5, 10, replace = TRUE)
#' idx
#' df[idx, ]  # Some rows may appear twice
#'
#' @export
srs <- function(n, N, replace = FALSE) {
  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("n must be a non-negative integer",
         call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1 || N < 1) {
    stop("N must be a positive integer",
         call. = FALSE)
  }
  n <- as.integer(n)
  N <- as.integer(N)

  if (!replace && n > N) {
    stop("n cannot exceed N when replace = FALSE",
         call. = FALSE)
  }

  if (n == 0) {
    return(integer(0))
  }
  sample.int(N, n, replace = replace)
}

#' Systematic Sampling
#'
#' Selects every k-th unit starting from a random position.
#'
#' @importFrom stats runif
#'
#' @param n Sample size.
#' @param N Population size.
#'
#' @return An integer vector of n selected indices (1 to N).
#'
#' @details
#' The sampling interval is k = N/n. A random start u is drawn from (0, k],
#' then units at positions ceiling(u), ceiling(u + k), ceiling(u + 2k), ...
#' are selected.
#'
#' Properties:
#' \itemize{
#'   \item Fixed sample size n
#'   \item Equal inclusion probabilities: \eqn{\pi = n/N}
#'   \item Very fast: O(n) time
#'   \item Implicit stratification based on unit ordering
#' }
#'
#' Note: Order of units matters. Sort by an auxiliary variable first
#' to achieve implicit stratification.
#'
#' @seealso [srs()] for simple random sampling,
#'   [up_systematic()] for systematic sampling with unequal probabilities
#'
#' @examples
#' # Basic usage
#' set.seed(42)
#' idx <- systematic(3, 12)
#' idx  # e.g., c(2, 6, 10) - interval of 4
#'
#' # Implicit stratification
#' df <- data.frame(
#'   id = 1:100,
#'   region = rep(1:10, each = 10)
#' )
#' df <- df[order(df$region), ]
#' idx <- systematic(10, 100)
#' table(df$region[idx])
#'
#' @export
systematic <- function(n, N) {
  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("n must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1 || N < 1) {
    stop("N must be a positive integer", call. = FALSE)
  }
  if (n > N) {
    stop("n cannot exceed N", call. = FALSE)
  }

  n <- as.integer(n)
  N <- as.integer(N)

  if (n == 0) {
    return(integer(0))
  }

  k <- N / n # Sampling interval
  u <- runif(1, 0, k) # Random start in (0, k]
  as.integer(ceiling(u + k * (0:(n - 1))))
}


#' Bernoulli Sampling
#'
#' Each unit is independently selected with the same probability p.
#' Sample size is random.
#'
#' @importFrom stats rbinom
#'
#' @param p Selection probability (same for all units), between 0 and 1.
#' @param N Population size.
#'
#' @return An integer vector of selected indices (1 to N).
#'   Length is random with expected value `N*p`.
#'
#' @details
#' Each unit is selected independently with probability p.
#'
#' Properties:
#' \itemize{
#'   \item Random sample size with expectation N*p and variance N*p*(1-p)
#'   \item Equal inclusion probabilities: \eqn{\pi = p}
#'   \item Very fast: O(N) time
#' }
#'
#' @seealso [srs()] for fixed sample size,
#'   [up_poisson()] for Bernoulli sampling with unequal probabilities
#'
#' @examples
#' # Basic usage
#' set.seed(42)
#' idx <- bernoulli(0.3, 100)
#' length(idx)  # Random, expected = 30
#'
#' # Sample size distribution
#' sizes <- replicate(1000, length(bernoulli(0.3, 100)))
#' mean(sizes)
#' var(sizes)
#'
#' @export
bernoulli <- function(p, N) {
  if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
    stop("p must be a probability between 0 and 1",
         call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1 || N < 1) {
    stop("N must be a positive integer",
         call. = FALSE)
  }
  N <- as.integer(N)
  selected <- as.logical(rbinom(N, 1, p))
  which(selected)
}
