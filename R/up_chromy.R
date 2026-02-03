#' Chromy's Sequential PPS Sampling
#'
#' Selects a sample using Chromy's (1979) sequential method with probability
#' proportional to size. This is the default METHOD=PPS_SEQ in SAS SURVEYSELECT.
#'
#' @param x A numeric vector of positive size measures (e.g., population,
#'   revenue, area). Must be non-negative with positive sum.
#' @param n The sample size (number of selections).
#'
#' @return An integer vector of selected indices (1 to `length(x)`).
#'   May contain repeated values when expected hits exceed 1 (minimum
#'   replacement).
#'
#' @details
#' Chromy's method is a strictly sequential algorithm that processes units
#' in order and makes an immediate selection decision for each. It achieves
#' spatial balancing similar to systematic sampling.
#'
#' Properties:
#' \itemize{
#'   \item Fixed sample size n
#'   \item Exact expected hits: \eqn{E[hits_k] = n \times x_k / \sum x}
#'   \item Spatial balance (sample spread throughout frame)
#'   \item O(N) time complexity (single pass)
#'   \item All joint inclusion probabilities > 0
#' }
#'
#' The method uses "minimum replacement": if expected hits for unit k is 2.3,
#' the unit appears exactly 2 or 3 times (never 0, 1, or 4+). This differs
#' from multinomial sampling where any count is possible.
#'
#' When all expected hits are < 1 (i.e., \code{n * max(x) / sum(x) <= 1}),
#' the method behaves as without replacement sampling.
#'
#' @references
#' Chromy, J.R. (1979). Sequential sample selection methods.
#' \emph{Proceedings of the Survey Research Methods Section, ASA}, 401-406.
#'
#' Chauvet, G. (2019). Properties of Chromy's sampling procedure.
#' \emph{arXiv:1912.10896}.
#'
#' @seealso [up_multinomial()] for PPS with replacement (any hit count),
#'   [up_systematic()] for systematic PPS,
#'   [up_brewer()] for Brewer's method (WOR only),
#'   [inclusion_prob()] for computing inclusion probabilities
#'
#' @examples
#' # Size measures
#' x <- c(40, 80, 50, 60, 70)
#'
#' # WOR case: n small relative to x
#' set.seed(42)
#' up_chromy(x, n = 3)  # No repeats
#'
#' # Minimum replacement: larger n causes repeats
#' up_chromy(x, n = 10)  # Some units appear multiple times
#'
#' # Verify expected hits
#' n <- 10
#' expected <- n * x / sum(x)
#' print(expected)  # 1.33, 2.67, 1.67, 2.00, 2.33
#'
#' # Simulate to check
#' hits <- table(factor(up_chromy(x, n), levels = 1:5))
#' print(hits)  # Should be close to floor or ceiling of expected
#'
#' @export
up_chromy <- function(x, n) {
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
  if (is.na(n) || n < 1) {
    stop("n must be at least 1", call. = FALSE)
  }
  .Call(C_up_chromy, as.double(x), as.integer(n))
}
