#' Compute Inclusion Probabilities from Measure of Size
#'
#' Converts a measure of size (MOS) into first-order inclusion probabilities
#' that sum to the desired sample size n.
#'
#' @param x A numeric vector of positive size measures (e.g., population,
#'   revenue, area). Negative values are treated as zero (with a warning).
#' @param n The desired sample size (sum of inclusion probabilities).
#'
#' @return A numeric vector of inclusion probabilities between 0 and 1 that
#'   sum to n. Units with very large size measures may have probability 1
#'   (certainty selections).
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Compute initial probabilities: \eqn{\pi_k = n \cdot x_k / \sum x_k}
#'   \item If any \eqn{\pi_k \geq 1}, set to 1 and redistribute remaining
#'     sample size to other units
#'   \item Iterate until all \eqn{\pi_k} are in the valid range
#' }
#'
#' This is the standard "probability proportional to size" (PPS) approach
#' with automatic handling of certainty selections.
#'
#' @seealso [up_maxent()], [up_brewer()], [up_systematic()] for sampling
#'   with these inclusion probabilities
#'
#' @examples
#' # Simple example
#' size <- c(10, 20, 30, 40)
#' pik <- inclusion_prob(size, n = 2)
#' pik
#' sum(pik)  # Should be 2
#'
#' # With certainty selections (large units)
#' size <- c(1, 1, 1, 100)  # Unit 4 is much larger
#' pik <- inclusion_prob(size, n = 2)
#' pik  # Unit 4 gets probability 1
#'
#' # Use with sampling
#' set.seed(42)
#' population_size <- c(500, 1200, 800, 3000, 600)
#' pik <- inclusion_prob(population_size, n = 3)
#' idx <- up_maxent(pik)
#' idx  # Selected indices
#'
#' @export
inclusion_prob <- function(x, n) {
  if (is.na(n) || n < 0) {
    stop("'n' must be non-negative and not NA", call. = FALSE)
  }

  if (!is.numeric(x) || length(n) != 1) {
    stop("'n' must be a single numeric value", call. = FALSE)
  }

  if (n > length(x)) {
    stop("'n' cannot exceed length of 'x'", call. = FALSE)
  }

  storage.mode(x) <- "double"
  neg <- x < 0 & !is.na(x)
  if (any(neg)) {
    warning("there are ", sum(neg), " negative value(s) shifted to zero")
  }
  .Call(C_inclusion_prob, x, as.double(n))
}
