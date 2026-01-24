#' Brewer's Method for Unequal Probability Sampling
#'
#' Selects a sample using Brewer's method with unequal inclusion probabilities.
#' Implements Algorithm 6.10 from Tille's "Sampling Algorithms".
#'
#' @param pik A numeric vector of inclusion probabilities. The sum should be
#'   an integer representing the desired sample size n.
#' @param eps A small threshold value. Units with `pik <= eps` are excluded
#'   and units with `pik >= 1-eps` are always included (certainty selections).
#'   Default is 1e-06.
#'
#' @return An integer vector of selected indices (1 to `length(pik)`).
#'
#' @details
#' Brewer's method is a draw-by-draw procedure that selects n units with
#' prescribed inclusion probabilities. At each draw i, unit k is selected
#' with probability proportional to:
#'
#' \deqn{p_k \propto \pi_k \cdot \frac{(n - a) - \pi_k}{(n - a) - \pi_k (n - i + 1)}}
#'
#' where \eqn{a = \sum \pi_\ell} for already selected units.
#'
#' Properties:
#' \itemize{
#'   \item Fixed sample size n = `round(sum(pik))`
#'   \item Exact inclusion probabilities
#'   \item All joint inclusion probabilities \eqn{\pi_{kl} > 0}
#'   \item Order invariant (result doesn't depend on unit ordering)
#' }
#'
#' @references
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer Series in Statistics.
#'
#' Brewer, K.R.W. (1963). A model of systematic sampling with unequal
#' probabilities. \emph{Australian Journal of Statistics}, 5, 5-13.
#'
#' @seealso [up_maxent()] for maximum entropy sampling,
#'   [up_systematic()] for systematic PPS sampling
#'
#' @examples
#' pik <- c(0.2, 0.4, 0.6, 0.8)  # sum = 2
#'
#' set.seed(42)
#' idx <- up_brewer(pik)
#' idx
#'
#' # Select from data frame
#' df <- data.frame(id = 1:4, x = c(10, 20, 30, 40))
#' df[idx, ]
#'
#' # Verify inclusion probabilities
#' samples <- replicate(5000, up_brewer(pik))
#' indicators <- sapply(samples, function(s) 1:4 %in% s)
#' rowMeans(indicators)  # Should be close to pik
#'
#' @export
up_brewer <- function(pik, eps = 1e-06) {
    if (any(is.na(pik))) {
      stop("there are missing values in the pik vector",
           call. = FALSE)
    }
    if (!is.numeric(pik)) {
      stop("pik must be a numeric vector",
           call. = FALSE)
    }
    if (length(pik) == 0) {
      stop("pik vector is empty",
           call. = FALSE)
    }
    .Call(C_up_brewer, as.double(pik), as.double(eps[1]))
}
