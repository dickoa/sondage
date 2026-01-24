#' Maximum Entropy Sampling (Conditional Poisson Sampling)
#'
#' Draws samples using the maximum entropy design, also known as Conditional
#' Poisson Sampling (CPS). This is the unique design that maximizes entropy
#' subject to fixed inclusion probabilities.
#'
#' @param pik A numeric vector of inclusion probabilities. The sum should be
#'   an integer representing the desired sample size.
#' @param nrep Number of sample replicates to draw. Default is 1.
#' @param eps A small threshold value for boundary cases. Default is 1e-06.
#'
#' @return If `nrep = 1`, an integer vector of selected indices.
#'   If `nrep > 1`, an integer matrix with n rows and `nrep` columns,
#'   where each column contains the indices for one replicate.
#'
#' @details
#' Maximum entropy sampling has several desirable properties:
#' \itemize{
#'   \item Fixed sample size: exactly `round(sum(pik))` units selected
#'   \item Exact inclusion probabilities: \eqn{E(I_k) = \pi_k}
#'   \item All joint inclusion probabilities are positive: \eqn{\pi_{kl} > 0}
#'   \item Maximum entropy among all designs with fixed \eqn{\pi_k}
#' }
#'
#' For repeated sampling (simulations), use the `nrep` parameter instead
#' of a loop for much better performance. The design is computed once
#' and reused for all replicates.
#'
#' @references
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer Series in Statistics.
#'
#' Chen, S. X., Dempster, A. P., & Liu, J. S. (1994). Weighted finite population
#'   sampling to maximize entropy. \emph{Biometrika}, 81(3), 457-469.
#'
#' @seealso [up_brewer()] for Brewer's method (also has positive joint probs),
#'   [up_systematic()] for systematic PPS (fastest, but some joint probs = 0),
#'   [inclusion_prob()] for computing inclusion probabilities from size measures
#'
#' @examples
#' pik <- c(0.2, 0.4, 0.6, 0.8)  # sum = 2
#'
#' # Single sample
#' set.seed(42)
#' idx <- up_maxent(pik)
#' idx
#'
#' # Select from data frame
#' df <- data.frame(id = 1:4, x = c(10, 20, 30, 40))
#' df[idx, ]
#'
#' # Multiple replicates for simulation
#' set.seed(42)
#' samples <- up_maxent(pik, nrep = 1000)
#' dim(samples)  # 2 x 1000 (n rows, nrep columns)
#'
#' # Verify inclusion probabilities
#' indicators <- apply(samples, 2, function(s) 1:4 %in% s)
#' rowMeans(indicators)  # Should be close to pik
#'
#' @export
up_maxent <- function(pik, nrep = 1L, eps = 1e-06) {
    # Input validation
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
    if (!is.numeric(nrep) || length(nrep) != 1 || nrep < 1) {
      stop("nrep must be a positive integer",
           call. = FALSE)
    }

    nrep <- as.integer(nrep)

    # Create design object (precomputes lambda parameters)
    design <- .Call(C_maxent_design_create, as.double(pik), as.double(eps))

    if (nrep == 1L) {
        # Single sample - return integer vector
        .Call(C_maxent_sample, design, as.integer(10000L))
    } else {
        # Multiple samples - return integer matrix
        .Call(C_maxent_sample_batch, design, nrep, as.integer(10000L))
    }
}
