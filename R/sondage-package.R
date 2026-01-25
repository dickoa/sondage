#' sondage: Survey Sampling Algorithms
#'
#' Fast implementations of survey sampling algorithms for drawing samples
#' from finite populations. All functions return indices for easy subsetting.
#'
#' @section Equal Probability Sampling:
#' \itemize{
#'   \item [srs()] - Simple random sampling (with/without replacement)
#'   \item [systematic()] - Systematic sampling
#'   \item [bernoulli()] - Bernoulli sampling (random size)
#' }
#'
#' @section Unequal Probability Sampling:
#' \itemize{
#'   \item [up_maxent()] - Maximum entropy / Conditional Poisson Sampling
#'   \item [up_brewer()] - Brewer's method
#'   \item [up_systematic()] - Systematic PPS
#'   \item [up_poisson()] - Poisson sampling (random size)
#'   \item [up_multinomial()] - PPS with replacement
#' }
#'
#' @section Utilities:
#' \itemize{
#'   \item [inclusion_prob()] - Compute inclusion probabilities from size measure
#' }
#'
#' @references
#' Tillé, Y. (2006). \emph{Sampling Algorithms}. Springer Series in Statistics.
#'
#' @useDynLib sondage, .registration = TRUE
"_PACKAGE"
