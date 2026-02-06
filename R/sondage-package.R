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
#' Functions prefixed with `up_` take either inclusion probabilities (`pik`,
#' values in \[0,1\] summing to n) or raw size measures (`x`, non-negative
#' values):
#'
#' **pik interface** (inclusion probabilities):
#' \itemize{
#'   \item [up_maxent()] - Maximum entropy / Conditional Poisson Sampling
#'   \item [up_brewer()] - Brewer's method
#'   \item [up_systematic()] - Systematic PPS
#'   \item [up_poisson()] - Poisson sampling (random size)
#' }
#'
#' **x interface** (raw size measures):
#' \itemize{
#'   \item [up_multinomial()] - PPS with replacement
#'   \item [up_chromy()] - PPS with minimum replacement
#' }
#'
#' Use [inclusion_prob()] to convert size measures to inclusion probabilities.
#'
#' @section Utilities:
#' \itemize{
#'   \item [inclusion_prob()] - Compute inclusion probabilities from size measure
#' }
#'
#' @section Joint Inclusion Probabilities:
#' \itemize{
#'   \item [up_maxent_jip()] - Exact CPS joint probabilities
#'   \item [up_brewer_jip()] - Brewer & Donadio approximation
#'   \item [up_systematic_jip()] - Exact systematic joint probabilities
#'   \item [up_poisson_jip()] - Independent selections
#' }
#'
#' @references
#' Tillé, Y. (2006). \emph{Sampling Algorithms}. Springer Series in Statistics.
#'
#' @useDynLib sondage, .registration = TRUE
"_PACKAGE"
