#' sondage: Survey Sampling Algorithms
#'
#' Fast implementations of survey sampling algorithms for drawing samples
#' from finite populations. All functions return design objects that carry
#' sample indices, inclusion probabilities (or expected hits), and design
#' metadata.
#'
#' @section Unequal Probability Sampling:
#' \itemize{
#'   \item [unequal_prob_wor()] - Without replacement: CPS (maximum entropy),
#'     Brewer, systematic PPS, Poisson, SPS (sequential Poisson),
#'     Pareto
#'   \item [unequal_prob_wr()] - With replacement: Chromy (minimum
#'     replacement), multinomial PPS
#' }
#'
#' @section Equal Probability Sampling:
#' \itemize{
#'   \item [equal_prob_wor()] - Without replacement: SRS, systematic,
#'     Bernoulli (random size)
#'   \item [equal_prob_wr()] - With replacement: SRS
#' }
#'
#' @section Design Queries:
#' All sampling functions return objects of class `"sondage_sample"`.
#' Use these generics to query the design:
#'
#' \itemize{
#'   \item [inclusion_prob()] - First-order inclusion probabilities
#'   \item [expected_hits()] - Expected number of selections (WR)
#'   \item [joint_inclusion_prob()] - Joint inclusion probabilities (WOR)
#'   \item [joint_expected_hits()] - Pairwise expectations (WR)
#'   \item [sampling_cov()] - Sampling covariance matrix
#' }
#'
#' @section Utilities:
#' \itemize{
#'   \item [inclusion_prob()] - Compute inclusion probabilities from size
#'     measures
#'   \item [expected_hits()] - Compute expected hits from size measures
#' }
#'
#' @references
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer Series in
#'   Statistics.
#'
#' Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
#'   estimator. \emph{Memorial JSM}.
#'
#' @importFrom stats runif rbinom
#'
#' @useDynLib sondage, .registration = TRUE
"_PACKAGE"
