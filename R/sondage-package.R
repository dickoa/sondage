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
#' @section Balanced Sampling:
#' \itemize{
#'   \item [balanced_wor()] - Cube method (Deville & Tille, 2004) for
#'     balanced sampling with unequal probabilities, with optional
#'     stratification (Chauvet & Tille, 2006)
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
#' @section Joint Probability Approximations:
#' Standalone approximation functions for joint inclusion
#' probabilities, useful as `joint_fn` arguments to
#' [register_method()]:
#'
#' \itemize{
#'   \item [he_jip()] - High-entropy approximation (Brewer &
#'     Donadio, 2003). Recommended default for most designs.
#'   \item [hajek_jip()] - Hajek (1964) approximation based on
#'     conditional Poisson (rejective) sampling theory. Simpler
#'     formula, slightly less accurate.
#' }
#'
#' For without-replacement designs, the stored `pik` vector is the
#' design-defining target. For most methods this equals the true
#' first-order inclusion probabilities. For order-sampling methods
#' (`"sps"`, `"pareto"`), the true probabilities are only approximately
#' equal to the target.
#'
#' @section Size-to-probability conversion:
#' \itemize{
#'   \item [inclusion_prob()] - Compute inclusion probabilities from size
#'     measures (with capping for certainty selections)
#'   \item [expected_hits()] - Compute expected hits from size measures
#'     (simple proportional allocation, no capping)
#' }
#'
#' @references
#' Tille, Y. (2006). \emph{Sampling Algorithms}. Springer Series in
#'   Statistics.
#'
#' Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
#'   estimator. \emph{Proceedings of the Survey Research Methods Section,
#'   American Statistical Association}.
#'
#' @importFrom stats runif rbinom
#'
#' @useDynLib sondage, .registration = TRUE
"_PACKAGE"
