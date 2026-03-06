#' Sampling Design Objects
#'
#' All sampling functions in `sondage` return objects inheriting from class
#' `"sondage_sample"`. These objects store the realized sample together with
#' the design-defining quantities needed by the query generics.
#'
#' @details
#' `sondage_sample` objects always contain these fields:
#' \describe{
#'   \item{`$sample`}{Realized sample indices. For `nrep = 1`, an integer
#'   vector. For `nrep > 1`, a matrix for fixed-size designs or a list of
#'   integer vectors for random-size designs.}
#'   \item{`$n`}{Target or expected sample size. Integer for fixed-size
#'   designs; may be non-integer for random-size designs such as Poisson,
#'   Bernoulli, or stratified cube with non-integer stratum totals.}
#'   \item{`$N`}{Population size.}
#'   \item{`$method`}{Method name string used to generate the design.}
#'   \item{`$fixed_size`}{Logical flag indicating whether the realized sample
#'   size is fixed by design.}
#' }
#'
#' Without-replacement (`"wor"`) objects additionally contain:
#' \describe{
#'   \item{`$pik`}{Design-defining inclusion probability vector. For methods
#'   with exact first-order guarantees, this equals the true first-order
#'   inclusion probabilities. For order-sampling methods such as `"sps"` and
#'   `"pareto"`, this is the stored target vector.}
#' }
#'
#' With-replacement (`"wr"`) objects additionally contain:
#' \describe{
#'   \item{`$prob`}{Selection probability vector for each draw.}
#'   \item{`$hits`}{Realized hit counts by unit. For `nrep = 1`, an integer
#'   vector of length `N`. For `nrep > 1`, an `N x nrep` integer matrix.}
#' }
#'
#' The class vector also records the design family:
#' \describe{
#'   \item{`c("equal_prob", "wor", "sondage_sample")`}{Equal-probability
#'   sampling without replacement.}
#'   \item{`c("equal_prob", "wr", "sondage_sample")`}{Equal-probability
#'   sampling with replacement.}
#'   \item{`c("unequal_prob", "wor", "sondage_sample")`}{Unequal-probability
#'   sampling without replacement, including balanced sampling.}
#'   \item{`c("unequal_prob", "wr", "sondage_sample")`}{Unequal-probability
#'   sampling with replacement or minimum replacement.}
#' }
#'
#' @seealso [equal_prob_wor()], [equal_prob_wr()], [unequal_prob_wor()],
#'   [unequal_prob_wr()], [balanced_wor()], [inclusion_prob()],
#'   [expected_hits()], [joint_inclusion_prob()], [joint_expected_hits()],
#'   [sampling_cov()], [print.sondage_sample()]
#'
#' @name sondage_sample
NULL

#' @keywords internal
#' @noRd
.new_wor_sample <- function(sample, pik, n, N, method, fixed_size, prob_class) {
  if (is.numeric(sample) && !is.matrix(sample)) {
    sample <- as.integer(sample)
  }
  structure(
    list(
      sample = sample,
      pik = pik,
      n = n,
      N = as.integer(N),
      method = method,
      fixed_size = fixed_size
    ),
    class = c(prob_class, "wor", "sondage_sample")
  )
}

#' @keywords internal
#' @noRd
.new_wr_sample <- function(
  sample,
  prob,
  hits,
  n,
  N,
  method,
  fixed_size,
  prob_class
) {
  if (is.numeric(sample) && !is.matrix(sample)) {
    sample <- as.integer(sample)
  }
  if (is.numeric(hits) && !is.matrix(hits)) {
    hits <- as.integer(hits)
  }
  structure(
    list(
      sample = sample,
      prob = prob,
      hits = hits,
      n = n,
      N = as.integer(N),
      method = method,
      fixed_size = fixed_size
    ),
    class = c(prob_class, "wr", "sondage_sample")
  )
}
