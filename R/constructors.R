#' Sampling Design Objects
#'
#' All sampling functions in `sondage` return objects inheriting from class
#' `"sondage_sample"`. These objects store the realized sample together with
#' the design-defining quantities needed by the query generics.
#'
#' @details
#' All `sondage_sample` objects contain:
#' \describe{
#'   \item{`$sample`}{Sample indices (integer vector, or matrix/list
#'   when `nrep > 1`).}
#'   \item{`$n`}{Target or expected sample size.}
#'   \item{`$N`}{Population size.}
#'   \item{`$method`}{Sampling method name.}
#'   \item{`$fixed_size`}{Whether the sample size is fixed by design.}
#' }
#'
#' Without-replacement (`"wor"`) objects also contain:
#' \describe{
#'   \item{`$pik`}{Inclusion probability vector. For most methods this
#'   equals the true first-order probabilities. For `"sps"` and
#'   `"pareto"`, this is the target vector.}
#' }
#'
#' With-replacement (`"wr"`) objects also contain:
#' \describe{
#'   \item{`$prob`}{Per-draw selection probability vector.}
#'   \item{`$hits`}{Realized selection counts (integer vector, or
#'   `N x nrep` matrix when `nrep > 1`).}
#' }
#'
#' The class vector records the design family:
#' \describe{
#'   \item{`c("equal_prob", "wor", "sondage_sample")`}{Equal probability, without replacement.}
#'   \item{`c("equal_prob", "wr", "sondage_sample")`}{Equal probability, with replacement.}
#'   \item{`c("unequal_prob", "wor", "sondage_sample")`}{Unequal probability, without replacement (includes balanced).}
#'   \item{`c("unequal_prob", "wr", "sondage_sample")`}{Unequal probability, with replacement.}
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
