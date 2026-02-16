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
