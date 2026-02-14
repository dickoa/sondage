#' @keywords internal
#' @noRd
.new_wor_sample <- function(sample, pik, n, N, method, fixed_size, prob_class) {
  structure(
    list(
      sample = as.integer(sample),
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
  structure(
    list(
      sample = as.integer(sample),
      prob = prob,
      hits = as.integer(hits),
      n = n,
      N = as.integer(N),
      method = method,
      fixed_size = fixed_size
    ),
    class = c(prob_class, "wr", "sondage_sample")
  )
}
