#' @noRd
inclusion_probs <- function(a, n) {
  storage.mode(a) <- "numeric"
  b <- a < 0
  if (any(b)) {
    warning("there are ", sum(b), " negative value(s) shifted to zero")
    a[b] <- 0
  }
  .Call(C_inclusion_probs, a, n)
}

#' @noRd
up_brewer <- function(pik, eps = sqrt(.Machine$double.eps)) {
  if (any(is.na(pik)))
    stop("there are missing values in the pik vector")
  storage.mode(pik) <- storage.mode(eps) <- "numeric"
  r <- .Call(C_up_brewer, pik, eps)
  r
}

#' @noRd
sample_pps <- function(n, size, prob, tolerance = sqrt(.Machine$double.eps)) {
  s <- sum(prob)
  sums_to_one <- isTRUE(all.equal(s, 1, tolerance = tolerance))
  sums_to_int <- isTRUE(all.equal(s, round(s), tolerance = tolerance))
  if (is.null(size)) {
    if(!sums_to_int)
      stop("sum(prob) must be an integer")
    size <- round(s)
  } else {
    size_is_sum <- isTRUE(all.equal(size, sum(prob), tolerance = tolerance))
    size_is_int <- isTRUE(all.equal(size, round(size), tolerance = tolerance))
    if (!size_is_int)
      stop("size must be NULL or an integer")
    if (sums_to_one && !size_is_sum) {
      warning("rescaling prob, which changes inclusion probabilities")
      prob <- inclusion_probs(prob * size, size)
    } else if (sums_to_int && !size_is_sum) {
      warning("sum(prob) is not equal to size or 1, rescaling")
      prob <- inclusion_probs(prob * size/s, size)
    }
  }
  up_brewer(prob)
}
