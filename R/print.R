#' @export
print.unequal_prob <- function(x, ...) {
  wr <- inherits(x, "wr")
  label <- if (wr) "Unequal prob WR" else "Unequal prob WOR"
  cat(sprintf("%s [%s] (n=%d, N=%d): ", label, x$method, x$n, x$N))
  if (length(x$sample) <= 20) {
    cat(x$sample, "\n")
  } else {
    cat(x$sample[1:10], "...\n")
  }
  invisible(x)
}

#' @export
print.equal_prob <- function(x, ...) {
  wr <- inherits(x, "wr")
  label <- if (wr) "Equal prob WR" else "Equal prob WOR"
  cat(sprintf("%s [%s] (n=%d, N=%d): ", label, x$method, x$n, x$N))
  if (length(x$sample) <= 20) {
    cat(x$sample, "\n")
  } else {
    cat(x$sample[1:10], "...\n")
  }
  invisible(x)
}
