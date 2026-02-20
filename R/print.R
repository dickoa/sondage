#' @noRd
.print_sample <- function(s) {
  if (is.matrix(s)) {
    nrep <- ncol(s)
    n <- nrow(s)
    cat(sprintf("%d replicates\n  rep 1: ", nrep))
    if (n <= 20) {
      cat(s[, 1], "\n")
    } else {
      cat(s[1:10, 1], "...\n")
    }
  } else if (is.list(s)) {
    nrep <- length(s)
    cat(sprintf("%d replicates\n  rep 1 (n=%d): ", nrep, length(s[[1]])))
    if (length(s[[1]]) <= 20) {
      cat(s[[1]], "\n")
    } else {
      cat(s[[1]][1:10], "...\n")
    }
  } else {
    if (length(s) <= 20) {
      cat(s, "\n")
    } else {
      cat(s[1:10], "...\n")
    }
  }
}

#' @export
print.unequal_prob <- function(x, ...) {
  wr <- inherits(x, "wr")
  label <- if (wr) "Unequal prob WR" else "Unequal prob WOR"
  cat(sprintf("%s [%s] (n=%g, N=%d): ", label, x$method, x$n, x$N))
  .print_sample(x$sample)
  invisible(x)
}

#' @export
print.equal_prob <- function(x, ...) {
  wr <- inherits(x, "wr")
  label <- if (wr) "Equal prob WR" else "Equal prob WOR"
  cat(sprintf("%s [%s] (n=%g, N=%d): ", label, x$method, x$n, x$N))
  .print_sample(x$sample)
  invisible(x)
}
