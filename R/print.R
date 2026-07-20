#' Print Sampling Design Objects
#'
#' @param x A sampling design object of class `"sondage_sample"`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return `invisible(x)`.
#'
#' @seealso [sondage_sample] for the documented object structure and fields.
#'
#' @name print.sondage_sample
NULL

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
    cat(sprintf(
      "%d replicates\n  rep 1 (realized n=%d): ",
      nrep,
      length(s[[1]])
    ))
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

#' @noRd
.print_design <- function(x, label) {
  if (!x$fixed_size && !is.list(x$sample)) {
    size <- sprintf(
      "expected n=%g, realized n=%d, N=%d",
      x$n,
      length(x$sample),
      x$N
    )
  } else if (!x$fixed_size) {
    size <- sprintf("expected n=%g, N=%d", x$n, x$N)
  } else {
    size <- sprintf("n=%g, N=%d", x$n, x$N)
  }

  cat(sprintf("%s [%s] (%s): ", label, x$method, size))
  .print_sample(x$sample)
  invisible(x)
}

#' @rdname print.sondage_sample
#' @export
print.unequal_prob <- function(x, ...) {
  label <- if (inherits(x, "balanced")) {
    "Balanced WOR"
  } else if (inherits(x, "wr")) {
    "Unequal prob WR"
  } else {
    "Unequal prob WOR"
  }
  .print_design(x, label)
}

#' @rdname print.sondage_sample
#' @export
print.equal_prob <- function(x, ...) {
  label <- if (inherits(x, "wr")) "Equal prob WR" else "Equal prob WOR"
  .print_design(x, label)
}
