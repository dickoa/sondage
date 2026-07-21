#' Validate a non-empty numeric vector
#'
#' Performs the shape and missing-value checks shared by the package's
#' domain-specific vector validators.
#'
#' @param x Vector to validate.
#' @param name Parameter name for error messages.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
.check_numeric_vector <- function(x, name) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(sprintf("'%s' must be a numeric vector", name), call. = FALSE)
  }
  if (length(x) == 0L) {
    stop(sprintf("'%s' vector is empty", name), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("there are missing values in '%s'", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate inclusion probability vector
#'
#' Checks that pik is a valid vector of inclusion probabilities.
#'
#' @param pik Vector to validate.
#' @param allow_zero If TRUE, allow pik values of exactly 0. Default TRUE.
#' @param allow_one If TRUE, allow pik values of exactly 1. Default TRUE.
#' @param fixed_size If TRUE, ensure sample size is fixed. Default FALSE.
#' @param tol Tolerance for integer check. The default (`NULL`) uses a
#'   strict floating-point tolerance scaled by the number and magnitude
#'   of the terms: large enough to absorb accumulation error from
#'   `inclusion_prob()` at large N, but far below any statistically
#'   meaningful deviation. An exact fixed-size design cannot have
#'   `sum(pik)` away from an integer, so looser sums are rejected
#'   rather than silently rounded.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
.check_pik <- function(
  pik,
  allow_zero = TRUE,
  allow_one = TRUE,
  fixed_size = FALSE,
  tol = NULL
) {
  .check_numeric_vector(pik, "pik")
  if (any(pik < 0 | pik > 1)) {
    stop("inclusion probabilities must be between 0 and 1", call. = FALSE)
  }
  if (!allow_zero && any(pik == 0)) {
    stop("'pik' values of exactly 0 are not allowed", call. = FALSE)
  }
  if (!allow_one && any(pik == 1)) {
    stop("'pik' values of exactly 1 are not allowed", call. = FALSE)
  }
  if (fixed_size) {
    s <- sum(pik)
    if (is.null(tol)) {
      tol <- max(1e-10, 64 * .Machine$double.eps * (length(pik) + s))
    }
    if (abs(s - round(s)) > tol) {
      stop(
        sprintf(
          "sum(pik) = %.4g is not close to an integer",
          s
        ),
        call. = FALSE
      )
    }
    if (round(s) < 1L) {
      stop(
        sprintf(
          "sum(pik) = %.4g rounds to 0; need at least n = 1 for fixed-size sampling",
          s
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Validate permanent random numbers
#'
#' Checks that prn is a valid vector of permanent random numbers.
#'
#' @param prn Vector to validate.
#' @param N Expected length (population size).
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
.check_prn <- function(prn, N) {
  .check_numeric_vector(prn, "prn")
  if (length(prn) != N) {
    stop(
      sprintf("'prn' must have length %d (same as population size)", N),
      call. = FALSE
    )
  }
  if (any(prn <= 0) || any(prn >= 1)) {
    stop("'prn' values must be in the open interval (0, 1)", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check whether a value is a non-empty scalar character string
#'
#' This predicate is intentionally non-throwing so dispatchers can distinguish
#' possible registered method names from built-in method choices.
#'
#' @keywords internal
#' @noRd
.is_method_name <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    nzchar(x)
}

#' Validate a method registry name
#'
#' @keywords internal
#' @noRd
.check_method_name <- function(name) {
  if (!.is_method_name(name)) {
    stop("'name' must be a non-empty character string", call. = FALSE)
  }
  name
}

#' Match a documented character choice without exposing match.arg internals
#'
#' Partial matching is deliberately preserved for backward compatibility.
#'
#' @keywords internal
#' @noRd
.match_choice <- function(x, choices, name) {
  if (is.null(x) || identical(x, choices)) {
    return(choices[[1L]])
  }
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    match <- pmatch(x, choices, nomatch = 0L)
    if (!is.na(match) && match > 0L) {
      return(choices[[match]])
    }
  }

  quoted <- sprintf('"%s"', choices)
  alternatives <- if (length(quoted) == 1L) {
    quoted
  } else if (length(quoted) == 2L) {
    paste(quoted, collapse = " or ")
  } else {
    paste0(
      paste(quoted[-length(quoted)], collapse = ", "),
      ", or ",
      quoted[[length(quoted)]]
    )
  }
  stop(
    sprintf("'%s' must be one of %s", name, alternatives),
    call. = FALSE
  )
}

#' Validate a scalar logical flag
#'
#' @keywords internal
#' @noRd
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x) || !is.null(dim(x))) {
    stop(sprintf("'%s' must be TRUE or FALSE", name), call. = FALSE)
  }
  x
}

#' Validate a finite scalar number
#'
#' @param x Value to validate.
#' @param name Parameter name for error messages.
#' @param integer If TRUE, require a value close to an integer and return an
#'   integer. Otherwise, return a double.
#' @param tol Tolerance used when `integer = TRUE`.
#'
#' @keywords internal
#' @noRd
.check_number <- function(x, name, integer = FALSE, tol = 1e-4) {
  if (length(x) != 1L || !is.null(dim(x))) {
    stop(sprintf("'%s' must be a single numeric value", name), call. = FALSE)
  }
  if (is.na(x)) {
    stop(sprintf("'%s' must not be NA", name), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be a single numeric value", name), call. = FALSE)
  }
  if (!is.finite(x)) {
    stop(sprintf("'%s' must be finite", name), call. = FALSE)
  }
  if (!integer) {
    return(as.double(x))
  }

  r <- round(x)
  if (abs(x - r) > tol) {
    stop(
      sprintf("%s (%.4g) is not close to an integer", name, x),
      call. = FALSE
    )
  }
  as.integer(r)
}

#' Validate indices used for a sampled joint-probability submatrix
#'
#' @keywords internal
#' @noRd
.check_sample_idx <- function(sample_idx, N) {
  if (!is.numeric(sample_idx) || !is.null(dim(sample_idx))) {
    stop("'sample_idx' must be an integer vector", call. = FALSE)
  }
  if (anyNA(sample_idx)) {
    stop("'sample_idx' must not contain missing values", call. = FALSE)
  }
  if (any(!is.finite(sample_idx)) || any(sample_idx != floor(sample_idx))) {
    stop("'sample_idx' must contain whole numbers", call. = FALSE)
  }
  if (any(sample_idx < 1 | sample_idx > N)) {
    stop(
      sprintf("'sample_idx' values must be between 1 and length(pik) = %d", N),
      call. = FALSE
    )
  }
  if (anyDuplicated(sample_idx)) {
    stop("'sample_idx' must not contain duplicate indices", call. = FALSE)
  }
  as.integer(sample_idx)
}

#' Reject arguments that a built-in method does not consume
#'
#' The caller supplies only dot metadata, so checking names never forces the
#' promises in `...`. Keep the empty-dots guard at the call site to avoid a
#' helper call on performance-sensitive paths.
#'
#' @param n Number of arguments in `...`.
#' @param names Names of the arguments in `...`.
#' @param allowed Names accepted by the selected built-in method.
#'
#' @keywords internal
#' @noRd
.check_dots <- function(n, names, allowed = character()) {
  if (is.null(names)) {
    names <- rep.int("", n)
  }

  bad <- is.na(names) | !nzchar(names) | !(names %in% allowed)
  if (!any(bad)) {
    return(invisible(TRUE))
  }

  labels <- ifelse(
    is.na(names[bad]) | !nzchar(names[bad]),
    "<unnamed>",
    sprintf("'%s'", names[bad])
  )
  labels <- unique(labels)
  stop(
    sprintf(
      "unused argument%s in '...': %s",
      if (length(labels) == 1L) "" else "s",
      paste(labels, collapse = ", ")
    ),
    call. = FALSE
  )
}

#' Report that a method cannot use permanent random numbers
#'
#' @param method Selected method.
#' @param supported Built-in alternatives that support `prn`, when useful.
#'
#' @keywords internal
#' @noRd
.stop_unsupported_prn <- function(method, supported = character()) {
  message <- sprintf("method '%s' does not support 'prn'", method)
  if (length(supported) > 0L) {
    alternatives <- paste(sprintf("'%s'", supported), collapse = ", ")
    message <- paste0(
      message,
      "; permanent random numbers are supported by ",
      alternatives
    )
  }
  stop(message, call. = FALSE)
}

#' Validate replicate count and permanent-random-number compatibility
#'
#' @keywords internal
#' @noRd
.check_nrep_prn <- function(
  nrep,
  prn = NULL,
  method = NULL,
  supports_prn = TRUE,
  supported = character()
) {
  nrep <- .check_number(nrep, "nrep", integer = TRUE)
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }
  if (!is.null(prn) && !supports_prn) {
    .stop_unsupported_prn(method, supported)
  }
  if (!is.null(prn) && nrep > 1L) {
    stop(
      "prn and nrep > 1 cannot be used together. ",
      "Permanent random numbers produce identical samples across replicates. ",
      "Use a loop with different prn vectors for coordinated repeated sampling.",
      call. = FALSE
    )
  }
  nrep
}

#' Reject the removed design-modifying 'eps' argument
#'
#' The sampling methods formerly excluded units with pik <= eps and
#' force-selected units with pik >= 1 - eps. That silently changed the
#' design (and could break the fixed-size contract), so it was removed:
#' only exact 0 and exact 1 receive special treatment now.
#'
#' @keywords internal
#' @noRd
.stop_eps_removed <- function() {
  stop(
    "'eps' no longer modifies the design. Units are excluded or selected ",
    "with certainty only when pik is exactly 0 or exactly 1; all other ",
    "values are sampled as given.",
    call. = FALSE
  )
}

#' Validate numerical tolerance parameter
#'
#' @param eps Numeric tolerance value.
#' @param name Parameter name for error messages.
#'
#' @return The tolerance as a double.
#'
#' @keywords internal
#' @noRd
.check_eps <- function(eps, name = "eps") {
  eps <- .check_number(eps, name)
  if (eps <= 0 || eps >= 0.5) {
    stop(
      sprintf("'%s' must be in the open interval (0, 0.5)", name),
      call. = FALSE
    )
  }
  as.double(eps)
}

#' Validate expected hits vector
#'
#' Checks that hits is a valid vector of expected hits for WR sampling.
#'
#' @param hits Vector to validate.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
.check_hits <- function(hits) {
  .check_numeric_vector(hits, "hits")
  if (any(!is.finite(hits))) {
    stop("'hits' values must be finite (no Inf or NaN)", call. = FALSE)
  }
  if (any(hits < 0)) {
    stop("'hits' values must be non-negative", call. = FALSE)
  }
  if (sum(hits) == 0) {
    stop("sum of 'hits' must be positive", call. = FALSE)
  }
  invisible(TRUE)
}
