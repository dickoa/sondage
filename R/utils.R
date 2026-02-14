#' Validate inclusion probability vector
#'
#' Checks that pik is a valid vector of inclusion probabilities.
#'
#' @param pik Vector to validate.
#' @param allow_zero If TRUE, allow pik values of exactly 0. Default TRUE.
#' @param allow_one If TRUE, allow pik values of exactly 1. Default TRUE.
#' @param tol Tolerance for integer check.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
check_pik <- function(
  pik,
  allow_zero = TRUE,
  allow_one = TRUE,
  fixed_size = FALSE,
  tol = 1e-4
) {
  if (!is.numeric(pik)) {
    stop("pik must be a numeric vector", call. = FALSE)
  }
  if (length(pik) == 0L) {
    stop("pik vector is empty", call. = FALSE)
  }
  if (anyNA(pik)) {
    stop("there are missing values in the pik vector", call. = FALSE)
  }
  if (any(pik < 0 | pik > 1)) {
    stop("inclusion probabilities must be between 0 and 1", call. = FALSE)
  }
  if (!allow_zero && any(pik == 0)) {
    stop("pik values of exactly 0 are not allowed", call. = FALSE)
  }
  if (!allow_one && any(pik == 1)) {
    stop("pik values of exactly 1 are not allowed", call. = FALSE)
  }
  if (fixed_size) {
    s <- sum(pik)
    if (abs(s - round(s)) > tol) {
      stop(
        sprintf(
          "sum(pik) = %.4g is not close to an integer",
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
check_prn <- function(prn, N) {
  if (!is.numeric(prn)) {
    stop("prn must be a numeric vector", call. = FALSE)
  }
  if (length(prn) != N) {
    stop(
      sprintf("prn must have length %d (same as population size)", N),
      call. = FALSE
    )
  }
  if (anyNA(prn)) {
    stop("there are missing values in prn", call. = FALSE)
  }
  if (any(prn <= 0) || any(prn >= 1)) {
    stop("prn values must be in the open interval (0, 1)", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate size measure vector
#'
#' Checks that x is a valid vector of non-negative size measures.
#'
#' @param x Vector to validate.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
check_mos <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be a numeric vector", call. = FALSE)
  }
  if (length(x) == 0L) {
    stop("x vector is empty", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("there are missing values in x", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x values must be finite (no Inf or NaN)", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("x values must be non-negative", call. = FALSE)
  }
  if (sum(x) == 0) {
    stop("sum of x must be positive", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that a value is close to an integer
#'
#' @param x Numeric value to check.
#' @param name Name of the parameter for the error message.
#' @param tol Tolerance for integer check.
#'
#' @return The value rounded to an integer.
#'
#' @keywords internal
#' @noRd
check_integer <- function(x, name = "n", tol = 1e-4) {
  if (is.na(x)) {
    stop(sprintf("%s must not be NA", name), call. = FALSE)
  }
  if (!is.finite(x)) {
    stop(sprintf("%s must be finite", name), call. = FALSE)
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
check_hits <- function(hits) {
  if (!is.numeric(hits)) {
    stop("hits must be a numeric vector", call. = FALSE)
  }
  if (length(hits) == 0L) {
    stop("hits vector is empty", call. = FALSE)
  }
  if (anyNA(hits)) {
    stop("there are missing values in hits", call. = FALSE)
  }
  if (any(!is.finite(hits))) {
    stop("hits values must be finite (no Inf or NaN)", call. = FALSE)
  }
  if (any(hits < 0)) {
    stop("hits values must be non-negative", call. = FALSE)
  }
  if (sum(hits) == 0) {
    stop("sum of hits must be positive", call. = FALSE)
  }
  invisible(TRUE)
}
