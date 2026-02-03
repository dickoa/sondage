#' Validate inclusion probability vector
#'
#' Checks that pik is a valid vector of inclusion probabilities.
#'
#' @param pik Vector to validate.
#' @param allow_zero If TRUE, allow pik values of exactly 0. Default TRUE.
#' @param allow_one If TRUE, allow pik values of exactly 1. Default TRUE.
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error.
#'
#' @keywords internal
#' @noRd
check_pik <- function(pik, allow_zero = TRUE, allow_one = TRUE) {
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
  invisible(TRUE)
}
