# Custom method registry
#
# An environment-based registry that lets users plug in custom unequal probability sampling
# methods (pure R) that flow through the existing dispatchers and generics.
.method_registry <- new.env(parent = emptyenv())

#' Register a Custom Unequal Probability Sampling Method
#'
#' Register a user-defined sampling method so it can be used through
#' [unequal_prob_wor()] or [unequal_prob_wr()] and their associated
#' generics.
#'
#' @param name A unique method name (character string). Must not
#'   collide with a built-in method name.
#' @param type `"wor"` (without replacement) or `"wr"` (with
#'   replacement).
#' @param sample_fn A function that draws a sample. See **Contracts**
#'   below.
#' @param joint_fn An optional function that computes joint inclusion
#'   probabilities (WOR) or joint expected hits (WR). If `NULL`,
#'   [joint_inclusion_prob()] / [joint_expected_hits()] will error for
#'   this method.
#' @param fixed_size Does this method always produce exactly `n` units?
#' @param supports_prn Does this method support permanent random
#'   numbers for sample coordination?
#'
#' @section Contracts:
#'
#' **`sample_fn(pik, n = NULL, prn = NULL, ...)`**
#' \describe{
#'   \item{`pik`}{Inclusion probabilities (WOR) or expected hits (WR),
#'     numeric vector of length N.}
#'   \item{`n`}{Target sample size (integer). For WOR this equals
#'     `round(sum(pik))`; for WR `round(sum(hits))`.}
#'   \item{`prn`}{Permanent random numbers (numeric vector length N,
#'     values in (0,1)), or `NULL`.}
#'   \item{Returns}{Integer vector of selected unit indices (1-based).
#'     For WOR: distinct indices of length `n` (fixed-size) or varying
#'     length (random-size). For WR: indices with possible repeats,
#'     of length `n`.}
#' }
#'
#' **`joint_fn(pik, sample_idx = NULL, ...)`** (optional)
#' \describe{
#'   \item{`pik`}{Same as above.}
#'   \item{`sample_idx`}{When non-NULL, an integer vector of sampled
#'     unit indices. Return only the submatrix for these units.}
#'   \item{Returns}{Symmetric matrix of joint inclusion probabilities
#'     (N x N when `sample_idx` is NULL, `length(sample_idx)` x
#'     `length(sample_idx)` otherwise).}
#' }
#'
#' @return Invisible `NULL`, called for its side effect.
#'
#' @seealso [registered_methods()], [unequal_prob_wor()],
#'   [unequal_prob_wr()]
#'
#' @examples
#' # Register a toy random sampler
#' my_sampler <- function(pik, n = NULL, prn = NULL, ...) {
#'   sample.int(length(pik), size = n, prob = pik)
#' }
#' register_method("toy", type = "wor", sample_fn = my_sampler)
#' s <- unequal_prob_wor(c(0.3, 0.3, 0.4), method = "toy")
#' s$method
#'
#' # Clean up
#' unregister_method("toy")
#'
#' @export
register_method <- function(
  name,
  type = c("wor", "wr"),
  sample_fn,
  joint_fn = NULL,
  fixed_size = TRUE,
  supports_prn = FALSE
) {
  if (!is.character(name) || length(name) != 1L || nchar(name) == 0L) {
    stop("'name' must be a non-empty character string", call. = FALSE)
  }
  type <- match.arg(type)
  if (!is.function(sample_fn)) {
    stop("'sample_fn' must be a function", call. = FALSE)
  }
  if (!is.null(joint_fn) && !is.function(joint_fn)) {
    stop("'joint_fn' must be a function or NULL", call. = FALSE)
  }
  if (!isTRUE(fixed_size) && !isFALSE(fixed_size)) {
    stop("'fixed_size' must be TRUE or FALSE", call. = FALSE)
  }
  if (!isTRUE(supports_prn) && !isFALSE(supports_prn)) {
    stop("'supports_prn' must be TRUE or FALSE", call. = FALSE)
  }

  # Reject collisions with built-in methods
  all_builtin <- unique(c(
    names(.wor_specs),
    names(.wr_specs),
    names(.ep_wor_specs),
    names(.ep_wr_specs)
  ))
  if (name %in% all_builtin) {
    stop(
      sprintf("'%s' is a built-in method and cannot be overridden", name),
      call. = FALSE
    )
  }

  .method_registry[[name]] <- list(
    name = name,
    type = type,
    sample_fn = sample_fn,
    joint_fn = joint_fn,
    fixed_size = fixed_size,
    supports_prn = supports_prn
  )

  invisible(NULL)
}

#' List Registered Custom Methods
#'
#' @return A character vector of registered method names (empty if none).
#'
#' @seealso [register_method()]
#'
#' @examples
#' registered_methods()
#'
#' @export
registered_methods <- function() {
  ls(.method_registry)
}

#' Check Whether a Method Is Registered
#'
#' @param name Method name (character string).
#' @return `TRUE` if the method has been registered via
#'   [register_method()], `FALSE` otherwise.
#'
#' @seealso [register_method()]
#'
#' @examples
#' is_registered_method("foo")
#'
#' @export
is_registered_method <- function(name) {
  exists(name, envir = .method_registry, inherits = FALSE)
}

#' Remove a Registered Method
#'
#' @param name Method name to unregister.
#' @return Invisible `TRUE` if the method was removed, `FALSE` if it
#'   was not registered.
#'
#' @seealso [register_method()]
#'
#' @examples
#' unregister_method("nonexistent")
#'
#' @export
unregister_method <- function(name) {
  if (is_registered_method(name)) {
    rm(list = name, envir = .method_registry)
    invisible(TRUE)
  } else {
    invisible(FALSE)
  }
}

#' @noRd
.dispatch_registered_wor <- function(pik, method, nrep, prn, ...) {
  reg <- .method_registry[[method]]
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn) && !reg$supports_prn) {
    warning(
      sprintf("prn is not used by method '%s' and will be ignored", method),
      call. = FALSE
    )
  }
  if (!is.null(prn) && nrep > 1L) {
    stop(
      "prn and nrep > 1 cannot be used together. ",
      "Permanent random numbers produce identical samples across replicates. ",
      "Use a loop with different prn vectors for coordinated repeated sampling.",
      call. = FALSE
    )
  }

  check_pik(pik, fixed_size = reg$fixed_size)
  N <- length(pik)
  n <- if (reg$fixed_size) as.integer(round(sum(pik))) else sum(pik)

  if (nrep == 1L) {
    idx <- as.integer(reg$sample_fn(pik, n = n, prn = prn, ...))
    .new_wor_sample(idx, pik, n, N, method, reg$fixed_size, "unequal_prob")
  } else if (reg$fixed_size) {
    n_int <- as.integer(round(sum(pik)))
    mat <- matrix(0L, n_int, nrep)
    for (i in seq_len(nrep)) {
      mat[, i] <- as.integer(reg$sample_fn(pik, n = n_int, prn = prn, ...))
    }
    .new_wor_sample(mat, pik, n_int, N, method, TRUE, "unequal_prob")
  } else {
    samples <- lapply(seq_len(nrep), function(i) {
      as.integer(reg$sample_fn(pik, n = n, prn = prn, ...))
    })
    .new_wor_sample(samples, pik, n, N, method, FALSE, "unequal_prob")
  }
}

#' @noRd
.dispatch_registered_wr <- function(hits, method, nrep, prn, ...) {
  reg <- .method_registry[[method]]
  nrep <- check_integer(nrep, "nrep")
  if (nrep < 1L) {
    stop("'nrep' must be at least 1", call. = FALSE)
  }

  if (!is.null(prn)) {
    warning(
      sprintf("prn is not used by method '%s' and will be ignored", method),
      call. = FALSE
    )
  }

  check_hits(hits)
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  prob <- hits / sum(hits)

  if (nrep == 1L) {
    idx <- as.integer(reg$sample_fn(hits, n = n, prn = prn, ...))
    realized_hits <- tabulate(idx, nbins = N)
    .new_wr_sample(
      idx,
      prob,
      realized_hits,
      n,
      N,
      method,
      reg$fixed_size,
      "unequal_prob"
    )
  } else {
    sample_mat <- matrix(0L, n, nrep)
    hits_mat <- matrix(0L, N, nrep)
    for (i in seq_len(nrep)) {
      idx <- as.integer(reg$sample_fn(hits, n = n, prn = prn, ...))
      sample_mat[, i] <- idx
      hits_mat[, i] <- tabulate(idx, nbins = N)
    }
    .new_wr_sample(
      sample_mat,
      prob,
      hits_mat,
      n,
      N,
      method,
      reg$fixed_size,
      "unequal_prob"
    )
  }
}

#' @noRd
.registered_joint_or_stop <- function(method, pik, sample_idx, quantity, ...) {
  if (!is_registered_method(method)) {
    .stop_no_joint(method, quantity)
  }
  reg <- .method_registry[[method]]
  if (is.null(reg$joint_fn)) {
    .stop_no_joint(method, quantity)
  }
  reg$joint_fn(pik, sample_idx = sample_idx, ...)
}
