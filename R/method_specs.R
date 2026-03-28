# Built-in method spec tables and internal helpers
#
# One source of truth for built-in method capabilities, split by dispatcher
# context to avoid cross-contamination (e.g. bernoulli never appears in
# unequal-prob PRN checks). These tables, the resolver, and the hook helpers
# are the future insertion points for custom method registration.
.wor_specs <- list(
  cps = list(fixed_size = TRUE, prn = FALSE),
  brewer = list(fixed_size = TRUE, prn = FALSE),
  systematic = list(fixed_size = TRUE, prn = FALSE),
  poisson = list(fixed_size = FALSE, prn = TRUE),
  sps = list(fixed_size = TRUE, prn = TRUE),
  pareto = list(fixed_size = TRUE, prn = TRUE)
)

.wr_specs <- list(
  chromy = list(fixed_size = TRUE, prn = FALSE),
  multinomial = list(fixed_size = TRUE, prn = FALSE)
)

.ep_wor_specs <- list(
  srs = list(fixed_size = TRUE, prn = FALSE),
  systematic = list(fixed_size = TRUE, prn = FALSE),
  bernoulli = list(fixed_size = FALSE, prn = TRUE)
)

.ep_wr_specs <- list(
  srs = list(fixed_size = TRUE, prn = FALSE)
)

.balanced_specs <- list(
  cube = list(fixed_size = TRUE, prn = FALSE)
)

# Methods using the high-entropy JIP approximation (Brewer & Donadio, 2003).
# Spans wor + balanced; used in joint_inclusion_prob.wor marginal defect check.
.he_jip_methods <- c("brewer", "sps", "pareto", "cube")

# Methods with C-level batch optimisation (workspace reuse).
.batch_optimised_methods <- c("cps")

#' Look up the spec for a built-in method.
#'
#' @param method Method name string.
#' @param context One of `"wor"`, `"wr"`, `"ep_wor"`, `"ep_wr"`.
#' @return The spec list, or `NULL` if `method` is not in the table.
#' @noRd
.get_builtin_spec <- function(method, context) {
  specs <- switch(
    context,
    wor = .wor_specs,
    wr = .wr_specs,
    ep_wor = .ep_wor_specs,
    ep_wr = .ep_wr_specs,
    stop("unknown context: ", context, call. = FALSE)
  )
  specs[[method]]
}

#' Query Method Metadata
#'
#' Return the capabilities of a sampling method. Works for built-in
#' methods and methods added via [register_method()].
#'
#' @param name Method name (character string), as used by the sondage
#'   dispatchers (e.g. `"brewer"`, `"cube"`, `"srs"`).
#'
#' @return A list with elements `type` (`"wor"` or `"wr"`),
#'   `fixed_size` (logical), and `supports_prn` (logical), or `NULL`
#'   if the method is unknown.
#'
#' @seealso [register_method()], [registered_methods()]
#'
#' @examples
#' method_spec("brewer")
#' method_spec("cube")
#' method_spec("nonexistent")
#'
#' @export
method_spec <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("'name' must be a single character string", call. = FALSE)
  }

  if (is_registered_method(name)) {
    reg <- .method_registry[[name]]
    return(list(
      type = reg$type,
      fixed_size = reg$fixed_size,
      supports_prn = reg$supports_prn
    ))
  }

  all_specs <- list(
    wor = .wor_specs,
    wr = .wr_specs,
    ep_wor = .ep_wor_specs,
    ep_wr = .ep_wr_specs,
    balanced = .balanced_specs
  )

  for (ctx in names(all_specs)) {
    spec <- all_specs[[ctx]][[name]]
    if (!is.null(spec)) {
      type <- if (ctx %in% c("wr", "ep_wr")) "wr" else "wor"
      return(list(
        type = type,
        fixed_size = isTRUE(spec$fixed_size),
        supports_prn = isTRUE(spec$prn)
      ))
    }
  }

  NULL
}

#' @noRd
.method_supports_prn <- function(method, context) {
  spec <- .get_builtin_spec(method, context)
  if (!is.null(spec)) return(isTRUE(spec$prn))
  reg <- .method_registry[[method]]
  if (!is.null(reg)) return(isTRUE(reg$supports_prn))
  FALSE
}

#' @noRd
.method_is_fixed_size <- function(method, context) {
  spec <- .get_builtin_spec(method, context)
  if (!is.null(spec)) return(isTRUE(spec$fixed_size))
  reg <- .method_registry[[method]]
  if (!is.null(reg)) return(isTRUE(reg$fixed_size))
  FALSE
}

# For sampling dispatch (dispatchers + batch helpers).
# New error path, switch currently returns NULL silently; match.arg fires
# first in practice, so this is unreachable for built-in methods.
#' @noRd
.stop_unknown_method <- function(method) {
  stop(
    sprintf("unknown sampling method '%s'", method),
    call. = FALSE
  )
}

# For joint prob / expected hits dispatch (generics.R).
# IMPORTANT: message format is byte-compatible with the existing inline
# stop() calls:
#   "joint_inclusion_prob not implemented for method '%s'"
#   "joint_expected_hits not implemented for method '%s'"
#' @noRd
.stop_no_joint <- function(method, quantity) {
  stop(
    sprintf("%s not implemented for method '%s'", quantity, method),
    call. = FALSE
  )
}
