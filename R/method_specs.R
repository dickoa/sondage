# Built-in method spec tables and internal helpers
#
# One source of truth for built-in method capabilities, split by dispatcher
# context to avoid cross-contamination (e.g. bernoulli never appears in
# unequal-prob PRN checks). These tables, the resolver, and the hook helpers
# are the future insertion points for custom method registration.
# variance_family names the design-based variance treatment a
# downstream consumer (e.g. a survey-export package) should apply.
# It names the estimator family, not a guarantee of exactness:
# equal-probability systematic gets the SRS treatment even though
# that variance is an approximation for it.
.variance_families <- c("srs", "pps_brewer", "poisson", "wr", "unsupported")

.wor_specs <- list(
  cps = list(fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer"),
  sampford = list(fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer"),
  brewer = list(fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer"),
  systematic = list(
    fixed_size = TRUE,
    prn = FALSE,
    variance_family = "pps_brewer"
  ),
  poisson = list(fixed_size = FALSE, prn = TRUE, variance_family = "poisson"),
  # The order-sampling pair honors pik to Rosen's documented
  # approximation rather than exactly; every other built-in is exact.
  sps = list(
    fixed_size = TRUE,
    prn = TRUE,
    variance_family = "pps_brewer",
    probabilities = "approximate"
  ),
  pareto = list(
    fixed_size = TRUE,
    prn = TRUE,
    variance_family = "pps_brewer",
    probabilities = "approximate"
  )
)

.wr_specs <- list(
  chromy = list(fixed_size = TRUE, prn = FALSE, variance_family = "wr"),
  multinomial = list(fixed_size = TRUE, prn = FALSE, variance_family = "wr")
)

.ep_wor_specs <- list(
  srs = list(fixed_size = TRUE, prn = FALSE, variance_family = "srs"),
  systematic = list(fixed_size = TRUE, prn = FALSE, variance_family = "srs"),
  bernoulli = list(fixed_size = FALSE, prn = TRUE, variance_family = "poisson")
)

.ep_wr_specs <- list(
  srs = list(fixed_size = TRUE, prn = FALSE, variance_family = "wr")
)

.balanced_specs <- list(
  cube = list(
    fixed_size = TRUE,
    prn = FALSE,
    aux = TRUE,
    strata = TRUE,
    spread = FALSE,
    variance_family = "pps_brewer"
  ),
  lpm2 = list(
    fixed_size = TRUE,
    prn = FALSE,
    aux = FALSE,
    strata = FALSE,
    spread = TRUE,
    variance_family = "unsupported"
  ),
  scps = list(
    fixed_size = TRUE,
    prn = FALSE,
    aux = FALSE,
    strata = FALSE,
    spread = TRUE,
    variance_family = "unsupported"
  )
)

# Methods using the high-entropy JIP approximation (Brewer & Donadio, 2003).
# Spans wor + balanced; used in joint_inclusion_prob.wor marginal defect check.
.he_jip_methods <- c("brewer", "sps", "pareto", "cube")

# Methods with C-level batch optimisation (workspace reuse).
.batch_optimised_methods <- c("cps")

#' Look up the spec for a built-in method.
#'
#' @param method Method name string.
#' @param context One of `"wor"`, `"wr"`, `"ep_wor"`, `"ep_wr"`,
#'   `"balanced"`.
#' @return The spec list, or `NULL` if `method` is not in the table.
#' @noRd
.get_builtin_spec <- function(method, context) {
  specs <- switch(
    context,
    wor = .wor_specs,
    wr = .wr_specs,
    ep_wor = .ep_wor_specs,
    ep_wr = .ep_wr_specs,
    balanced = .balanced_specs,
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
#' @return A list with elements `type` (`"wor"`, `"wr"`, or
#'   `"balanced"`), `fixed_size` (logical), `variance_family` (one of
#'   `"srs"`, `"pps_brewer"`, `"poisson"`, `"wr"`, `"unsupported"`, or
#'   `NULL` for a registered method that did not declare one; see
#'   [register_method()]), `probabilities` (where the method sits in the
#'   first-order probability taxonomy: `"exact"` for every built-in
#'   except `"sps"` and `"pareto"`, which honor `pik` to a documented
#'   approximation and report `"approximate"`; for a registered
#'   method the declared tier, `"unknown"` when the author did not
#'   establish one), `supports_prn` (logical), `supports_aux`
#'   (logical), `supports_strata` (logical), and `supports_spread`
#'   (logical), or `NULL` if the method is unknown. The
#'   aux/strata/spread capabilities are only `TRUE` for balanced
#'   methods. For the two names shared by an equal- and an
#'   unequal-probability built-in, the lookup resolves as it does for
#'   `type`: `"systematic"` reports the unequal-probability variant
#'   and `"srs"` the without-replacement variant.
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
      # NULL when undeclared: list(x = NULL) keeps the element, so the
      # field is always present in the returned spec
      variance_family = reg$variance_family,
      probabilities = reg$probabilities,
      supports_prn = reg$supports_prn,
      supports_aux = isTRUE(reg$supports_aux),
      supports_strata = isTRUE(reg$supports_strata),
      supports_spread = isTRUE(reg$supports_spread)
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
      type <- switch(
        ctx,
        wr = ,
        ep_wr = "wr",
        balanced = "balanced",
        "wor"
      )
      return(list(
        type = type,
        fixed_size = isTRUE(spec$fixed_size),
        variance_family = spec$variance_family,
        # Exact is the built-in norm; only the order-sampling pair
        # (sps, pareto) carries "approximate" in its spec entry.
        probabilities = if (is.null(spec$probabilities)) {
          "exact"
        } else {
          spec$probabilities
        },
        supports_prn = isTRUE(spec$prn),
        supports_aux = isTRUE(spec$aux),
        supports_strata = isTRUE(spec$strata),
        supports_spread = isTRUE(spec$spread)
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
