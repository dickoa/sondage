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
  cps = list(
    fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer",
    draw = function(pik, prn = NULL) .Call(C_cps_single, pik)
  ),
  sampford = list(
    fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer",
    draw = function(pik, prn = NULL) .Call(C_up_sampford, pik)
  ),
  brewer = list(
    fixed_size = TRUE, prn = FALSE, variance_family = "pps_brewer",
    draw = function(pik, prn = NULL) .Call(C_up_brewer, pik)
  ),
  systematic = list(
    fixed_size = TRUE,
    prn = FALSE,
    variance_family = "pps_brewer",
    draw = function(pik, prn = NULL) .Call(C_up_systematic, pik)
  ),
  poisson = list(
    fixed_size = FALSE, prn = TRUE, variance_family = "poisson",
    draw = function(pik, prn = NULL) .Call(C_up_poisson, pik, prn)
  ),
  # The order-sampling pair honors pik to Rosen's documented
  # approximation rather than exactly; every other built-in is exact.
  sps = list(
    fixed_size = TRUE,
    prn = TRUE,
    variance_family = "pps_brewer",
    probabilities = "approximate",
    draw = function(pik, prn = NULL) .Call(C_up_sps, pik, prn)
  ),
  pareto = list(
    fixed_size = TRUE,
    prn = TRUE,
    variance_family = "pps_brewer",
    probabilities = "approximate",
    draw = function(pik, prn = NULL) .Call(C_up_pareto, pik, prn)
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

.dispatcher_specs <- list(
  equal_prob_wor = .ep_wor_specs,
  equal_prob_wr = .ep_wr_specs,
  unequal_prob_wor = .wor_specs,
  unequal_prob_wr = .wr_specs,
  balanced_wor = .balanced_specs
)

.dispatcher_types <- c(
  equal_prob_wor = "wor",
  equal_prob_wr = "wr",
  unequal_prob_wor = "wor",
  unequal_prob_wr = "wr",
  balanced_wor = "balanced"
)

#' @noRd
.format_dispatcher_choices <- function(dispatchers) {
  quoted <- paste0('"', dispatchers, '"')
  if (length(quoted) == 1L) {
    return(quoted)
  }
  if (length(quoted) == 2L) {
    return(paste(quoted, collapse = " or "))
  }
  paste0(
    paste(quoted[-length(quoted)], collapse = ", "),
    ", or ",
    quoted[[length(quoted)]]
  )
}

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
#' @param dispatcher Optional sampling entry point. One of
#'   `"equal_prob_wor"`, `"equal_prob_wr"`, `"unequal_prob_wor"`,
#'   `"unequal_prob_wr"`, or `"balanced_wor"`. Required when `name`
#'   is available through more than one entry point, as with `"srs"`
#'   and `"systematic"`. Values must match exactly.
#'
#' @return A list with elements `dispatcher` (the sampling entry point),
#'   `type` (`"wor"`, `"wr"`, or `"balanced"`), `fixed_size` (logical),
#'   `variance_family` (one of
#'   `"srs"`, `"pps_brewer"`, `"poisson"`, `"wr"`, `"unsupported"`, or
#'   `NULL` for a registered method that did not declare one; see
#'   [register_method()]), `supports_prn` (logical), `supports_aux`
#'   (logical), `supports_strata` (logical), `supports_spread`
#'   (logical), and `probabilities` (where the method sits in the
#'   first-order probability taxonomy: `"exact"` for every built-in
#'   except `"sps"` and `"pareto"`, which honor `pik` to a documented
#'   approximation and report `"approximate"`; for a registered
#'   method the declared tier, `"unknown"` when the author did not
#'   establish one), plus `sample_fn` and `joint_fn` (the registered
#'   implementation functions for a registered method, `NULL` for
#'   built-ins, whose implementations are internal dispatch paths).
#'   Returns `NULL` if the method is unknown. The
#'   aux/strata/spread capabilities are only `TRUE` for balanced
#'   methods. An ambiguous built-in name without `dispatcher` is an
#'   error rather than silently selecting one variant.
#'
#' @seealso [register_method()], [registered_methods()]
#'
#' @examples
#' method_spec("brewer")
#' method_spec("cube")
#' method_spec("srs", dispatcher = "equal_prob_wr")
#' method_spec("systematic", dispatcher = "equal_prob_wor")
#' method_spec("nonexistent")
#'
#' @export
method_spec <- function(name, dispatcher = NULL) {
  name <- .check_method_name(name)
  dispatchers <- names(.dispatcher_specs)
  if (
    !is.null(dispatcher) &&
      (!is.character(dispatcher) ||
        length(dispatcher) != 1L ||
        !is.null(dim(dispatcher)) ||
        is.na(dispatcher) ||
        !(dispatcher %in% dispatchers))
  ) {
    stop(
      "'dispatcher' must be NULL or exactly one of ",
      .format_dispatcher_choices(dispatchers),
      call. = FALSE
    )
  }

  if (is_registered_method(name)) {
    reg <- .method_registry[[name]]
    registered_dispatcher <- switch(
      reg$type,
      wor = "unequal_prob_wor",
      wr = "unequal_prob_wr",
      balanced = "balanced_wor"
    )
    if (!is.null(dispatcher) && dispatcher != registered_dispatcher) {
      stop(
        sprintf(
          "registered method '%s' uses dispatcher \"%s\", not \"%s\"",
          name,
          registered_dispatcher,
          dispatcher
        ),
        call. = FALSE
      )
    }
    return(list(
      dispatcher = registered_dispatcher,
      type = reg$type,
      fixed_size = reg$fixed_size,
      # NULL when undeclared: list(x = NULL) keeps the element, so the
      # field is always present in the returned spec
      variance_family = reg$variance_family,
      supports_prn = reg$supports_prn,
      supports_aux = isTRUE(reg$supports_aux),
      supports_strata = isTRUE(reg$supports_strata),
      supports_spread = isTRUE(reg$supports_spread),
      probabilities = reg$probabilities,
      sample_fn = reg$sample_fn,
      joint_fn = reg$joint_fn
    ))
  }

  available <- dispatchers[vapply(
    .dispatcher_specs,
    function(specs) !is.null(specs[[name]]),
    logical(1)
  )]
  if (length(available) == 0L) {
    return(NULL)
  }
  if (is.null(dispatcher) && length(available) > 1L) {
    stop(
      sprintf(
        "method '%s' is ambiguous; supply 'dispatcher' as %s",
        name,
        .format_dispatcher_choices(available)
      ),
      call. = FALSE
    )
  }
  if (!is.null(dispatcher) && !(dispatcher %in% available)) {
    stop(
      sprintf(
        "method '%s' is not available through dispatcher \"%s\"; use %s",
        name,
        dispatcher,
        .format_dispatcher_choices(available)
      ),
      call. = FALSE
    )
  }
  if (is.null(dispatcher)) {
    dispatcher <- available[[1L]]
  }

  spec <- .dispatcher_specs[[dispatcher]][[name]]
  list(
    dispatcher = dispatcher,
    type = unname(.dispatcher_types[[dispatcher]]),
    fixed_size = isTRUE(spec$fixed_size),
    variance_family = spec$variance_family,
    supports_prn = isTRUE(spec$prn),
    supports_aux = isTRUE(spec$aux),
    supports_strata = isTRUE(spec$strata),
    supports_spread = isTRUE(spec$spread),
    # Exact is the built-in norm; only the order-sampling pair
    # (sps, pareto) carries "approximate" in its spec entry.
    probabilities = if (is.null(spec$probabilities)) {
      "exact"
    } else {
      spec$probabilities
    },
    # Built-in implementations are internal dispatch paths, not
    # registry entries; only registered methods expose functions.
    sample_fn = NULL,
    joint_fn = NULL
  )
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
