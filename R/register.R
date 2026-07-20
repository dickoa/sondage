# Custom method registry
#
# An environment-based registry that lets users plug in custom sampling
# methods (pure R) that flow through the existing dispatchers and generics.
.method_registry <- new.env(parent = emptyenv())

#' Register a Custom Sampling Method
#'
#' Register a user-defined sampling method so it can be used through
#' [unequal_prob_wor()], [unequal_prob_wr()], or [balanced_wor()] and
#' their associated generics.
#'
#' Registrations last for the current R session. Registering an existing
#' custom method name is an error; call [unregister_method()] first when
#' deliberately replacing a method.
#'
#' @param name A unique method name (character string). Must not
#'   collide with a built-in method name.
#' @param type `"wor"` (without replacement), `"wr"` (with
#'   replacement), or `"balanced"` (balanced without replacement,
#'   dispatched through [balanced_wor()]).
#' @param sample_fn A function that draws a sample. See **Contracts**
#'   below.
#' @param joint_fn An optional function that computes joint inclusion
#'   probabilities (WOR) or joint expected hits (WR). If `NULL`,
#'   [joint_inclusion_prob()] / [joint_expected_hits()] will error for
#'   this method.
#' @param fixed_size Does this method always produce exactly `n` units?
#'   Must be `TRUE` for `type = "wr"`.
#' @param variance_family Optional declaration of how design-based
#'   variance should be estimated for this method, for downstream
#'   packages that export designs for variance estimation. One of
#'   `"srs"`, `"pps_brewer"`, `"poisson"`, `"wr"`, or
#'   `"unsupported"`; see **Variance families** below. `NULL` (the
#'   default) means undeclared: consumers fall back on inferring a
#'   treatment from `type` and `fixed_size`.
#' @param supports_prn Does this method support permanent random
#'   numbers for sample coordination? Set to `TRUE` or `FALSE` for
#'   `"wor"` and `"wr"` methods. `NULL` (the default) resolves to
#'   `FALSE`. Leave `NULL` for `"balanced"` methods; [balanced_wor()]
#'   has no `prn` argument.
#' @param supports_aux Does this method use auxiliary balancing
#'   variables? Set to `TRUE` or `FALSE` for `type = "balanced"`.
#'   `NULL` (the default) resolves to `TRUE` for balanced methods and
#'   `FALSE` otherwise. Set it to `FALSE` for spread-only (spatially
#'   balanced) methods such as the local pivotal method, so that passing
#'   `aux` to [balanced_wor()] is an error instead of being silently
#'   ignored. Leave `NULL` for `"wor"` and `"wr"` methods.
#' @param supports_strata Does this method support stratified balanced
#'   sampling? Set to `TRUE` or `FALSE` for `type = "balanced"`.
#'   `NULL` (the default) resolves to `FALSE`. When `FALSE`, passing
#'   `strata` to [balanced_wor()] with this method is an error, and
#'   `sample_fn` does not need a `strata` argument. Leave `NULL` for
#'   `"wor"` and `"wr"` methods.
#' @param supports_spread Does this method support spatial spreading
#'   (well-spread / spatially balanced sampling)? Set to `TRUE` or
#'   `FALSE` for `type = "balanced"`. `NULL` (the default) resolves to
#'   `FALSE`. When `FALSE`, passing `spread` to [balanced_wor()] with
#'   this method is an error, and `sample_fn` does not need a `spread`
#'   argument. Leave `NULL` for `"wor"` and `"wr"` methods.
#' @param probabilities Where the method sits in the first-order
#'   probability taxonomy, for downstream packages that weight or
#'   record per-unit selection probabilities:
#'   - `"exact"`: the true first-order inclusion probabilities
#'     (types `"wor"` and `"balanced"`) or expected hits (type
#'     `"wr"`) equal the `pik` or `hits` vector passed to `sample_fn`,
#'     as for Sampford or the cube method.
#'   - `"approximate"`: `pik` or `hits` is the method's first-order target,
#'     achieved up to a documented approximation, as for Pareto or
#'     sequential Poisson order sampling. Design weights `1/pik`
#'     remain standard practice.
#'   - `"unknown"` (the default): `pik` or `hits` is an input weight or
#'     preference only, so the corresponding first-order quantities
#'     are not known. The toy WOR sampler in the examples is such a method:
#'     successive sampling with `prob = pik` does not yield inclusion
#'     probabilities equal to `pik`, so its design weights `1/pik`
#'     would be systematically biased, not merely noisy. (The same
#'     draw \emph{with} replacement does honor expected hits, so a
#'     multinomial-style `type = "wr"` method declares `"exact"`.)
#'
#'   The default is deliberately strict: if you have not established
#'   which tier your method is in, its selection probabilities are
#'   unknown, and downstream packages that weight estimation by
#'   `1/pik` may refuse to draw with it rather than produce biased
#'   weights. Sampling through sondage itself is unaffected: the
#'   declaration describes the method, it never disables it. Like
#'   `variance_family`, the declaration is an assertion by the method
#'   author that sondage cannot verify; the package vignette shows
#'   how to check a first-order contract by simulation.
#'
#' @section Contracts:
#'
#' **`sample_fn(pik, n = NULL, prn = NULL, ...)`** (type `"wor"`)
#' \describe{
#'   \item{`pik`}{Inclusion probabilities, numeric vector of length N.}
#'   \item{`n`}{Target sample size, equal to `round(sum(pik))`.}
#'   \item{`prn`}{Permanent random numbers (numeric vector length N,
#'     values in (0,1)), or `NULL`. Supplying `prn` is an error when the
#'     method is registered with `supports_prn = FALSE`, and `sample_fn`
#'     is not called.}
#'   \item{Returns}{Integer vector of selected unit indices (1-based).
#'     Indices are distinct and have length `n` for fixed-size methods,
#'     or varying length for random-size methods. The dispatcher validates
#'     the type, range, size, and replacement rules before constructing
#'     the sample object.}
#' }
#'
#' **`sample_fn(hits, n = NULL, prn = NULL, ...)`** (type `"wr"`)
#' \describe{
#'   \item{`hits`}{Expected hits, numeric vector of length N.}
#'   \item{`n`}{Target sample size, equal to `round(sum(hits))`.}
#'   \item{`prn`}{Permanent random numbers (numeric vector length N,
#'     values in (0,1)), or `NULL`. Supplying `prn` is an error when the
#'     method is registered with `supports_prn = FALSE`, and `sample_fn`
#'     is not called.}
#'   \item{Returns}{Integer vector of `n` selected unit indices (1-based),
#'     with possible repeats. The dispatcher validates the type, range,
#'     size, and replacement rules before constructing the sample object.}
#' }
#'
#' **`sample_fn(pik, n = NULL, aux = NULL, ...)`** (type `"balanced"`)
#' \describe{
#'   \item{`pik`}{Inclusion probabilities, numeric vector of length N.}
#'   \item{`n`}{Target sample size (integer when `fixed_size`,
#'     otherwise `sum(pik)`).}
#'   \item{`aux`}{Auxiliary balancing matrix (N x p, double), or
#'     `NULL`. Passed through as supplied to [balanced_wor()] after
#'     validation; the sample-size constraint is \strong{not}
#'     prepended, so add it yourself if your algorithm needs it
#'     (e.g. `cbind(pik, aux)`). Methods registered with
#'     `supports_aux = FALSE` always receive `aux = NULL`.}
#'   \item{`strata`}{Only when registered with
#'     `supports_strata = TRUE` and the caller supplies `strata`:
#'     an integer vector (length N) of dense stratum labels `1:H`.
#'     Declare it as `strata = NULL` in your function signature.}
#'   \item{`spread`}{Only when registered with
#'     `supports_spread = TRUE` and the caller supplies `spread`:
#'     a numeric matrix (N x d, double) of spatial coordinates (or
#'     other spreading variables). Declare it as `spread = NULL` in
#'     your function signature.}
#'   \item{Returns}{Integer vector of distinct selected unit indices
#'     (1-based). The dispatcher validates the returned indices.}
#' }
#'
#' **`joint_fn(pik, sample_idx = NULL, ...)`** (optional; types `"wor"`
#' and `"balanced"`)
#' \describe{
#'   \item{`pik`}{Inclusion probabilities, numeric vector of length N.}
#'   \item{`eps`}{If the function explicitly declares an `eps` formal,
#'     it receives the value supplied to [joint_inclusion_prob()].}
#'   \item{`sample_idx`}{When non-NULL, an integer vector of sampled
#'     unit indices. Return only the submatrix for these units.}
#'   \item{Returns}{Symmetric matrix of joint inclusion probabilities
#'     (N x N when `sample_idx` is NULL, `length(sample_idx)` x
#'     `length(sample_idx)` otherwise). The dispatcher validates that
#'     the matrix has the required dimensions and contains finite,
#'     symmetric numeric values.}
#' }
#'
#' **`joint_fn(hits, sample_idx = NULL, ...)`** (optional; type `"wr"`)
#' \describe{
#'   \item{`hits`}{Expected hits, numeric vector of length N, renormalized
#'     to sum exactly to the integer sample size `n`.}
#'   \item{`nsim`}{If the function explicitly declares an `nsim` formal,
#'     it receives the value supplied to [joint_expected_hits()].}
#'   \item{`sample_idx`}{When non-NULL, an integer vector of distinct
#'     sampled unit indices. Return only the submatrix for these units.}
#'   \item{Returns}{Symmetric matrix of joint expected hits
#'     (N x N when `sample_idx` is NULL, `length(sample_idx)` x
#'     `length(sample_idx)` otherwise). The dispatcher validates that
#'     the matrix has the required dimensions and contains finite,
#'     symmetric numeric values.}
#' }
#'
#' @section Variance families:
#'
#' `variance_family` names the estimator treatment a variance consumer
#' (such as a survey-export package) should apply to samples drawn
#' with this method. Selection metadata alone cannot determine it: a
#' fixed-size WOR method may need Brewer's unequal-probability
#' approximation or an SRS-style variance, and for a random-size WOR
#' method no safe inference exists at all, because a Poisson-type
#' method (independent selections) and a correlated random-size
#' scheme need different estimators.
#'
#' \describe{
#'   \item{`"srs"`}{Equal-probability fixed-size WOR. SRS-style
#'     variance with a finite population correction. Requires
#'     `fixed_size = TRUE`.}
#'   \item{`"pps_brewer"`}{Fixed-size unequal-probability WOR.
#'     Brewer's approximation from the marginal inclusion
#'     probabilities. Requires `fixed_size = TRUE`.}
#'   \item{`"poisson"`}{Random-size WOR with \strong{independent}
#'     selections (Poisson-type). Exact Poisson linearization.
#'     Requires `type = "wor"` and `fixed_size = FALSE`.}
#'   \item{`"wr"`}{With-replacement (or minimum-replacement)
#'     selection. Hansen-Hurwitz variance, no finite population
#'     correction. Requires `type = "wr"`.}
#'   \item{`"unsupported"`}{No linearization treatment is valid;
#'     consumers should refuse to linearize and point to replicate
#'     methods instead. The correct declaration for correlated
#'     random-size schemes. Always allowed.}
#' }
#'
#' Balanced methods allow only `"pps_brewer"` or `"unsupported"`:
#' balancing constraints couple selections across units, so
#' `"poisson"` can never hold for them.
#'
#' The declaration is an assertion by the method author, not
#' something sondage can verify; `"poisson"` in particular asserts
#' that units are selected independently, and a wrong declaration
#' produces silently wrong variance estimates for every user of the
#' method. The package vignette
#' (`vignette("custom-methods", package = "sondage")`) shows how to
#' check a declared family by simulation.
#'
#' @return Invisible `NULL`, called for its side effect.
#'
#' @seealso [registered_methods()], [method_spec()],
#'   [unequal_prob_wor()], [unequal_prob_wr()], [balanced_wor()]
#'
#' @examples
#' # Register a toy random sampler. Successive sampling with
#' # prob = pik does not give first-order inclusion probabilities
#' # equal to pik, so its probabilities stay "unknown" (the default).
#' my_sampler <- function(pik, n = NULL, prn = NULL, ...) {
#'   sample.int(length(pik), size = n, prob = pik)
#' }
#' register_method(
#'   "toy",
#'   type = "wor",
#'   sample_fn = my_sampler
#' )
#' s <- unequal_prob_wor(c(0.3, 0.3, 0.4), method = "toy")
#' s$method
#'
#' # Register a toy balanced sampler (ignores aux, keeps size fixed)
#' my_balanced <- function(pik, n = NULL, aux = NULL, ...) {
#'   sample.int(length(pik), size = n, prob = pik)
#' }
#' register_method("toy_bal", type = "balanced", sample_fn = my_balanced)
#' s <- balanced_wor(c(0.3, 0.3, 0.4), method = "toy_bal")
#' s$method
#'
#' # Clean up
#' unregister_method("toy")
#' unregister_method("toy_bal")
#'
#' @export
register_method <- function(
  name,
  type = c("wor", "wr", "balanced"),
  sample_fn,
  joint_fn = NULL,
  fixed_size = TRUE,
  variance_family = NULL,
  supports_prn = NULL,
  supports_aux = NULL,
  supports_strata = NULL,
  supports_spread = NULL,
  probabilities = c("unknown", "exact", "approximate")
) {
  name <- .check_method_name(name)
  type <- .match_choice(type, c("wor", "wr", "balanced"), "type")
  if (!is.function(sample_fn)) {
    stop("'sample_fn' must be a function", call. = FALSE)
  }
  if (!is.null(joint_fn) && !is.function(joint_fn)) {
    stop("'joint_fn' must be a function or NULL", call. = FALSE)
  }
  if (!isTRUE(fixed_size) && !isFALSE(fixed_size)) {
    stop("'fixed_size' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(variance_family)) {
    if (
      !is.character(variance_family) ||
        length(variance_family) != 1L ||
        !(variance_family %in% .variance_families)
    ) {
      stop(
        "'variance_family' must be NULL or one of: ",
        paste0('"', .variance_families, '"', collapse = ", "),
        call. = FALSE
      )
    }
    if (variance_family %in% c("srs", "pps_brewer") && !fixed_size) {
      stop(
        sprintf(
          "variance_family \"%s\" requires fixed_size = TRUE",
          variance_family
        ),
        call. = FALSE
      )
    }
    if (variance_family == "poisson" && (type != "wor" || fixed_size)) {
      stop(
        "variance_family \"poisson\" requires type = \"wor\" and ",
        "fixed_size = FALSE: it asserts independent selections with ",
        "a random sample size",
        call. = FALSE
      )
    }
    if (variance_family == "wr" && type != "wr") {
      stop("variance_family \"wr\" requires type = \"wr\"", call. = FALSE)
    }
    if (type == "wr" && !(variance_family %in% c("wr", "unsupported"))) {
      stop(
        sprintf(
          paste0(
            "type = \"wr\" methods only allow variance_family \"wr\" ",
            "or \"unsupported\", not \"%s\""
          ),
          variance_family
        ),
        call. = FALSE
      )
    }
    if (
      type == "balanced" &&
        !(variance_family %in% c("pps_brewer", "unsupported"))
    ) {
      stop(
        sprintf(
          paste0(
            "type = \"balanced\" methods only allow variance_family ",
            "\"pps_brewer\" or \"unsupported\", not \"%s\": balancing ",
            "constraints couple selections across units"
          ),
          variance_family
        ),
        call. = FALSE
      )
    }
  }
  probabilities <- .match_choice(
    probabilities,
    c("unknown", "exact", "approximate"),
    "probabilities"
  )
  if (type == "wr" && !fixed_size) {
    stop(
      "type = \"wr\" methods require fixed_size = TRUE",
      call. = FALSE
    )
  }
  if (type == "balanced" && !is.null(supports_prn)) {
    stop(
      "'supports_prn' only applies to type = \"wor\" or \"wr\"; ",
      "balanced_wor() has no 'prn' argument",
      call. = FALSE
    )
  }
  if (type != "balanced") {
    if (!is.null(supports_aux)) {
      stop(
        "'supports_aux' only applies to type = \"balanced\"",
        call. = FALSE
      )
    }
    if (!is.null(supports_strata)) {
      stop(
        "'supports_strata' only applies to type = \"balanced\"",
        call. = FALSE
      )
    }
    if (!is.null(supports_spread)) {
      stop(
        "'supports_spread' only applies to type = \"balanced\"",
        call. = FALSE
      )
    }
  }
  capability_flags <- list(
    supports_prn = supports_prn,
    supports_aux = supports_aux,
    supports_strata = supports_strata,
    supports_spread = supports_spread
  )
  for (flag in names(capability_flags)) {
    value <- capability_flags[[flag]]
    if (!is.null(value) && !isTRUE(value) && !isFALSE(value)) {
      stop(
        sprintf("'%s' must be TRUE or FALSE when supplied", flag),
        call. = FALSE
      )
    }
  }
  if (is.null(supports_prn)) {
    supports_prn <- FALSE
  }
  if (is.null(supports_aux)) {
    supports_aux <- type == "balanced"
  }
  if (is.null(supports_strata)) {
    supports_strata <- FALSE
  }
  if (is.null(supports_spread)) {
    supports_spread <- FALSE
  }

  # Reject collisions with built-in methods
  all_builtin <- unique(c(
    names(.wor_specs),
    names(.wr_specs),
    names(.ep_wor_specs),
    names(.ep_wr_specs),
    names(.balanced_specs)
  ))
  if (name %in% all_builtin) {
    stop(
      sprintf("'%s' is a built-in method and cannot be overridden", name),
      call. = FALSE
    )
  }
  if (is_registered_method(name)) {
    stop(
      sprintf(
        "method '%s' is already registered; call unregister_method(\"%s\") ",
        name,
        name
      ),
      "before registering a replacement",
      call. = FALSE
    )
  }

  .method_registry[[name]] <- list(
    name = name,
    type = type,
    sample_fn = sample_fn,
    joint_fn = joint_fn,
    fixed_size = fixed_size,
    variance_family = variance_family,
    supports_prn = supports_prn,
    supports_aux = supports_aux,
    supports_strata = supports_strata,
    supports_spread = supports_spread,
    probabilities = probabilities
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
  name <- .check_method_name(name)
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

#' Stop when a registered method is used through the wrong dispatcher.
#' @noRd
.check_registered_type <- function(reg, expected) {
  if (reg$type != expected) {
    entry_point <- switch(
      reg$type,
      wor = "unequal_prob_wor()",
      wr = "unequal_prob_wr()",
      balanced = "balanced_wor()"
    )
    stop(
      sprintf(
        "method '%s' is registered as type '%s'; use %s",
        reg$name,
        reg$type,
        entry_point
      ),
      call. = FALSE
    )
  }
}

#' @noRd
.check_registered_sample <- function(
  sample,
  N,
  n,
  method,
  fixed_size,
  replace
) {
  prefix <- sprintf("registered method '%s' returned", method)
  if (!is.numeric(sample) || !is.null(dim(sample))) {
    stop(prefix, " a sample that is not a numeric vector", call. = FALSE)
  }
  if (any(!is.finite(sample))) {
    stop(prefix, " non-finite sample indices", call. = FALSE)
  }
  if (any(sample != trunc(sample))) {
    stop(prefix, " non-integer sample indices", call. = FALSE)
  }
  if (any(sample < 1 | sample > N)) {
    stop(
      prefix,
      sprintf(" sample indices outside 1:%d", N),
      call. = FALSE
    )
  }
  if (!replace && anyDuplicated(sample)) {
    stop(prefix, " duplicate sample indices", call. = FALSE)
  }
  if (fixed_size && length(sample) != n) {
    stop(
      prefix,
      sprintf(" %d sample indices; expected %d", length(sample), n),
      call. = FALSE
    )
  }
  as.integer(sample)
}

#' @noRd
.dispatch_registered_wor <- function(pik, method, nrep, prn, ...) {
  reg <- .method_registry[[method]]
  .check_registered_type(reg, "wor")
  nrep <- .check_nrep_prn(
    nrep, prn, method, supports_prn = reg$supports_prn
  )

  check_pik(pik, fixed_size = reg$fixed_size)
  N <- length(pik)
  if (!is.null(prn)) {
    check_prn(prn, N)
  }
  n <- if (reg$fixed_size) as.integer(round(sum(pik))) else sum(pik)

  draw <- function() {
    .check_registered_sample(
      reg$sample_fn(pik, n = n, prn = prn, ...),
      N,
      n,
      method,
      reg$fixed_size,
      replace = FALSE
    )
  }

  if (nrep == 1L) {
    idx <- draw()
    .new_wor_sample(idx, pik, n, N, method, reg$fixed_size, "unequal_prob")
  } else if (reg$fixed_size) {
    n_int <- as.integer(round(sum(pik)))
    mat <- matrix(0L, n_int, nrep)
    for (i in seq_len(nrep)) {
      mat[, i] <- draw()
    }
    .new_wor_sample(mat, pik, n_int, N, method, TRUE, "unequal_prob")
  } else {
    samples <- lapply(seq_len(nrep), function(i) draw())
    .new_wor_sample(samples, pik, n, N, method, FALSE, "unequal_prob")
  }
}

#' @noRd
.dispatch_registered_wr <- function(hits, method, nrep, prn, ...) {
  reg <- .method_registry[[method]]
  .check_registered_type(reg, "wr")
  nrep <- .check_nrep_prn(
    nrep, prn, method, supports_prn = reg$supports_prn
  )

  check_hits(hits)
  n <- check_integer(sum(hits), "sum(hits)")
  N <- length(hits)
  if (!is.null(prn)) {
    check_prn(prn, N)
  }
  prob <- hits / sum(hits)

  draw <- function() {
    .check_registered_sample(
      reg$sample_fn(hits, n = n, prn = prn, ...),
      N,
      n,
      method,
      reg$fixed_size,
      replace = TRUE
    )
  }

  if (nrep == 1L) {
    idx <- draw()
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
      idx <- draw()
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
.dispatch_registered_balanced <- function(
  pik,
  aux,
  strata,
  spread,
  method,
  nrep,
  ...
) {
  reg <- .method_registry[[method]]
  .check_registered_type(reg, "balanced")
  nrep <- .check_nrep_prn(nrep)

  check_pik(pik, fixed_size = reg$fixed_size)
  N <- length(pik)

  if (!is.null(aux)) {
    if (!reg$supports_aux) {
      stop(
        sprintf(
          paste0(
            "method '%s' does not use auxiliary balancing variables ",
            "(registered with supports_aux = FALSE)"
          ),
          method
        ),
        call. = FALSE
      )
    }
    aux <- .check_cube_aux(aux, N)
  }

  strata_int <- NULL
  fixed_size <- reg$fixed_size
  if (!is.null(strata)) {
    if (!reg$supports_strata) {
      stop(
        sprintf(
          paste0(
            "method '%s' does not support stratified balanced sampling ",
            "(registered with supports_strata = FALSE)"
          ),
          method
        ),
        call. = FALSE
      )
    }
    strata_int <- .check_strata(strata, N)
    fixed_size <- fixed_size && .check_stratum_sizes(pik, strata_int)
  }

  if (!is.null(spread)) {
    if (!reg$supports_spread) {
      stop(
        sprintf(
          paste0(
            "method '%s' does not support spatial spreading ",
            "(registered with supports_spread = FALSE)"
          ),
          method
        ),
        call. = FALSE
      )
    }
    spread <- .check_cube_aux(spread, N, what = "spread")
    if (ncol(spread) == 0L) {
      stop("'spread' must have at least one column", call. = FALSE)
    }
  }

  n <- if (fixed_size) as.integer(round(sum(pik))) else sum(pik)

  # Optional design inputs are passed only when supplied, so minimal
  # sample_fn signatures (pik, n, aux, ...) keep working.
  sample_args <- list(pik, n = n, aux = aux)
  if (!is.null(strata_int)) {
    sample_args$strata <- strata_int
  }
  if (!is.null(spread)) {
    sample_args$spread <- spread
  }
  sample_args <- c(sample_args, list(...))

  draw <- function() {
    .check_registered_sample(
      do.call(reg$sample_fn, sample_args),
      N,
      n,
      method,
      fixed_size,
      replace = FALSE
    )
  }

  if (nrep == 1L) {
    sample_data <- draw()
  } else if (fixed_size) {
    sample_data <- matrix(0L, n, nrep)
    for (i in seq_len(nrep)) {
      sample_data[, i] <- draw()
    }
  } else {
    sample_data <- lapply(seq_len(nrep), function(i) draw())
  }

  .new_wor_sample(
    sample = sample_data,
    pik = pik,
    n = n,
    N = N,
    method = method,
    fixed_size = fixed_size,
    prob_class = "unequal_prob",
    extra_class = "balanced"
  )
}

#' @noRd
.registered_joint_or_stop <- function(
  method,
  marginals,
  sample_idx,
  quantity,
  generic_args,
  ...
) {
  if (!is_registered_method(method)) {
    .stop_no_joint(method, quantity)
  }
  reg <- .method_registry[[method]]
  if (is.null(reg$joint_fn)) {
    .stop_no_joint(method, quantity)
  }
  declared <- names(formals(reg$joint_fn))
  generic_args <- generic_args[names(generic_args) %in% declared]
  joint <- do.call(
    reg$joint_fn,
    c(list(marginals, sample_idx = sample_idx), list(...), generic_args)
  )
  expected <- if (is.null(sample_idx)) length(marginals) else length(sample_idx)
  prefix <- sprintf("registered method '%s' returned", method)
  labels <- switch(
    quantity,
    joint_inclusion_prob = list(
      result = "a joint-inclusion-probability result",
      values = "joint inclusion probabilities",
      matrix = "joint-inclusion-probability matrix"
    ),
    joint_expected_hits = list(
      result = "a joint-expected-hits result",
      values = "joint expected hits",
      matrix = "joint-expected-hits matrix"
    )
  )
  if (!is.matrix(joint) || !is.numeric(joint) || is.complex(joint)) {
    stop(
      prefix,
      " ",
      labels$result,
      " that is not a numeric matrix",
      call. = FALSE
    )
  }
  if (!identical(dim(joint), c(expected, expected))) {
    stop(
      prefix,
      sprintf(
        " a %d x %d %s; expected %d x %d",
        nrow(joint),
        ncol(joint),
        labels$matrix,
        expected,
        expected
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(joint))) {
    stop(
      prefix,
      " ",
      labels$values,
      " containing non-finite values",
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(
    joint,
    t(joint),
    tolerance = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  ))) {
    stop(prefix, " a non-symmetric ", labels$matrix, call. = FALSE)
  }
  joint
}
