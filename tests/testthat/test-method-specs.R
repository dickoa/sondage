# Regression tests for method spec tables, resolver, and hook helpers.
# These lock in the internal contracts introduced by the registration-readiness
# refactor and will catch drift if the spec tables are modified incorrectly.

# Spec table correctness
test_that("WOR spec table has correct PRN capabilities", {
  expect_true(.method_supports_prn("sps", "wor"))
  expect_true(.method_supports_prn("pareto", "wor"))
  expect_true(.method_supports_prn("poisson", "wor"))
  expect_false(.method_supports_prn("cps", "wor"))
  expect_false(.method_supports_prn("sampford", "wor"))
  expect_false(.method_supports_prn("brewer", "wor"))
  expect_false(.method_supports_prn("systematic", "wor"))
})

test_that("WR spec table has correct PRN capabilities", {
  expect_false(.method_supports_prn("chromy", "wr"))
  expect_false(.method_supports_prn("multinomial", "wr"))
})

test_that("EP WOR spec table has correct PRN capabilities", {
  expect_true(.method_supports_prn("bernoulli", "ep_wor"))
  expect_false(.method_supports_prn("srs", "ep_wor"))
  expect_false(.method_supports_prn("systematic", "ep_wor"))
})

test_that("fixed_size is correct for WOR methods", {
  expect_true(.method_is_fixed_size("cps", "wor"))
  expect_true(.method_is_fixed_size("sampford", "wor"))
  expect_true(.method_is_fixed_size("brewer", "wor"))
  expect_true(.method_is_fixed_size("systematic", "wor"))
  expect_true(.method_is_fixed_size("sps", "wor"))
  expect_true(.method_is_fixed_size("pareto", "wor"))
  expect_false(.method_is_fixed_size("poisson", "wor"))
})

test_that("fixed_size is correct for WR methods", {
  expect_true(.method_is_fixed_size("chromy", "wr"))
  expect_true(.method_is_fixed_size("multinomial", "wr"))
})

test_that("fixed_size is correct for EP WOR methods", {
  expect_true(.method_is_fixed_size("srs", "ep_wor"))
  expect_true(.method_is_fixed_size("systematic", "ep_wor"))
  expect_false(.method_is_fixed_size("bernoulli", "ep_wor"))
})

test_that("HE JIP methods list is complete", {
  expect_true(all(c("brewer", "sps", "pareto", "cube") %in% .he_jip_methods))
  expect_false("cps" %in% .he_jip_methods)
  expect_false("sampford" %in% .he_jip_methods)
  expect_false("systematic" %in% .he_jip_methods)
  expect_false("poisson" %in% .he_jip_methods)
})

test_that(".get_builtin_spec returns spec for known methods", {
  spec <- .get_builtin_spec("cps", "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$prn)

  spec <- .get_builtin_spec("poisson", "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$prn)

  spec <- .get_builtin_spec("bernoulli", "ep_wor")
  expect_false(spec$fixed_size)
  expect_true(spec$prn)
})

test_that(".get_builtin_spec returns NULL for unknown methods", {
  expect_null(.get_builtin_spec("lpm", "wor"))
  expect_null(.get_builtin_spec("bernoulli", "wor"))
  expect_null(.get_builtin_spec("cps", "ep_wor"))
})

test_that(".get_builtin_spec rejects invalid context", {
  expect_error(.get_builtin_spec("cps", "nonsense"), "unknown context")
})

# Hook helper error messages (byte-compatible with originals)
test_that(".stop_no_joint preserves exact generics.R error format", {
  expect_error(
    .stop_no_joint("lpm", "joint_inclusion_prob"),
    "^joint_inclusion_prob not implemented for method 'lpm'$"
  )
  expect_error(
    .stop_no_joint("lpm", "joint_expected_hits"),
    "^joint_expected_hits not implemented for method 'lpm'$"
  )
})

test_that(".stop_unknown_method produces informative error", {
  expect_error(
    .stop_unknown_method("lpm"),
    "^unknown sampling method 'lpm'$"
  )
})

test_that("method_spec returns correct metadata for WOR methods", {
  spec <- method_spec("brewer")
  expect_equal(spec$dispatcher, "unequal_prob_wor")
  expect_equal(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)

  spec <- method_spec("poisson")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)

  spec <- method_spec("sps")
  expect_equal(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_true(spec$supports_prn)
})

test_that("method_spec returns correct metadata for WR methods", {
  spec <- method_spec("chromy")
  expect_equal(spec$dispatcher, "unequal_prob_wr")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)

  spec <- method_spec("multinomial")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("method_spec returns correct metadata for EP methods", {
  wor <- method_spec("srs", dispatcher = "equal_prob_wor")
  expect_equal(wor$dispatcher, "equal_prob_wor")
  expect_equal(wor$type, "wor")
  expect_true(wor$fixed_size)
  expect_equal(wor$variance_family, "srs")

  wr <- method_spec("srs", dispatcher = "equal_prob_wr")
  expect_equal(wr$dispatcher, "equal_prob_wr")
  expect_equal(wr$type, "wr")
  expect_true(wr$fixed_size)
  expect_equal(wr$variance_family, "wr")

  spec <- method_spec("bernoulli")
  expect_equal(spec$dispatcher, "equal_prob_wor")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)
})

test_that("method_spec returns correct metadata for balanced methods", {
  spec <- method_spec("cube")
  expect_equal(spec$dispatcher, "balanced_wor")
  expect_equal(spec$type, "balanced")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
  expect_true(spec$supports_aux)
  expect_true(spec$supports_strata)
  expect_false(spec$supports_spread)
})

test_that("method_spec returns NULL for unknown methods", {
  expect_null(method_spec("nonexistent"))
  expect_null(method_spec("lpm"))
  expect_null(method_spec("nonexistent", dispatcher = "equal_prob_wor"))
})

test_that("method_spec validates input", {
  expect_error(method_spec(42), "non-empty character string")
  expect_error(method_spec(c("a", "b")), "non-empty character string")
  expect_error(method_spec("brewer", dispatcher = "unequal"), "exactly one of")
  expect_error(
    method_spec("brewer", dispatcher = c("unequal_prob_wor", "equal_prob_wor")),
    "exactly one of"
  )
  expect_error(
    method_spec("brewer", dispatcher = NA_character_),
    "exactly one of"
  )
  expect_error(
    method_spec("brewer", dispatcher = matrix("unequal_prob_wor")),
    "exactly one of"
  )
})

test_that("method_spec requires a dispatcher for ambiguous names", {
  expect_error(
    method_spec("srs"),
    'ambiguous.*"equal_prob_wor" or "equal_prob_wr"'
  )
  expect_error(
    method_spec("systematic"),
    'ambiguous.*"equal_prob_wor" or "unequal_prob_wor"'
  )
})

test_that("method_spec rejects incompatible dispatchers", {
  expect_error(
    method_spec("brewer", dispatcher = "equal_prob_wor"),
    'not available.*use "unequal_prob_wor"'
  )
  expect_error(
    method_spec("srs", dispatcher = "balanced_wor"),
    'not available.*"equal_prob_wor" or "equal_prob_wr"'
  )
})

test_that("method_spec distinguishes systematic variants", {
  equal <- method_spec("systematic", dispatcher = "equal_prob_wor")
  unequal <- method_spec("systematic", dispatcher = "unequal_prob_wor")

  expect_equal(equal$type, "wor")
  expect_equal(unequal$type, "wor")
  expect_equal(equal$variance_family, "srs")
  expect_equal(unequal$variance_family, "pps_brewer")
})

test_that("method_spec returns correct metadata for registered methods", {
  on.exit(unregister_method("test_ms"), add = TRUE)
  register_method(
    "test_ms",
    "wor",
    sample_fn = identity,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

  spec <- method_spec("test_ms")
  expect_equal(spec$dispatcher, "unequal_prob_wor")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)

  scoped <- method_spec("test_ms", dispatcher = "unequal_prob_wor")
  expect_identical(scoped, spec)
  expect_error(
    method_spec("test_ms", dispatcher = "balanced_wor"),
    'uses dispatcher "unequal_prob_wor", not "balanced_wor"'
  )
})

test_that("method_spec picks up registered WR methods", {
  on.exit(unregister_method("test_wr"), add = TRUE)
  register_method("test_wr", "wr", sample_fn = identity)

  spec <- method_spec("test_wr")
  expect_equal(spec$dispatcher, "unequal_prob_wr")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("dispatchers reject unknown methods", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "lpm"))
  expect_error(unequal_prob_wr(c(1, 1), method = "lpm"))
  expect_error(equal_prob_wor(N = 10, n = 2, method = "lpm"))
})

test_that("joint_inclusion_prob errors for unknown WOR method", {
  fake <- structure(
    list(
      sample = 1L,
      pik = c(0.5, 0.5),
      n = 1L,
      N = 2L,
      method = "fake",
      fixed_size = TRUE
    ),
    class = c("unequal_prob", "wor", "sondage_sample")
  )
  expect_error(
    joint_inclusion_prob(fake),
    "^joint_inclusion_prob not implemented for method 'fake'$"
  )
  expect_error(
    joint_inclusion_prob(fake, sampled_only = TRUE),
    "^joint_inclusion_prob not implemented for method 'fake'$"
  )
})

test_that("joint_expected_hits errors for unknown WR method", {
  fake <- structure(
    list(
      sample = 1L,
      prob = c(0.5, 0.5),
      hits = c(1L, 0L),
      n = 1L,
      N = 2L,
      method = "fake",
      fixed_size = TRUE
    ),
    class = c("unequal_prob", "wr", "sondage_sample")
  )
  expect_error(
    joint_expected_hits(fake),
    "^joint_expected_hits not implemented for method 'fake'$"
  )
  expect_error(
    joint_expected_hits(fake, sampled_only = TRUE),
    "^joint_expected_hits not implemented for method 'fake'$"
  )
})

test_that(".get_builtin_spec resolves balanced context", {
  spec <- sondage:::.get_builtin_spec("cube", "balanced")
  expect_type(spec, "list")
  expect_true(isTRUE(spec$fixed_size))
  expect_false(isTRUE(spec$prn))
  # Non-cube method in balanced context returns NULL
  expect_null(sondage:::.get_builtin_spec("srs", "balanced"))
})

# variance_family in the spec tables
test_that("built-in spec tables declare correct variance families", {
  # Fixed-size unequal-probability WOR: Brewer treatment
  for (m in c("cps", "sampford", "brewer", "systematic", "sps", "pareto")) {
    expect_identical(
      .get_builtin_spec(m, "wor")$variance_family,
      "pps_brewer",
      info = m
    )
  }
  # Random-size with independent selections: Poisson treatment
  expect_identical(.get_builtin_spec("poisson", "wor")$variance_family, "poisson")
  expect_identical(
    .get_builtin_spec("bernoulli", "ep_wor")$variance_family,
    "poisson"
  )
  # Equal-probability fixed-size WOR: SRS treatment
  expect_identical(.get_builtin_spec("srs", "ep_wor")$variance_family, "srs")
  expect_identical(
    .get_builtin_spec("systematic", "ep_wor")$variance_family,
    "srs"
  )
  # With-replacement / PMR: Hansen-Hurwitz treatment
  expect_identical(.get_builtin_spec("chromy", "wr")$variance_family, "wr")
  expect_identical(.get_builtin_spec("multinomial", "wr")$variance_family, "wr")
  expect_identical(.get_builtin_spec("srs", "ep_wr")$variance_family, "wr")
  # Balanced: Brewer treatment
  expect_identical(
    .get_builtin_spec("cube", "balanced")$variance_family,
    "pps_brewer"
  )
})

test_that("all built-in variance families are in the allowed set", {
  all_specs <- c(
    sondage:::.wor_specs,
    sondage:::.wr_specs,
    sondage:::.ep_wor_specs,
    sondage:::.ep_wr_specs,
    sondage:::.balanced_specs
  )
  for (m in names(all_specs)) {
    expect_true(
      all_specs[[m]]$variance_family %in% sondage:::.variance_families,
      info = m
    )
  }
})

test_that("method_spec reports variance_family for built-ins", {
  expect_identical(method_spec("brewer")$variance_family, "pps_brewer")
  expect_identical(method_spec("poisson")$variance_family, "poisson")
  expect_identical(method_spec("bernoulli")$variance_family, "poisson")
  expect_identical(method_spec("chromy")$variance_family, "wr")
  expect_identical(method_spec("multinomial")$variance_family, "wr")
  expect_identical(method_spec("cube")$variance_family, "pps_brewer")
  expect_identical(
    method_spec("systematic", "equal_prob_wor")$variance_family,
    "srs"
  )
  expect_identical(
    method_spec("systematic", "unequal_prob_wor")$variance_family,
    "pps_brewer"
  )
  expect_identical(method_spec("srs", "equal_prob_wor")$variance_family, "srs")
  expect_identical(method_spec("srs", "equal_prob_wr")$variance_family, "wr")
})
