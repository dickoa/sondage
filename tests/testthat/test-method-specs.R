# Regression tests for method spec tables, resolver, and hook helpers.
# These lock in the internal contracts introduced by the registration-readiness
# refactor and will catch drift if the spec tables are modified incorrectly.

# --- Spec table correctness ---

test_that("WOR spec table has correct PRN capabilities", {
  expect_true(.method_supports_prn("sps", "wor"))
  expect_true(.method_supports_prn("pareto", "wor"))
  expect_true(.method_supports_prn("poisson", "wor"))
  expect_false(.method_supports_prn("cps", "wor"))
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
  expect_false("systematic" %in% .he_jip_methods)
  expect_false("poisson" %in% .he_jip_methods)
})

# --- Resolver ---

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

# --- Hook helper error messages (byte-compatible with originals) ---

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

# --- method_spec() ---

test_that("method_spec returns correct metadata for WOR methods", {
  spec <- method_spec("brewer")
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
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)

  spec <- method_spec("multinomial")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("method_spec returns correct metadata for EP methods", {
  spec <- method_spec("srs")
  expect_equal(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)

  spec <- method_spec("bernoulli")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)
})

test_that("method_spec returns correct metadata for balanced methods", {
  spec <- method_spec("cube")
  expect_equal(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("method_spec returns NULL for unknown methods", {
  expect_null(method_spec("nonexistent"))
  expect_null(method_spec("lpm"))
})

test_that("method_spec validates input", {
  expect_error(method_spec(42), "single character string")
  expect_error(method_spec(c("a", "b")), "single character string")
})

test_that("method_spec returns correct metadata for registered methods", {
  on.exit(unregister_method("test_ms"), add = TRUE)
  register_method("test_ms", "wor", sample_fn = identity,
                   fixed_size = FALSE, supports_prn = TRUE)

  spec <- method_spec("test_ms")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)
})

test_that("method_spec picks up registered WR methods", {
  on.exit(unregister_method("test_wr"), add = TRUE)
  register_method("test_wr", "wr", sample_fn = identity)

  spec <- method_spec("test_wr")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

# --- End-to-end: dispatchers reject unknown methods ---

test_that("dispatchers reject unknown methods", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "lpm"))
  expect_error(unequal_prob_wr(c(1, 1), method = "lpm"))
  expect_error(equal_prob_wor(N = 10, n = 2, method = "lpm"))
})

# --- End-to-end: .stop_no_joint fires via fake design objects ---

test_that("joint_inclusion_prob errors for unknown WOR method", {
  fake <- structure(
    list(
      sample = 1L, pik = c(0.5, 0.5), n = 1L, N = 2L,
      method = "fake", fixed_size = TRUE
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
      sample = 1L, prob = c(0.5, 0.5), hits = c(1L, 0L),
      n = 1L, N = 2L, method = "fake", fixed_size = TRUE
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
