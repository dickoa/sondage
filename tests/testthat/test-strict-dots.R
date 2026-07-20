test_that("built-in samplers reject unused arguments in dots", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  hits <- c(0.4, 0.8, 0.5, 0.6, 0.7)
  spread <- cbind(seq_along(pik), rev(seq_along(pik)))

  expect_error(equal_prob_wor(10, 3, nrepp = 2), "nrepp")
  expect_error(equal_prob_wor(10, 3, nrep = 2, typo = TRUE), "typo")
  expect_error(equal_prob_wr(10, 3, typo = TRUE), "typo")
  expect_error(equal_prob_wr(10, 3, nrep = 2, typo = TRUE), "typo")

  expect_error(unequal_prob_wor(pik, typo = TRUE), "typo")
  expect_error(
    unequal_prob_wor(pik, method = "cps", nrep = 2, typo = TRUE),
    "typo"
  )
  expect_error(unequal_prob_wr(hits, typo = TRUE), "typo")
  expect_error(
    unequal_prob_wr(hits, method = "multinomial", nrep = 2, typo = TRUE),
    "typo"
  )

  expect_error(balanced_wor(pik, typo = TRUE), "'typo'")
  expect_error(
    balanced_wor(pik, spread = spread, method = "lpm2", qr_tol = 1e-8),
    "'qr_tol'"
  )
  expect_error(
    equal_prob_wor(10, 3, method = "srs", nrep = 1, prn = NULL, TRUE),
    "unused argument"
  )
})

test_that("strict dots checks do not force rejected arguments", {
  forced <- FALSE
  expect_error(
    equal_prob_wor(10, 3, typo = {
      forced <<- TRUE
      1
    }),
    "typo"
  )
  expect_false(forced)
})

test_that("documented balanced-method dots remain supported", {
  pik <- rep(0.5, 8)
  aux <- cbind(seq_along(pik), 2 * seq_along(pik))
  spread <- cbind(seq_along(pik), rev(seq_along(pik)))

  cube <- balanced_wor(
    pik,
    aux = aux,
    eps = 1e-9,
    condition_aux = TRUE,
    qr_tol = 1e-8
  )
  lpm2 <- balanced_wor(pik, spread = spread, method = "lpm2", eps = 1e-9)
  scps <- balanced_wor(pik, spread = spread, method = "scps", eps = 1e-9)

  expect_s3_class(cube, "sondage_sample")
  expect_s3_class(lpm2, "sondage_sample")
  expect_s3_class(scps, "sondage_sample")
})

test_that("removed unequal-probability eps keeps its tailored error", {
  pik <- c(0.2, 0.4, 0.6, 0.8)

  expect_error(unequal_prob_wor(pik, eps = 1e-6), "no longer modifies")
  expect_error(
    unequal_prob_wor(pik, method = "poisson", eps = 1e-6),
    "no longer modifies"
  )
  expect_s3_class(unequal_prob_wor(pik, eps = NULL), "sondage_sample")
})

test_that("built-in design queries reject unused arguments", {
  wor <- equal_prob_wor(10, 3)
  wr <- equal_prob_wr(10, 3)

  expect_error(inclusion_prob(1:5, n = 2, typo = TRUE), "typo")
  expect_error(inclusion_prob(wor, n = 2), "'n'")
  expect_error(expected_hits(1:5, n = 2, typo = TRUE), "typo")
  expect_error(expected_hits(wr, n = 2), "'n'")
  expect_error(joint_inclusion_prob(wor, typo = TRUE), "'typo'")
  expect_error(joint_expected_hits(wr, typo = TRUE), "'typo'")
  expect_error(sampling_cov(wor, typo = TRUE), "'typo'")
  expect_error(sampling_cov(wr, typo = TRUE), "'typo'")
})

test_that("registered methods continue to own their dots", {
  on.exit(unregister_method("dots_wor"), add = TRUE)
  sample_token <- NULL
  joint_token <- NULL

  sampler <- function(pik, n = NULL, prn = NULL, token = NULL, ...) {
    sample_token <<- token
    seq_len(n)
  }
  joint <- function(pik, sample_idx = NULL, token = NULL, ...) {
    joint_token <<- token
    if (!is.null(sample_idx)) {
      pik <- pik[sample_idx]
    }
    out <- outer(pik, pik)
    diag(out) <- pik
    out
  }
  register_method(
    "dots_wor",
    "wor",
    sample_fn = sampler,
    joint_fn = joint
  )

  design <- unequal_prob_wor(
    rep(0.5, 4),
    method = "dots_wor",
    token = "sample"
  )
  joint_inclusion_prob(design, token = "joint")

  expect_identical(sample_token, "sample")
  expect_identical(joint_token, "joint")
})
