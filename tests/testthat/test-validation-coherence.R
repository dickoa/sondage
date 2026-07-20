capture_error <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = identity
  )
}

test_that("choice errors identify the public argument and hide internal calls", {
  err <- capture_error(equal_prob_wor(10, 2, method = "not-a-method"))
  expect_match(conditionMessage(err), "'method' must be one of", fixed = TRUE)
  expect_null(conditionCall(err))

  err <- capture_error(
    register_method("bad-choice", type = "other", sample_fn = identity)
  )
  expect_match(conditionMessage(err), "'type' must be one of", fixed = TRUE)
  expect_null(conditionCall(err))
})

test_that("choice validation preserves unambiguous partial matching", {
  expect_identical(equal_prob_wor(10, 2, method = "sys")$method, "systematic")
  expect_identical(
    unequal_prob_wr(c(1, 1), method = "multi")$method,
    "multinomial"
  )
})

test_that("documented numeric vectors reject dimensional input", {
  expect_error(
    unequal_prob_wor(matrix(c(0.5, 0.5), ncol = 1), method = "cps"),
    "'pik' must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    unequal_prob_wr(matrix(c(1, 1), ncol = 1), method = "multinomial"),
    "'hits' must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    inclusion_prob(matrix(1:4, ncol = 2), 2),
    "'x' must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    expected_hits(matrix(1:4, ncol = 2), 2),
    "'x' must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    equal_prob_wor(4, 2, method = "bernoulli", prn = matrix(rep(0.5, 4))),
    "'prn' must be a numeric vector",
    fixed = TRUE
  )
})

test_that("documented scalar numbers reject dimensional input", {
  expect_error(equal_prob_wor(matrix(4), 2), "single number")
  expect_error(equal_prob_wor(4, matrix(2)), "single numeric value")
  expect_error(equal_prob_wor(4, 2, nrep = matrix(1)), "single number")
  expect_error(he_jip(c(0.5, 0.5), eps = matrix(1e-6)), "single numeric value")

  s <- unequal_prob_wr(c(1, 1), method = "chromy")
  expect_error(joint_expected_hits(s, nsim = matrix(10)), "single number")
})

test_that("auxiliary inputs retain useful coercions and reject non-numeric data", {
  pik <- rep(0.5, 4)
  expect_s3_class(balanced_wor(pik, aux = 1:4), "balanced")
  expect_s3_class(
    balanced_wor(pik, aux = data.frame(x = 1:4)),
    "balanced"
  )

  err <- capture_error(balanced_wor(pik, aux = list(1:4)))
  expect_match(conditionMessage(err), "numeric vector or matrix")
  expect_null(conditionCall(err))

  expect_error(
    balanced_wor(
      pik,
      bounds = list(B = matrix(1, nrow = 4), lower = matrix(0), upper = 4)
    ),
    "numeric vector of length"
  )
})

test_that("strata must be a positive integer-valued numeric vector", {
  pik <- rep(0.5, 4)
  expect_error(
    balanced_wor(pik, strata = c(1.9, 1.9, 2.1, 2.1)),
    "positive integers"
  )
  expect_error(
    balanced_wor(pik, strata = matrix(c(1, 1, 2, 2), ncol = 1)),
    "numeric vector"
  )
  expect_error(
    balanced_wor(pik, strata = factor(c(1, 1, 2, 2))),
    "numeric vector"
  )
  expect_s3_class(
    balanced_wor(pik, strata = c(10, 10, 42, 42)),
    "balanced"
  )
})

test_that("logical flags and QR tolerance are validated explicitly", {
  wor <- equal_prob_wor(4, 2)
  expect_error(
    joint_inclusion_prob(wor, sampled_only = NA),
    "'sampled_only' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    sampling_cov(wor, weighted = "yes"),
    "'weighted' must be TRUE or FALSE",
    fixed = TRUE
  )

  pik <- rep(0.5, 4)
  expect_error(
    balanced_wor(pik, aux = 1:4, condition_aux = "yes"),
    "'condition_aux' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    balanced_wor(pik, aux = 1:4, qr_tol = -1),
    "'qr_tol' must be non-negative",
    fixed = TRUE
  )
  expect_error(
    balanced_wor(pik, aux = 1:4, qr_tol = c(0, 1e-7)),
    "single numeric value"
  )
})

test_that("sample_idx is validated before joint-probability computation", {
  pik <- c(0.2, 0.3, 0.5)
  functions <- list(he_jip, hajek_jip)

  for (fn in functions) {
    expect_error(fn(pik, c(0, 2)), "between 1 and length")
    expect_error(fn(pik, c(1.9, 2)), "whole numbers")
    expect_error(fn(pik, c(NA, 2)), "missing values")
    expect_error(fn(pik, c(1, 1)), "duplicate indices")
    expect_equal(dim(fn(pik, integer())), c(0L, 0L))
  }
})

test_that("registry operations share one non-empty name contract", {
  operations <- list(method_spec, is_registered_method, unregister_method)

  for (fn in operations) {
    for (bad in list(NA_character_, "", 1, c("a", "b"))) {
      err <- capture_error(fn(bad))
      expect_match(conditionMessage(err), "non-empty character string")
      expect_null(conditionCall(err))
    }
  }
})
