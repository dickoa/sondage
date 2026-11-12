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
  expect_error(equal_prob_wor(matrix(4), 2), "single numeric value")
  expect_error(equal_prob_wor(4, matrix(2)), "single numeric value")
  expect_error(equal_prob_wor(4, 2, nrep = matrix(1)), "single numeric value")
  expect_error(he_jip(c(0.5, 0.5), eps = matrix(1e-6)), "single numeric value")

  s <- unequal_prob_wr(c(1, 1), method = "chromy")
  expect_error(joint_expected_hits(s, nsim = matrix(10)), "single numeric value")
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

## Sizes beyond the integer range

test_that("a sample size above the integer range is refused before C", {
  big <- .Machine$integer.max + 1
  for (m in c("chromy", "multinomial", "systematic_mr")) {
    expect_no_warning(expect_error(
      unequal_prob_wr(c(big, 0), method = m),
      "exceeds the largest supported sample size"
    ))
  }
  for (m in c("srs", "systematic", "bernoulli")) {
    expect_no_warning(expect_error(
      equal_prob_wor(big, 2, method = m),
      "'N' = 2.14748e\\+09 exceeds the largest supported integer"
    ))
  }
  expect_no_warning(expect_error(
    equal_prob_wr(10, big),
    "exceeds the largest supported integer"
  ))
})

## Edited design objects before a sampled-only query

test_that("sampled-only queries refuse units the compiled code cannot index", {
  # The submatrix code indexes the population by $sample. An index out of
  # range read outside pik, and a duplicate left a CPS pair uncomputed.
  p <- c(0.3, 0.3 + 1e-9, 0.4, 0.5, 0.5)
  p <- p * 2 / sum(p)
  for (m in c("cps", "sampford", "systematic", "brewer", "poisson")) {
    s <- unequal_prob_wor(p, method = m)
    for (bad in list(c(1L, 1L), c(0L, 2L), c(2L, 7L), c(1, NA))) {
      s$sample <- bad
      expect_error(
        joint_inclusion_prob(s, sampled_only = TRUE),
        "distinct unit indices between 1 and N = 5",
        info = paste(m, deparse(bad))
      )
    }
  }
  s <- unequal_prob_wor(p, method = "cps")
  s$sample <- integer(0)
  expect_identical(dim(joint_inclusion_prob(s, sampled_only = TRUE)), c(0L, 0L))

  w <- unequal_prob_wr(c(0.5, 1.5, 1), method = "systematic_mr")
  w$hits <- c(w$hits, 1L, 1L)
  expect_error(
    joint_expected_hits(w, sampled_only = TRUE),
    "one count per unit, N = 3"
  )
})

test_that("the systematic_mr kernel refuses expected hits it cannot split", {
  # The R entry points refuse these first. Called directly, the kernel
  # used to cast them to int, which is undefined.
  expect_error(
    .Call(C_up_wr_systematic, c(2147483648, 0), 1L),
    "do not add up"
  )
  expect_error(.Call(C_up_wr_systematic, c(NaN, 1), 1L), "do not add up")
  expect_error(.Call(C_wr_systematic_jeh, c(Inf, 1), 1L), "do not add up")
})
