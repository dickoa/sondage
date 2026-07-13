test_that("fixed-size WOR methods reject the removed eps argument", {
  pik <- c(0.2, 0.3, 0.5)
  for (method in c("cps", "sampford", "brewer", "systematic", "sps", "pareto")) {
    expect_error(
      unequal_prob_wor(pik, method = method, eps = 1e-6),
      "no longer modifies the design",
      info = method
    )
  }
  expect_error(
    unequal_prob_wor(pik, method = "cps", nrep = 2, eps = 1e-6),
    "no longer modifies the design"
  )
})

test_that("batch cps validates fixed-size pik constraints", {
  expect_error(
    unequal_prob_wor(c(0.49, 0.49, 0.49), method = "cps", nrep = 2),
    "not close to an integer"
  )
})

test_that("balanced_wor validates eps in single and batch mode", {
  pik <- c(0.2, 0.3, 0.5)
  expect_error(balanced_wor(pik, eps = 0), "open interval")
  expect_error(balanced_wor(pik, nrep = 2, eps = 0.5), "open interval")
  expect_error(balanced_wor(pik, eps = "x"), "single numeric")
  expect_error(balanced_wor(pik, eps = NA_real_), "must not be NA")
})

test_that("balanced_wor eps does not reclassify input pik", {
  # interior values within eps of the boundary are rejected, not snapped
  # (sums stay integer to fp accuracy so the size check passes first)
  expect_error(
    balanced_wor(c(1e-8, 0.5, 0.5 - 1e-8), eps = 1e-6),
    "exactly 0 or exactly 1"
  )
  expect_error(
    balanced_wor(c(1 - 1e-8, 0.5, 0.5 + 1e-8), eps = 1e-6, nrep = 2),
    "exactly 0 or exactly 1"
  )
})

test_that("joint_inclusion_prob validates eps", {
  s <- unequal_prob_wor(c(0.2, 0.3, 0.5), method = "brewer")
  expect_error(joint_inclusion_prob(s, eps = 0), "open interval")
  expect_error(joint_inclusion_prob(s, eps = NA_real_), "must not be NA")
})

test_that("bernoulli accepts non-integer expected sample size", {
  expect_no_error(equal_prob_wor(10, 2.5, method = "bernoulli"))
})
