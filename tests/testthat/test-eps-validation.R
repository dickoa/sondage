test_that("fixed-size WOR methods validate eps in (0, 0.5)", {
  pik <- c(0.2, 0.3, 0.5)
  for (method in c("cps", "brewer", "systematic", "sps", "pareto")) {
    expect_error(
      unequal_prob_wor(pik, method = method, eps = 0),
      "open interval",
      info = method
    )
    expect_error(
      unequal_prob_wor(pik, method = method, eps = 0.5),
      "open interval",
      info = method
    )
  }
})

test_that("batch cps validates eps and fixed-size pik constraints", {
  pik <- c(0.2, 0.3, 0.5)
  expect_error(
    unequal_prob_wor(pik, method = "cps", nrep = 2, eps = Inf),
    "finite"
  )

  expect_error(
    unequal_prob_wor(c(0.49, 0.49, 0.49), method = "cps", nrep = 2),
    "not close to an integer"
  )
})

test_that("balanced_wor validates eps in single and batch mode", {
  pik <- c(0.2, 0.3, 0.5)
  expect_error(balanced_wor(pik, eps = 0), "open interval")
  expect_error(balanced_wor(pik, nrep = 2, eps = 0.5), "open interval")
})

test_that("joint_inclusion_prob validates eps", {
  s <- unequal_prob_wor(c(0.2, 0.3, 0.5), method = "cps")
  expect_error(joint_inclusion_prob(s, eps = 0), "open interval")
  expect_error(joint_inclusion_prob(s, eps = NA_real_), "must not be NA")
})

test_that("bernoulli accepts non-integer expected sample size", {
  expect_no_error(equal_prob_wor(10, 2.5, method = "bernoulli"))
})
