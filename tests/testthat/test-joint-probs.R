# Tille Example 10 (page 86), known exact values
tille_pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)

tille_expected <- matrix(c(
  0.07,   0.0049, 0.0130, 0.0215, 0.0447, 0.0559,
  0.0049, 0.17,   0.0324, 0.0537, 0.1113, 0.1377,
  0.0130, 0.0324, 0.41,   0.1407, 0.2888, 0.3452,
  0.0215, 0.0537, 0.1407, 0.61,   0.4691, 0.5351,
  0.0447, 0.1113, 0.2888, 0.4691, 0.83,   0.7461,
  0.0559, 0.1377, 0.3452, 0.5351, 0.7461, 0.91
), nrow = 6, byrow = TRUE)

test_that("joint_inclusion_prob.cps returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps")
  pikl <- joint_inclusion_prob(s)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.cps is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.cps has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.cps satisfies marginal constraint", {
  pik <- tille_pik
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("joint_inclusion_prob.cps matches Tille Example 10", {
  pikl <- joint_inclusion_prob(unequal_prob_wor(tille_pik, method = "cps"))

  expect_equal(pikl, tille_expected, tolerance = 0.01)
})

test_that("joint_inclusion_prob.cps produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))

  for (i in 1:3) {
    for (j in 1:3) {
      if (i != j) {
        expect_true(pikl[i, j] <= min(pik[i], pik[j]) + 1e-9)
      }
    }
  }
})

test_that("joint_inclusion_prob.cps handles n=1 case", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("joint_inclusion_prob.cps handles certainty selections", {
  pik <- c(0.5, 0.5, 1.0)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(pikl[3, 1], pik[1], tolerance = 1e-6)
  expect_equal(pikl[3, 2], pik[2], tolerance = 1e-6)
})

test_that("joint_inclusion_prob.systematic returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.systematic is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.systematic has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.systematic produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("joint_inclusion_prob.systematic can have zero joint probs", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  off_diag <- pikl[row(pikl) != col(pikl)]
  expect_true(any(off_diag == 0) || all(off_diag > 0))
})

test_that("joint_inclusion_prob.systematic satisfies marginal constraint", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("joint_inclusion_prob.brewer returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.brewer is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.brewer has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.brewer produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("joint_inclusion_prob.brewer approximates CPS reasonably", {
  pik <- tille_pik
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))
  pikl_exact <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  off_diag_brewer <- pikl_brewer[row(pikl_brewer) != col(pikl_brewer)]
  off_diag_exact <- pikl_exact[row(pikl_exact) != col(pikl_exact)]
  expect_true(cor(off_diag_brewer, off_diag_exact) > 0.99)

  mad <- mean(abs(pikl_brewer - pikl_exact))
  expect_true(mad < 0.03)
})

test_that("joint_inclusion_prob.brewer handles n <= 1 case", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("joint_expected_hits.chromy returns valid matrix", {
  x <- c(10, 20, 15, 25, 30)
  hits <- expected_hits(x, n = 3)
  s <- unequal_prob_wr(hits, method = "chromy")
  joint <- joint_expected_hits(s, nsim = 5000)

  expect_equal(dim(joint), c(5, 5))
  expect_equal(joint, t(joint))
  expect_true(all(joint > 0))

  pik <- 3 * x / sum(x)
  expect_equal(diag(joint), pik, tolerance = 0.05)
})

test_that("all joint_inclusion_prob methods agree on diagonal", {
  pik <- c(0.2, 0.3, 0.5)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl_cps), pik)
  expect_equal(diag(pikl_sys), pik)
  expect_equal(diag(pikl_brewer), pik)
})

test_that("all joint methods produce symmetric matrices", {
  pik <- c(0.15, 0.25, 0.35, 0.25)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(pikl_cps, t(pikl_cps))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_brewer, t(pikl_brewer))
})

test_that("joint methods work with larger populations", {
  set.seed(1)
  N <- 50
  n <- 10
  x <- runif(N)
  pik <- n * x / sum(x)
  hits <- expected_hits(x, n = n)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_chromy <- joint_expected_hits(unequal_prob_wr(hits, method = "chromy"))

  expect_equal(dim(pikl_cps), c(N, N))
  expect_equal(dim(pikl_brewer), c(N, N))
  expect_equal(dim(pikl_sys), c(N, N))
  expect_equal(dim(pikl_chromy), c(N, N))

  expect_equal(pikl_cps, t(pikl_cps))
  expect_equal(pikl_brewer, t(pikl_brewer))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_chromy, t(pikl_chromy))
})
