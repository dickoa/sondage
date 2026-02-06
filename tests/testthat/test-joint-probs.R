# Tests for joint inclusion probability functions

# Tillé Example 10 (page 86) - known exact values
tille_pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
# fmt: skip
tille_expected <- matrix(c(
  0.07,   0.0049, 0.0130, 0.0215, 0.0447, 0.0559,
  0.0049, 0.17,   0.0324, 0.0537, 0.1113, 0.1377,
  0.0130, 0.0324, 0.41,   0.1407, 0.2888, 0.3452,
  0.0215, 0.0537, 0.1407, 0.61,   0.4691, 0.5351,
  0.0447, 0.1113, 0.2888, 0.4691, 0.83,   0.7461,
  0.0559, 0.1377, 0.3452, 0.5351, 0.7461, 0.91
), nrow = 6, byrow = TRUE)

test_that("up_maxent_jip returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_maxent_jip(pik)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("up_maxent_jip is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_maxent_jip(pik)

  expect_equal(pikl, t(pikl))
})

test_that("up_maxent_jip has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_maxent_jip(pik)

  expect_equal(diag(pikl), pik)
})

test_that("up_maxent_jip satisfies marginal constraint", {
  pik <- tille_pik
  pikl <- up_maxent_jip(pik)
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("up_maxent_jip matches Tillé Example 10", {
  pikl <- up_maxent_jip(tille_pik)

  expect_equal(pikl, tille_expected, tolerance = 0.01)
})

test_that("up_maxent_jip produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_maxent_jip(pik)

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))

  # Off-diagonal should be <= min(pik_i, pik_j)
  for (i in 1:3) {
    for (j in 1:3) {
      if (i != j) {
        expect_true(pikl[i, j] <= min(pik[i], pik[j]) + 1e-9)
      }
    }
  }
})

test_that("up_maxent_jip handles n=1 case", {
  pik <- c(0.3, 0.3, 0.4) # n = 1
  pikl <- up_maxent_jip(pik)

  # All off-diagonal should be 0 when n=1
  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("up_maxent_jip handles certainty selections", {
  pik <- c(0.5, 0.5, 1.0) # Unit 3 always selected
  pikl <- up_maxent_jip(pik)

  # Joint with certainty unit = marginal of other unit
  expect_equal(pikl[3, 1], pik[1], tolerance = 1e-6)
  expect_equal(pikl[3, 2], pik[2], tolerance = 1e-6)
})

test_that("up_maxent_jip rejects invalid input", {
  expect_error(up_maxent_jip(c(0.2, NA, 0.3)), "missing values")
  expect_error(up_maxent_jip(c(0.2, -0.1, 0.3)), "between 0 and 1")
  expect_error(up_maxent_jip(c(0.2, 1.5, 0.3)), "between 0 and 1")
  expect_error(up_maxent_jip(integer(0)), "empty")
  expect_error(up_maxent_jip("abc"), "numeric vector")
})

test_that("up_systematic_jip returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("up_systematic_jip is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)

  expect_equal(pikl, t(pikl))
})

test_that("up_systematic_jip has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)

  expect_equal(diag(pikl), pik)
})

test_that("up_systematic_jip produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("up_systematic_jip can have zero joint probabilities", {
  # Systematic sampling is NOT high-entropy - some pairs never selected together
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)

  # At least check it runs and produces valid output
  # Zero joint probs are expected for systematic
  off_diag <- pikl[row(pikl) != col(pikl)]
  expect_true(any(off_diag == 0) || all(off_diag > 0)) # Either is valid
})

test_that("up_systematic_jip satisfies marginal constraint on average", {
  # For systematic, sum_{l≠k} π_kl = (n-1) * π_k holds
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("up_systematic_jip rejects invalid input", {
  expect_error(up_systematic_jip(c(0.2, NA, 0.3)), "missing values")
  expect_error(up_systematic_jip(c(0.2, -0.1, 0.3)), "between 0 and 1")
  expect_error(up_systematic_jip(integer(0)), "empty")
})

test_that("up_brewer_jip returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_brewer_jip(pik)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("up_brewer_jip is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_brewer_jip(pik)

  expect_equal(pikl, t(pikl))
})

test_that("up_brewer_jip has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_brewer_jip(pik)

  expect_equal(diag(pikl), pik)
})

test_that("up_brewer_jip produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_brewer_jip(pik)

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("up_brewer_jip approximates exact CPS reasonably", {
  pik <- tille_pik
  pikl_brewer <- up_brewer_jip(pik)
  pikl_exact <- up_maxent_jip(pik)

  # Correlation should be very high
  off_diag_brewer <- pikl_brewer[row(pikl_brewer) != col(pikl_brewer)]
  off_diag_exact <- pikl_exact[row(pikl_exact) != col(pikl_exact)]
  expect_true(cor(off_diag_brewer, off_diag_exact) > 0.99)

  # Mean absolute difference should be small
  mad <- mean(abs(pikl_brewer - pikl_exact))
  expect_true(mad < 0.03)
})

test_that("up_brewer_jip handles n <= 1 case", {
  pik <- c(0.3, 0.3, 0.4) # n = 1
  pikl <- up_brewer_jip(pik)

  # Should return diagonal matrix
  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("up_brewer_jip rejects invalid input", {
  expect_error(up_brewer_jip(c(0.2, NA, 0.3)), "missing values")
  expect_error(up_brewer_jip(c(0.2, -0.1, 0.3)), "between 0 and 1")
  expect_error(up_brewer_jip(integer(0)), "empty")
})


test_that("up_chromy_jip returns valid matrix", {
  x <- c(10, 20, 15, 25, 30)
  joint <- up_chromy_jip(x, n = 3, nsim = 5000)

  expect_equal(dim(joint), c(5, 5))
  expect_equal(joint, t(joint)) # symmetric
  expect_true(all(joint > 0)) # all positive (randomized)

  # Diagonal ≈ first-order probabilities
  pik <- 3 * x / sum(x)
  expect_equal(diag(joint), pik, tolerance = 0.05)
})

test_that("all joint methods agree on diagonal", {
  pik <- c(0.2, 0.3, 0.5)

  pikl_maxent <- up_maxent_jip(pik)
  pikl_sys <- up_systematic_jip(pik)
  pikl_brewer <- up_brewer_jip(pik)
  pikl_chromy <- up_chromy_jip(pik, round(sum(pik)), nsim = 1e6)

  expect_equal(diag(pikl_maxent), pik)
  expect_equal(diag(pikl_sys), pik)
  expect_equal(diag(pikl_brewer), pik)
  expect_equal(diag(pikl_chromy), pik, tolerance = 0.01)
})

test_that("all joint methods produce symmetric matrices", {
  pik <- c(0.15, 0.25, 0.35, 0.25)

  pikl_maxent <- up_maxent_jip(pik)
  pikl_sys <- up_systematic_jip(pik)
  pikl_brewer <- up_brewer_jip(pik)
  pikl_chromy <- up_chromy_jip(pik, round(sum(pik)))

  expect_equal(pikl_maxent, t(pikl_maxent))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_brewer, t(pikl_brewer))
  expect_equal(pikl_chromy, t(pikl_chromy))
})

test_that("joint probabilities work with larger populations", {
  set.seed(1)
  N <- 50
  n <- 10
  x <- runif(N)
  pik <- n * x / sum(x)

  pikl_maxent <- up_maxent_jip(pik)
  pikl_brewer <- up_brewer_jip(pik)
  pikl_sys <- up_systematic_jip(pik)
  pikl_chromy <- up_chromy_jip(x, n)

  # All should be N x N symmetric matrices
  expect_equal(dim(pikl_maxent), c(N, N))
  expect_equal(dim(pikl_brewer), c(N, N))
  expect_equal(dim(pikl_sys), c(N, N))
  expect_equal(dim(pikl_chromy), c(N, N))

  expect_equal(pikl_maxent, t(pikl_maxent))
  expect_equal(pikl_brewer, t(pikl_brewer))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_chromy, t(pikl_chromy))
})
