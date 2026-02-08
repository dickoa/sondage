test_that("up_poisson returns indices in valid range", {
  pik <- c(0.2, 0.5, 0.8)
  idx <- up_poisson(pik)
  expect_true(all(idx >= 1 & idx <= 3))
})

test_that("up_poisson has no duplicates", {
  pik <- c(0.2, 0.5, 0.8)
  idx <- up_poisson(pik)
  expect_equal(length(unique(idx)), length(idx))
})

test_that("up_poisson achieves correct inclusion probabilities", {
  pik <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  N <- length(pik)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- up_poisson(pik)
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("up_poisson handles boundary cases", {
  # All 0 - select nothing
  idx <- up_poisson(c(0, 0, 0))
  expect_length(idx, 0)

  # All 1 - select everything
  idx <- up_poisson(c(1, 1, 1))
  expect_equal(sort(idx), 1:3)
})

test_that("up_poisson is reproducible with set.seed", {
  pik <- c(0.2, 0.5, 0.8)

  set.seed(123)
  idx1 <- up_poisson(pik)
  set.seed(123)
  idx2 <- up_poisson(pik)

  expect_identical(idx1, idx2)
})

test_that("up_poisson rejects invalid input", {
  expect_error(up_poisson(c(0.5, NA)), "missing values")
  expect_error(up_poisson(c(-0.1, 0.5)), "between 0 and 1")
  expect_error(up_poisson(c(0.5, 1.1)), "between 0 and 1")
})

test_that("up_multinomial returns correct number of indices", {
  x <- c(10, 20, 30, 40)
  idx <- up_multinomial(x, n = 5)
  expect_length(idx, 5)

  idx <- up_multinomial(x, n = 1)
  expect_length(idx, 1)

  idx <- up_multinomial(x, n = 0)
  expect_length(idx, 0)
})

test_that("up_multinomial indices are in valid range", {
  x <- c(1, 2, 3, 4)
  idx <- up_multinomial(x, n = 10)
  expect_true(all(idx >= 1 & idx <= 4))
})

test_that("up_multinomial can have duplicates", {
  set.seed(42)
  x <- c(1, 2, 3, 4)
  idx <- up_multinomial(x, n = 10)
  # With 10 draws from 4 units, we expect duplicates
  expect_true(length(unique(idx)) < 10)
})

test_that("up_multinomial achieves correct proportions", {
  x <- c(1, 2, 3, 4) # Should get 10%, 20%, 30%, 40%
  n_sim <- 5000
  n_draws <- 10
  N <- length(x)

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- up_multinomial(x, n = n_draws)
    tab <- tabulate(idx, nbins = N)
    counts <- counts + tab
  }

  proportions <- counts / (n_sim * n_draws)
  expected <- x / sum(x)
  expect_equal(proportions, expected, tolerance = 0.02)
})

test_that("up_multinomial is reproducible with set.seed", {
  x <- c(1, 2, 3, 4)

  set.seed(123)
  idx1 <- up_multinomial(x, n = 10)

  set.seed(123)
  idx2 <- up_multinomial(x, n = 10)

  expect_identical(idx1, idx2)
})

test_that("up_multinomial rejects invalid x", {
  expect_error(up_multinomial(c(1, NA, 2), n = 5), "missing values")
  expect_error(up_multinomial(c(1, -1, 2), n = 5), "non-negative")
  expect_error(up_multinomial(c(0, 0, 0), n = 5), "sum of x must be positive")
  expect_error(up_multinomial(integer(0), n = 5), "empty")
  expect_error(up_multinomial("abc", n = 5), "numeric vector")
})

test_that("up_multinomial rejects invalid n", {
  x <- c(1, 2, 3, 4)
  expect_error(up_multinomial(x, n = NA), "single numeric value")
  expect_error(up_multinomial(x, n = -1), "non-negative")
  expect_error(up_multinomial(x, n = c(1, 2)), "single numeric value")
  expect_error(up_multinomial(x, n = "5"), "single numeric value")
})

test_that("up_systematic returns correct number of indices", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6) # n = 2
  idx <- up_systematic(pik)
  expect_length(idx, 2)
})

test_that("up_systematic indices are in valid range", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  idx <- up_systematic(pik)
  expect_true(all(idx >= 1 & idx <= 5))
})

test_that("up_systematic has no duplicates", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  idx <- up_systematic(pik)
  expect_equal(length(unique(idx)), length(idx))
})

test_that("up_systematic achieves correct inclusion probabilities", {
  pik <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5) # n = 2
  N <- length(pik)
  n_sim <- 10000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- up_systematic(pik)
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.02)
})

test_that("up_systematic handles certainty selections", {
  # pik near 1 should always be selected
  pik <- c(0.999999, 0.5, 0.5)
  set.seed(123)
  for (i in 1:50) {
    idx <- up_systematic(pik)
    expect_true(1 %in% idx)
  }
})

test_that("up_systematic is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)

  set.seed(123)
  idx1 <- up_systematic(pik)
  set.seed(123)
  idx2 <- up_systematic(pik)

  expect_identical(idx1, idx2)
})

test_that("up_systematic works with larger populations", {
  set.seed(42)
  N <- 500
  n <- 100
  pik <- runif(N)
  pik <- pik / sum(pik) * n

  idx <- up_systematic(pik)

  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("up_systematic rejects invalid input", {
  expect_error(up_systematic(c(0.5, NA)), "missing values")
  expect_error(up_systematic(c("a", "b")), "numeric")
  expect_error(up_systematic(numeric(0)), "empty")
})

# Non-integer parameter errors

test_that("up_multinomial rejects non-integer n", {
  expect_error(up_multinomial(c(10, 20, 30), 2.9), "not close to an integer")
})

test_that("up_multinomial silently accepts integer-like n", {
  expect_no_error(up_multinomial(c(10, 20, 30), 3.0))
})

test_that("up_multinomial rejects Inf in x", {
  expect_error(up_multinomial(c(10, Inf, 30), n = 2), "finite")
})

test_that("up_multinomial rejects NaN in x", {
  expect_error(up_multinomial(c(10, NaN, 30), n = 2), "missing values")
})

# Non-integer sum(pik) errors for up_systematic

test_that("up_systematic rejects non-integer sum(pik)", {
  expect_error(up_systematic(c(0.49, 0.49, 0.49)), "not close to an integer")
})

test_that("up_systematic rejects slightly non-integer sum(pik)", {
  expect_error(up_systematic(c(0.501, 0.25, 0.25)), "not close to an integer")
})

test_that("up_systematic silently accepts integer sum(pik)", {
  expect_no_error(up_systematic(c(0.2, 0.3, 0.5)))
})
