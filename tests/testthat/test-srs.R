test_that("srs returns correct number of indices", {
  idx <- srs(5, 100)
  expect_length(idx, 5)
})

test_that("srs indices are in valid range", {
  idx <- srs(10, 50)
  expect_true(all(idx >= 1 & idx <= 50))
})

test_that("srs without replacement has no duplicates", {
  idx <- srs(20, 100)
  expect_equal(length(unique(idx)), 20)
})

test_that("srs with replacement can have duplicates", {
  set.seed(42)
  # With high probability, 50 draws from 10 will have duplicates
  idx <- srs(50, 10, replace = TRUE)
  expect_true(length(unique(idx)) < 50)
})

test_that("srs achieves uniform selection", {
  n_sim <- 5000
  N <- 10
  n <- 3

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- srs(n, N)
    counts[idx] <- counts[idx] + 1
  }

  expected <- n_sim * n / N
  expect_true(all(abs(counts - expected) < 0.1 * expected))
})

test_that("srs is reproducible with set.seed", {
  set.seed(123)
  idx1 <- srs(5, 100)
  set.seed(123)
  idx2 <- srs(5, 100)
  expect_identical(idx1, idx2)
})

test_that("srs rejects invalid input", {
  expect_error(srs(10, 5), "cannot exceed") # n > N without replacement
  expect_error(srs(-1, 10), "non-negative")
  expect_error(srs(5, 0), "positive")
})

test_that("srs handles edge cases", {
  # n = 0
  idx <- srs(0, 10)
  expect_length(idx, 0)

  # n = N (select all)
  idx <- srs(5, 5)
  expect_equal(sort(idx), 1:5)
})

test_that("systematic returns correct number of indices", {
  idx <- systematic(5, 100)
  expect_length(idx, 5)
})

test_that("systematic indices are in valid range", {
  idx <- systematic(10, 50)
  expect_true(all(idx >= 1 & idx <= 50))
})

test_that("sys has no duplicates", {
  idx <- systematic(20, 100)
  expect_equal(length(unique(idx)), 20)
})

test_that("systematic produces evenly spaced indices", {
  idx <- systematic(5, 100) # interval = 20
  diffs <- diff(idx)
  # All differences should be close to 20
  expect_true(all(diffs >= 19 & diffs <= 21))
})

test_that("systematic achieves uniform selection", {
  n_sim <- 5000
  N <- 12
  n <- 3

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- systematic(n, N)
    counts[idx] <- counts[idx] + 1
  }

  expected <- n_sim * n / N
  expect_true(all(abs(counts - expected) < 0.15 * expected))
})

test_that("systematic is reproducible with set.seed", {
  set.seed(123)
  idx1 <- systematic(5, 100)
  set.seed(123)
  idx2 <- systematic(5, 100)
  expect_identical(idx1, idx2)
})

test_that("systematic rejects invalid input", {
  expect_error(systematic(10, 5), "cannot exceed")
  expect_error(systematic(-1, 10), "non-negative")
})

test_that("bernoulli returns indices in valid range", {
  idx <- bernoulli(0.5, 100)
  expect_true(all(idx >= 1 & idx <= 100))
})

test_that("bernoulli has no duplicates", {
  idx <- bernoulli(0.5, 100)
  expect_equal(length(unique(idx)), length(idx))
})

test_that("bernoulli achieves correct expected size", {
  n_sim <- 1000
  N <- 100
  p <- 0.3

  set.seed(42)
  sizes <- replicate(n_sim, length(bernoulli(p, N)))

  expect_equal(mean(sizes), N * p, tolerance = 2)
})

test_that("bernoulli achieves correct inclusion probability", {
  n_sim <- 5000
  N <- 20
  p <- 0.4

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- bernoulli(p, N)
    counts[idx] <- counts[idx] + 1
  }

  pi_hat <- counts / n_sim
  expect_equal(pi_hat, rep(p, N), tolerance = 0.03)
})

test_that("bernoulli is reproducible with set.seed", {
  set.seed(123)
  idx1 <- bernoulli(0.5, 100)
  set.seed(123)
  idx2 <- bernoulli(0.5, 100)
  expect_identical(idx1, idx2)
})

test_that("bernoulli handles boundary cases", {
  idx <- bernoulli(0, 100)
  expect_length(idx, 0)

  idx <- bernoulli(1, 100)
  expect_equal(sort(idx), 1:100)
})

test_that("bernoulli rejects invalid input", {
  expect_error(bernoulli(-0.1, 100), "probability")
  expect_error(bernoulli(1.1, 100), "probability")
})

# NA input handling

test_that("srs rejects Inf inputs", {
  expect_error(srs(Inf, 10), "finite")
  expect_error(srs(3, Inf), "finite")
  expect_error(srs(-Inf, 10), "non-negative")
})

test_that("systematic rejects Inf inputs", {
  expect_error(systematic(Inf, 10), "finite")
  expect_error(systematic(3, Inf), "finite")
})

test_that("bernoulli rejects Inf N", {
  expect_error(bernoulli(0.5, Inf), "finite")
})

test_that("srs rejects NA_real_ inputs", {
  expect_error(srs(NA_real_, 10), "non-negative")
  expect_error(srs(NA, 10), "non-negative")
  expect_error(srs(3, NA_real_), "positive")
  expect_error(srs(3, NA), "positive")
})

test_that("systematic rejects NA_real_ inputs", {
  expect_error(systematic(NA_real_, 10), "non-negative")
  expect_error(systematic(NA, 10), "non-negative")
  expect_error(systematic(3, NA_real_), "positive")
  expect_error(systematic(3, NA), "positive")
})

test_that("bernoulli rejects NA_real_ inputs", {
  expect_error(bernoulli(NA_real_, 100), "probability")
  expect_error(bernoulli(NA, 100), "probability")
  expect_error(bernoulli(0.5, NA_real_), "positive")
  expect_error(bernoulli(0.5, NA), "positive")
})

# Non-integer parameter errors

test_that("srs rejects non-integer n", {
  expect_error(srs(2.9, 10), "not close to an integer")
})

test_that("srs rejects non-integer N", {
  expect_error(srs(2, 10.7), "not close to an integer")
})

test_that("srs silently accepts integer-like doubles", {
  expect_no_error(srs(3.0, 10))
  expect_no_error(srs(3, 10.0))
})

test_that("systematic rejects non-integer n", {
  expect_error(systematic(2.9, 10), "not close to an integer")
})

test_that("bernoulli rejects non-integer N", {
  expect_error(bernoulli(0.5, 10.7), "not close to an integer")
})
