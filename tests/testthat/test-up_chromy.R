test_that("up_chromy returns correct sample size", {
  x <- c(40, 80, 50, 60, 70)
  idx <- up_chromy(x, n = 3)
  expect_length(idx, 3)
})

test_that("up_chromy indices are in valid range", {
  x <- c(40, 80, 50, 60, 70)
  idx <- up_chromy(x, n = 3)
  expect_true(all(idx >= 1 & idx <= 5))
})

test_that("up_chromy has no duplicates when WOR", {
  # When n * max(x) / sum(x) < 1, no repeats
  x <- c(40, 80, 50, 60, 70)  # max expected = 3 * 80/300 = 0.8
  idx <- up_chromy(x, n = 3)
  expect_equal(length(unique(idx)), length(idx))
})

test_that("up_chromy achieves correct expected hits", {
  x <- c(40, 80, 50, 60, 70)
  n <- 3
  N <- length(x)
  expected <- n * x / sum(x)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- up_chromy(x, n)
    for (j in idx) {
      counts[j] <- counts[j] + 1
    }
  }
  empirical <- counts / n_sim

  expect_equal(empirical, expected, tolerance = 0.03)
})

test_that("up_chromy returns sorted indices", {
  x <- c(20, 30, 50, 40, 60)
  for (i in 1:20) {
    idx <- up_chromy(x, n = 2)
    expect_equal(idx, sort(idx))
  }
})

test_that("up_chromy is reproducible with set.seed", {
  x <- c(20, 40, 40)

  set.seed(999)
  idx1 <- up_chromy(x, n = 1)
  set.seed(999)
  idx2 <- up_chromy(x, n = 1)

  expect_identical(idx1, idx2)
})

test_that("up_chromy works with larger populations", {
  set.seed(456)
  N <- 100
  n <- 20
  x <- runif(N, 1, 10)

  idx <- up_chromy(x, n)
  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("up_chromy handles minimum replacement", {
  # Large n causes expected hits > 1
  x <- c(40, 80, 50, 60, 70)
  n <- 10
  expected <- n * x / sum(x)
  # expected = c(1.33, 2.67, 1.67, 2.0, 2.33)

  set.seed(123)
  idx <- up_chromy(x, n)
  expect_length(idx, n)

  # Check hits are floor or ceiling of expected
  hits <- tabulate(idx, nbins = 5)
  for (k in 1:5) {
    expect_true(hits[k] %in% c(floor(expected[k]), ceiling(expected[k])),
      label = sprintf("Unit %d: hits=%d, expected=%.2f", k, hits[k], expected[k]))
  }
})

test_that("up_chromy minimum replacement achieves correct expected hits", {
  skip_on_cran()

  x <- c(40, 80, 50, 60, 70)
  n <- 10
  N <- length(x)
  expected <- n * x / sum(x)
  n_sim <- 3000

  set.seed(42)
  total_hits <- integer(N)
  for (i in 1:n_sim) {
    idx <- up_chromy(x, n)
    for (j in idx) {
      total_hits[j] <- total_hits[j] + 1
    }
  }
  empirical <- total_hits / n_sim

  expect_equal(empirical, expected, tolerance = 0.05)
})

test_that("up_chromy has all positive joint probabilities", {
  skip_on_cran()

  x <- c(40, 80, 50, 60, 70)
  N <- length(x)
  n_sim <- 5000

  pair_seen <- matrix(FALSE, N, N)

  set.seed(42)
  for (i in 1:n_sim) {
    idx <- unique(up_chromy(x, n = 3))  # unique for WOR
    for (a in idx) {
      for (b in idx) {
        pair_seen[a, b] <- TRUE
      }
    }
  }

  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        expect_true(pair_seen[i, j],
          label = sprintf("Pair (%d,%d) never selected", i, j))
      }
    }
  }
})

test_that("up_chromy rejects invalid input", {
  expect_error(up_chromy(c(10, NA, 20), n = 1), "missing values")
  expect_error(up_chromy(c("a", "b"), n = 1), "numeric")
  expect_error(up_chromy(numeric(0), n = 1), "empty")
  expect_error(up_chromy(c(-1, 2, 3), n = 1), "non-negative")
  expect_error(up_chromy(c(0, 0, 0), n = 1), "positive")
  expect_error(up_chromy(c(1, 2, 3), n = 0), "at least 1")
})
