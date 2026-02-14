test_that("chromy returns sondage_sample object", {
  x <- c(40, 80, 50, 60, 70)
  hits <- expected_hits(x, n = 3)
  s <- unequal_prob_wr(hits, method = "chromy")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wr")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "chromy")
  expect_equal(s$n, 3L)
  expect_equal(s$N, 5L)
  expect_equal(s$prob, x / sum(x))
})

test_that("chromy returns correct sample size", {
  x <- c(40, 80, 50, 60, 70)
  hits <- expected_hits(x, n = 3)
  expect_length(unequal_prob_wr(hits, method = "chromy")$sample, 3)
})

test_that("chromy indices are in valid range", {
  x <- c(40, 80, 50, 60, 70)
  hits <- expected_hits(x, n = 3)
  idx <- unequal_prob_wr(hits, method = "chromy")$sample
  expect_true(all(idx >= 1 & idx <= 5))
})

test_that("chromy has no duplicates when WOR", {
  x <- c(40, 80, 50, 60, 70)
  hits <- expected_hits(x, n = 3)
  idx <- unequal_prob_wr(hits, method = "chromy")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("chromy achieves correct expected hits", {
  x <- c(40, 80, 50, 60, 70)
  n <- 3
  N <- length(x)
  expected <- n * x / sum(x)
  hits <- expected_hits(x, n = n)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wr(hits, method = "chromy")$sample
    for (j in idx) {
      counts[j] <- counts[j] + 1
    }
  }
  empirical <- counts / n_sim

  expect_equal(empirical, expected, tolerance = 0.03)
})

test_that("chromy returns sorted indices", {
  x <- c(20, 30, 50, 40, 60)
  hits <- expected_hits(x, n = 2)
  for (i in 1:20) {
    idx <- unequal_prob_wr(hits, method = "chromy")$sample
    expect_equal(idx, sort(idx))
  }
})

test_that("chromy is reproducible with set.seed", {
  x <- c(20, 40, 40)
  hits <- expected_hits(x, n = 1)

  set.seed(999)
  idx1 <- unequal_prob_wr(hits, method = "chromy")$sample
  set.seed(999)
  idx2 <- unequal_prob_wr(hits, method = "chromy")$sample

  expect_identical(idx1, idx2)
})

test_that("chromy works with larger populations", {
  set.seed(456)
  N <- 100
  n <- 20
  x <- runif(N, 1, 10)
  hits <- expected_hits(x, n = n)

  idx <- unequal_prob_wr(hits, method = "chromy")$sample
  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("chromy handles minimum replacement", {
  x <- c(40, 80, 50, 60, 70)
  n <- 10
  expected <- n * x / sum(x)
  hits <- expected_hits(x, n = n)

  set.seed(123)
  idx <- unequal_prob_wr(hits, method = "chromy")$sample
  expect_length(idx, n)

  realized <- tabulate(idx, nbins = 5)
  for (k in 1:5) {
    expect_true(realized[k] %in% c(floor(expected[k]), ceiling(expected[k])),
      label = sprintf("Unit %d: hits=%d, expected=%.2f",
                      k, realized[k], expected[k]))
  }
})

test_that("chromy minimum replacement achieves correct expected hits", {
  skip_on_cran()

  x <- c(40, 80, 50, 60, 70)
  n <- 10
  N <- length(x)
  expected <- n * x / sum(x)
  hits <- expected_hits(x, n = n)
  n_sim <- 3000

  set.seed(42)
  total_hits <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wr(hits, method = "chromy")$sample
    for (j in idx) {
      total_hits[j] <- total_hits[j] + 1
    }
  }
  empirical <- total_hits / n_sim

  expect_equal(empirical, expected, tolerance = 0.05)
})

test_that("chromy has all positive joint probabilities", {
  skip_on_cran()

  x <- c(40, 80, 50, 60, 70)
  N <- length(x)
  hits <- expected_hits(x, n = 3)
  n_sim <- 5000

  pair_seen <- matrix(FALSE, N, N)

  set.seed(42)
  for (i in 1:n_sim) {
    idx <- unique(unequal_prob_wr(hits, method = "chromy")$sample)
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

test_that("chromy rejects invalid input", {
  expect_error(unequal_prob_wr(c(NA, 0.5, 0.5), method = "chromy"),
               "missing values")
  expect_error(unequal_prob_wr(c("a", "b"), method = "chromy"), "numeric")
  expect_error(unequal_prob_wr(numeric(0), method = "chromy"), "empty")
  expect_error(unequal_prob_wr(c(-1, 2, 3), method = "chromy"), "non-negative")
  expect_error(unequal_prob_wr(c(0, 0, 0), method = "chromy"), "positive")
})

test_that("chromy rejects Inf in hits", {
  expect_error(unequal_prob_wr(c(1, Inf), method = "chromy"), "finite")
})

test_that("chromy rejects NaN in hits", {
  expect_error(unequal_prob_wr(c(1, NaN), method = "chromy"), "missing values")
})

test_that("chromy rejects non-integer sum(hits)", {
  expect_error(unequal_prob_wr(c(0.5, 0.6, 0.8), method = "chromy"),
               "not close to an integer")
})

test_that("chromy silently accepts integer-like sum(hits)", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 2)
  expect_no_error(unequal_prob_wr(hits, method = "chromy"))
})

test_that("unequal_prob_wr dispatches to chromy by default", {
  set.seed(42)
  hits <- c(0.4, 0.8, 0.8)
  s <- unequal_prob_wr(hits)
  expect_equal(s$method, "chromy")
  expect_s3_class(s, "wr")
})
