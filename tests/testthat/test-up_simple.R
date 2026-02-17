# ---- Poisson PPS ----

test_that("poisson returns sondage_sample object", {
  pik <- c(0.2, 0.5, 0.8)
  s <- unequal_prob_wor(pik, method = "poisson")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "poisson")
})

test_that("poisson returns indices in valid range", {
  pik <- c(0.2, 0.5, 0.8)
  idx <- unequal_prob_wor(pik, method = "poisson")$sample
  expect_true(all(idx >= 1 & idx <= 3))
})

test_that("poisson has no duplicates", {
  pik <- c(0.2, 0.5, 0.8)
  idx <- unequal_prob_wor(pik, method = "poisson")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("poisson achieves correct inclusion probabilities", {
  pik <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  N <- length(pik)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wor(pik, method = "poisson")$sample
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("poisson handles boundary cases", {
  idx <- unequal_prob_wor(c(0, 0, 0), method = "poisson")$sample
  expect_length(idx, 0)

  idx <- unequal_prob_wor(c(1, 1, 1), method = "poisson")$sample
  expect_equal(sort(idx), 1:3)
})

test_that("poisson is reproducible with set.seed", {
  pik <- c(0.2, 0.5, 0.8)

  set.seed(123)
  idx1 <- unequal_prob_wor(pik, method = "poisson")$sample
  set.seed(123)
  idx2 <- unequal_prob_wor(pik, method = "poisson")$sample

  expect_identical(idx1, idx2)
})

test_that("poisson rejects invalid input", {
  expect_error(unequal_prob_wor(c(0.5, NA), method = "poisson"),
               "missing values")
  expect_error(unequal_prob_wor(c(-0.1, 0.5), method = "poisson"),
               "between 0 and 1")
  expect_error(unequal_prob_wor(c(0.5, 1.1), method = "poisson"),
               "between 0 and 1")
})

# ---- Systematic PPS ----

test_that("systematic_pps returns sondage_sample object", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  s <- unequal_prob_wor(pik, method = "systematic")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "systematic")
})

test_that("systematic_pps returns correct number of indices", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  expect_length(unequal_prob_wor(pik, method = "systematic")$sample, 2)
})

test_that("systematic_pps indices are in valid range", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  idx <- unequal_prob_wor(pik, method = "systematic")$sample
  expect_true(all(idx >= 1 & idx <= 5))
})

test_that("systematic_pps has no duplicates", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
  idx <- unequal_prob_wor(pik, method = "systematic")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("systematic_pps achieves correct inclusion probabilities", {
  pik <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5)
  N <- length(pik)
  n_sim <- 10000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wor(pik, method = "systematic")$sample
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.02)
})

test_that("systematic_pps handles certainty selections", {
  pik <- c(0.999999, 0.5, 0.5)
  set.seed(123)
  for (i in 1:50) {
    idx <- unequal_prob_wor(pik, method = "systematic")$sample
    expect_true(1 %in% idx)
  }
})

test_that("systematic_pps is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)

  set.seed(123)
  idx1 <- unequal_prob_wor(pik, method = "systematic")$sample
  set.seed(123)
  idx2 <- unequal_prob_wor(pik, method = "systematic")$sample

  expect_identical(idx1, idx2)
})

test_that("systematic_pps works with larger populations", {
  set.seed(42)
  N <- 500
  n <- 100
  pik <- runif(N)
  pik <- pik / sum(pik) * n

  idx <- unequal_prob_wor(pik, method = "systematic")$sample

  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("systematic_pps rejects invalid input", {
  expect_error(unequal_prob_wor(c(0.5, NA), method = "systematic"),
               "missing values")
  expect_error(unequal_prob_wor(c("a", "b"), method = "systematic"),
               "numeric")
  expect_error(unequal_prob_wor(numeric(0), method = "systematic"), "empty")
})

test_that("systematic_pps rejects non-integer sum(pik)", {
  expect_error(unequal_prob_wor(c(0.49, 0.49, 0.49), method = "systematic"),
               "not close to an integer")
})

test_that("systematic_pps silently accepts integer sum(pik)", {
  expect_no_error(unequal_prob_wor(c(0.2, 0.3, 0.5), method = "systematic"))
})

# ---- Multinomial PPS ----

test_that("multinomial returns sondage_sample object", {
  x <- c(10, 20, 30, 40)
  hits <- expected_hits(x, n = 5)
  s <- unequal_prob_wr(hits, method = "multinomial")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wr")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "multinomial")
})

test_that("multinomial returns correct number of indices", {
  x <- c(10, 20, 30, 40)
  expect_length(unequal_prob_wr(expected_hits(x, n = 5),
                                method = "multinomial")$sample, 5)
  expect_length(unequal_prob_wr(expected_hits(x, n = 1),
                                method = "multinomial")$sample, 1)
})

test_that("multinomial indices are in valid range", {
  x <- c(1, 2, 3, 4)
  hits <- expected_hits(x, n = 10)
  idx <- unequal_prob_wr(hits, method = "multinomial")$sample
  expect_true(all(idx >= 1 & idx <= 4))
})

test_that("multinomial can have duplicates", {
  set.seed(42)
  x <- c(1, 2, 3, 4)
  hits <- expected_hits(x, n = 10)
  idx <- unequal_prob_wr(hits, method = "multinomial")$sample
  expect_true(length(unique(idx)) < 10)
})

test_that("multinomial achieves correct proportions", {
  x <- c(1, 2, 3, 4)
  n_sim <- 5000
  n_draws <- 10
  N <- length(x)
  hits <- expected_hits(x, n = n_draws)

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wr(hits, method = "multinomial")$sample
    tab <- tabulate(idx, nbins = N)
    counts <- counts + tab
  }

  proportions <- counts / (n_sim * n_draws)
  expected <- x / sum(x)
  expect_equal(proportions, expected, tolerance = 0.02)
})

test_that("multinomial is reproducible with set.seed", {
  x <- c(1, 2, 3, 4)
  hits <- expected_hits(x, n = 10)

  set.seed(123)
  idx1 <- unequal_prob_wr(hits, method = "multinomial")$sample
  set.seed(123)
  idx2 <- unequal_prob_wr(hits, method = "multinomial")$sample

  expect_identical(idx1, idx2)
})

test_that("multinomial rejects invalid hits", {
  expect_error(unequal_prob_wr(c(1, NA, 2), method = "multinomial"),
               "missing values")
  expect_error(unequal_prob_wr(c(1, -1, 2), method = "multinomial"),
               "non-negative")
  expect_error(unequal_prob_wr(c(0, 0, 0), method = "multinomial"),
               "sum of 'hits' must be positive")
  expect_error(unequal_prob_wr(integer(0), method = "multinomial"), "empty")
  expect_error(unequal_prob_wr("abc", method = "multinomial"),
               "numeric vector")
})

test_that("multinomial rejects non-integer sum(hits)", {
  expect_error(unequal_prob_wr(c(0.5, 0.6, 0.8), method = "multinomial"),
               "not close to an integer")
})

test_that("multinomial silently accepts integer-like sum(hits)", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 3)
  expect_no_error(unequal_prob_wr(hits, method = "multinomial"))
})

test_that("multinomial rejects Inf in hits", {
  expect_error(unequal_prob_wr(c(10, Inf, 30), method = "multinomial"),
               "finite")
})

test_that("multinomial rejects NaN in hits", {
  expect_error(unequal_prob_wr(c(10, NaN, 30), method = "multinomial"),
               "missing values")
})
