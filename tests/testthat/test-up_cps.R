test_that("cps returns sondage_sample object", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "cps")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "cps")
  expect_equal(s$n, 2L)
  expect_equal(s$N, 4L)
  expect_equal(s$pik, pik)
})

test_that("cps returns correct number of indices", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  expect_length(unequal_prob_wor(pik, method = "cps")$sample, 2)
})

test_that("cps indices are in valid range", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  idx <- unequal_prob_wor(pik, method = "cps")$sample
  expect_true(all(idx >= 1 & idx <= 4))
})

test_that("cps has no duplicates", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  idx <- unequal_prob_wor(pik, method = "cps")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("cps achieves correct inclusion probabilities", {
  pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  N <- length(pik)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wor(pik, method = "cps")$sample
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("cps handles certainty selections", {
  pik <- c(0.999999, 0.5, 0.5)
  set.seed(123)
  for (i in 1:50) {
    idx <- unequal_prob_wor(pik, method = "cps")$sample
    expect_true(1 %in% idx)
  }
})

test_that("cps is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5)

  set.seed(999)
  idx1 <- unequal_prob_wor(pik, method = "cps")$sample
  set.seed(999)
  idx2 <- unequal_prob_wor(pik, method = "cps")$sample

  expect_identical(idx1, idx2)
})

test_that("cps batch mode returns design object with matrix sample", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "cps", nrep = 100)

  expect_s3_class(s, "sondage_sample")
  expect_true(is.matrix(s$sample))
  expect_equal(nrow(s$sample), 2)
  expect_equal(ncol(s$sample), 100)
  expect_equal(s$pik, pik)
  expect_equal(s$method, "cps")
})

test_that("cps batch achieves correct inclusion probabilities", {
  skip_on_cran()
  pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  N <- length(pik)

  set.seed(42)
  s <- unequal_prob_wor(pik, method = "cps", nrep = 5000)
  samples <- s$sample

  counts <- integer(N)
  for (j in 1:ncol(samples)) {
    counts[samples[, j]] <- counts[samples[, j]] + 1
  }
  pi_hat <- counts / ncol(samples)

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("cps batch is faster than loop", {
  skip_on_cran()
  pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  nrep <- 500

  t_batch <- system.time({
    set.seed(42)
    unequal_prob_wor(pik, method = "cps", nrep = nrep)
  })[3]

  t_loop <- system.time({
    set.seed(42)
    replicate(nrep, unequal_prob_wor(pik, method = "cps")$sample)
  })[3]

  expect_true(t_batch < t_loop / 2 || t_batch < 0.1)
})

test_that("cps works with larger populations", {
  set.seed(123)
  N <- 500
  n <- 100
  pik <- runif(N)
  pik <- pik / sum(pik) * n

  idx <- unequal_prob_wor(pik, method = "cps")$sample
  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("cps rejects invalid input", {
  expect_error(unequal_prob_wor(c(0.5, NA, 0.5), method = "cps"),
               "missing values")
  expect_error(unequal_prob_wor(c("a", "b"), method = "cps"), "numeric")
  expect_error(unequal_prob_wor(numeric(0), method = "cps"), "empty")
})

test_that("cps rejects invalid nrep", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "cps", nrep = 0),
               "at least 1")
})

test_that("cps rejects non-integer sum(pik)", {
  expect_error(unequal_prob_wor(c(0.49, 0.49, 0.49), method = "cps"),
               "not close to an integer")
})

test_that("cps rejects slightly non-integer sum(pik)", {
  expect_error(unequal_prob_wor(c(0.501, 0.25, 0.25), method = "cps"),
               "not close to an integer")
})

test_that("cps silently accepts integer sum(pik)", {
  expect_no_error(unequal_prob_wor(c(0.2, 0.4, 0.6, 0.8), method = "cps"))
})

test_that("cps rejects non-integer nrep", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "cps", nrep = 2.9),
               "not close to an integer")
})

test_that("cps rejects Inf nrep", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "cps", nrep = Inf),
               "finite")
})

test_that("cps handles N=1 census", {
  s <- unequal_prob_wor(1.0, method = "cps")
  expect_equal(s$sample, 1L)
})

test_that("cps handles all-certainty pik", {
  s <- unequal_prob_wor(c(1.0, 1.0, 1.0), method = "cps")
  expect_equal(sort(s$sample), 1:3)
})

test_that("unequal_prob_wor dispatches to cps by default", {
  set.seed(42)
  s <- unequal_prob_wor(c(0.2, 0.3, 0.5))
  expect_equal(s$method, "cps")
  expect_s3_class(s, "wor")
})

test_that("cps handles high sampling fraction via complement path", {
  skip_on_cran()
  set.seed(321)

  N <- 400L
  n <- 360L
  pik <- inclusion_prob(rgamma(N, shape = 0.7, rate = 1), n = n)

  # Single draws
  for (i in 1:20) {
    s <- unequal_prob_wor(pik, method = "cps")
    expect_length(s$sample, n)
    expect_true(all(s$sample >= 1L & s$sample <= N))
    expect_false(anyDuplicated(s$sample) > 0)
  }

  # Batch draws
  sb <- unequal_prob_wor(pik, method = "cps", nrep = 100)
  expect_true(is.matrix(sb$sample))
  expect_equal(dim(sb$sample), c(n, 100))
  for (j in 1:100) {
    col <- sb$sample[, j]
    expect_true(all(col >= 1L & col <= N))
    expect_false(anyDuplicated(col) > 0)
  }
})

test_that("cps complement path achieves correct inclusion probs", {
  skip_on_cran()
  set.seed(42)

  N <- 100L
  n <- 80L
  pik <- inclusion_prob(rgamma(N, shape = 0.7, rate = 1), n = n)
  n_sim <- 5000

  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wor(pik, method = "cps")$sample
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("cps complement path is faster for high sampling fractions", {
  skip_on_cran()
  set.seed(123)

  N <- 1000L
  n_high <- 900L
  n_low <- 100L
  pik_high <- inclusion_prob(rgamma(N, shape = 0.7, rate = 1), n = n_high)
  pik_low <- inclusion_prob(rgamma(N, shape = 0.7, rate = 1), n = n_low)

  # High fraction (n=900, complement=100) should be comparable to low (n=100)
  t_high <- system.time({
    set.seed(42)
    unequal_prob_wor(pik_high, method = "cps", nrep = 50)
  })[3]

  t_low <- system.time({
    set.seed(42)
    unequal_prob_wor(pik_low, method = "cps", nrep = 50)
  })[3]

  # High fraction should take at most 3x the time of low (without complement
  # path it would be 9x+ due to O(N*n) DP). Allow generous margin for CI.
  expect_true(t_high < t_low * 4 || t_high < 0.5)
})

test_that("cps stress: repeated medium-large draws do not crash", {
  skip_on_cran()
  set.seed(123)
  N <- 2500L
  n <- 1250L
  pik <- inclusion_prob(rgamma(N, shape = 0.7, rate = 1), n = n)

  # Single draws
  for (i in 1:40) {
    s <- unequal_prob_wor(pik, method = "cps")
    expect_length(s$sample, n)
    expect_true(all(s$sample >= 1L & s$sample <= N))
    expect_false(anyDuplicated(s$sample) > 0)
  }

  # Batch draw
  sb <- unequal_prob_wor(pik, method = "cps", nrep = 40)
  expect_true(is.matrix(sb$sample))
  expect_equal(nrow(sb$sample), n)
  expect_equal(ncol(sb$sample), 40L)
  for (j in 1:40) {
    col <- sb$sample[, j]
    expect_true(all(col >= 1L & col <= N))
    expect_false(anyDuplicated(col) > 0)
  }
})
