test_that("brewer returns sondage_sample object", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "brewer")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "brewer")
  expect_equal(s$n, 2L)
  expect_equal(s$N, 4L)
})

test_that("brewer returns correct number of indices", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  expect_length(unequal_prob_wor(pik, method = "brewer")$sample, 2)
})

test_that("brewer indices are in valid range", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  idx <- unequal_prob_wor(pik, method = "brewer")$sample
  expect_true(all(idx >= 1 & idx <= 4))
})

test_that("brewer has no duplicates", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  idx <- unequal_prob_wor(pik, method = "brewer")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("brewer achieves correct inclusion probabilities", {
  pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  N <- length(pik)
  n_sim <- 5000

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- unequal_prob_wor(pik, method = "brewer")$sample
    counts[idx] <- counts[idx] + 1
  }
  pi_hat <- counts / n_sim

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("brewer handles certainty selections", {
  pik <- c(0.999999, 0.5, 0.5)
  set.seed(123)
  for (i in 1:50) {
    idx <- unequal_prob_wor(pik, method = "brewer")$sample
    expect_true(1 %in% idx)
  }
})

test_that("brewer excludes near-zero units", {
  pik <- c(1e-8, 0.5, 0.5)
  set.seed(123)
  for (i in 1:50) {
    idx <- unequal_prob_wor(pik, method = "brewer")$sample
    expect_false(1 %in% idx)
  }
})

test_that("brewer is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5)

  set.seed(999)
  idx1 <- unequal_prob_wor(pik, method = "brewer")$sample
  set.seed(999)
  idx2 <- unequal_prob_wor(pik, method = "brewer")$sample

  expect_identical(idx1, idx2)
})

test_that("brewer works with larger populations", {
  set.seed(456)
  N <- 100
  n <- 20
  x <- runif(N)
  pik <- n * x / sum(x)

  idx <- unequal_prob_wor(pik, method = "brewer")$sample
  expect_length(idx, n)
  expect_true(all(idx >= 1 & idx <= N))
})

test_that("brewer rejects invalid input", {
  expect_error(unequal_prob_wor(c(0.5, NA, 0.5), method = "brewer"),
               "missing values")
  expect_error(unequal_prob_wor(c("a", "b"), method = "brewer"), "numeric")
  expect_error(unequal_prob_wor(numeric(0), method = "brewer"), "empty")
})

test_that("brewer rejects non-integer sum(pik)", {
  expect_error(unequal_prob_wor(c(0.49, 0.49, 0.49), method = "brewer"),
               "not close to an integer")
})

test_that("brewer rejects slightly non-integer sum(pik)", {
  expect_error(unequal_prob_wor(c(0.501, 0.25, 0.25), method = "brewer"),
               "not close to an integer")
})

test_that("brewer silently accepts integer sum(pik)", {
  expect_no_error(unequal_prob_wor(c(0.2, 0.4, 0.6, 0.8), method = "brewer"))
})

test_that("brewer handles N=1 census", {
  s <- unequal_prob_wor(1.0, method = "brewer")
  expect_equal(s$sample, 1L)
})

test_that("brewer handles all-certainty pik", {
  s <- unequal_prob_wor(c(1.0, 1.0, 1.0), method = "brewer")
  expect_equal(sort(s$sample), 1:3)
})

test_that("brewer denom clamp handles near-certainty units with tight eps", {
  # With small eps, high-pik units enter the active pool and the Brewer

  # denominator (m - pk*r) can become extremely small. A previous guard
  # that zeroed the draw probability when denom < 1e-10 would incorrectly
  # skip these units. The clamp to 1e-15 ensures they dominate selection.
  delta <- 2e-11
  pik <- c(1 - delta, 1 - delta, delta / 2, delta / 2 + delta)

  set.seed(42)
  counts <- integer(4)
  nrep <- 200
  for (i in seq_len(nrep)) {
    s <- unequal_prob_wor(pik, method = "brewer", eps = 1e-12)
    counts[s$sample] <- counts[s$sample] + 1L
  }
  # Near-certainty units must be selected every time
  expect_equal(counts[1], nrep)
  expect_equal(counts[2], nrep)
})
