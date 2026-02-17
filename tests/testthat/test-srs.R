test_that("srs wor returns sondage_sample object", {
  s <- equal_prob_wor(100, 5)
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "equal_prob")
  expect_equal(s$method, "srs")
  expect_equal(s$n, 5L)
  expect_equal(s$N, 100L)
  expect_equal(s$pik, rep(5 / 100, 100))
})

test_that("srs wor returns correct number of indices", {
  s <- equal_prob_wor(100, 5)
  expect_length(s$sample, 5)
})

test_that("srs wor indices are in valid range", {
  idx <- equal_prob_wor(50, 10)$sample
  expect_true(all(idx >= 1 & idx <= 50))
})

test_that("srs wor has no duplicates", {
  idx <- equal_prob_wor(100, 20)$sample
  expect_equal(length(unique(idx)), 20)
})

test_that("srs wr can have duplicates", {
  set.seed(42)
  idx <- equal_prob_wr(10, 50)$sample
  expect_true(length(unique(idx)) < 50)
})

test_that("srs wr returns correct object", {
  s <- equal_prob_wr(10, 3)
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wr")
  expect_s3_class(s, "equal_prob")
  expect_equal(s$method, "srs")
  expect_equal(s$n, 3L)
  expect_equal(s$N, 10L)
  expect_equal(s$prob, rep(1 / 10, 10))
})

test_that("srs achieves uniform selection", {
  n_sim <- 5000
  N <- 10
  n <- 3

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- equal_prob_wor(N, n)$sample
    counts[idx] <- counts[idx] + 1
  }

  expected <- n_sim * n / N
  expect_true(all(abs(counts - expected) < 0.1 * expected))
})

test_that("srs is reproducible with set.seed", {
  set.seed(123)
  idx1 <- equal_prob_wor(100, 5)$sample
  set.seed(123)
  idx2 <- equal_prob_wor(100, 5)$sample
  expect_identical(idx1, idx2)
})

test_that("srs rejects invalid input", {
  expect_error(equal_prob_wor(5, 10), "cannot exceed")
  expect_error(equal_prob_wor(10, -1), "non-negative")
  expect_error(equal_prob_wor(0, 5), "positive")
})

test_that("srs handles edge cases", {
  idx <- equal_prob_wor(10, 0)$sample
  expect_length(idx, 0)

  idx <- equal_prob_wor(5, 5)$sample
  expect_equal(sort(idx), 1:5)
})

test_that("systematic returns sondage_sample object", {
  s <- equal_prob_wor(100, 5, method = "systematic")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "equal_prob")
  expect_equal(s$method, "systematic")
})

test_that("systematic returns correct number of indices", {
  expect_length(equal_prob_wor(100, 5, method = "systematic")$sample, 5)
})

test_that("systematic indices are in valid range", {
  idx <- equal_prob_wor(50, 10, method = "systematic")$sample
  expect_true(all(idx >= 1 & idx <= 50))
})

test_that("systematic has no duplicates", {
  idx <- equal_prob_wor(100, 20, method = "systematic")$sample
  expect_equal(length(unique(idx)), 20)
})

test_that("systematic produces evenly spaced indices", {
  idx <- equal_prob_wor(100, 5, method = "systematic")$sample
  diffs <- diff(idx)
  expect_true(all(diffs >= 19 & diffs <= 21))
})

test_that("systematic achieves uniform selection", {
  n_sim <- 5000
  N <- 12
  n <- 3

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- equal_prob_wor(N, n, method = "systematic")$sample
    counts[idx] <- counts[idx] + 1
  }

  expected <- n_sim * n / N
  expect_true(all(abs(counts - expected) < 0.15 * expected))
})

test_that("systematic is reproducible with set.seed", {
  set.seed(123)
  idx1 <- equal_prob_wor(100, 5, method = "systematic")$sample
  set.seed(123)
  idx2 <- equal_prob_wor(100, 5, method = "systematic")$sample
  expect_identical(idx1, idx2)
})

test_that("systematic rejects invalid input", {
  expect_error(equal_prob_wor(5, 10, method = "systematic"), "cannot exceed")
  expect_error(equal_prob_wor(10, -1, method = "systematic"), "non-negative")
})

test_that("bernoulli returns sondage_sample object", {
  s <- equal_prob_wor(100, 50, method = "bernoulli")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "equal_prob")
  expect_equal(s$method, "bernoulli")
  expect_equal(s$pik, rep(0.5, 100))
})

test_that("bernoulli returns indices in valid range", {
  idx <- equal_prob_wor(100, 50, method = "bernoulli")$sample
  expect_true(all(idx >= 1 & idx <= 100))
})

test_that("bernoulli has no duplicates", {
  idx <- equal_prob_wor(100, 50, method = "bernoulli")$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("bernoulli achieves correct expected size", {
  n_sim <- 1000
  N <- 100
  n <- 30

  set.seed(42)
  sizes <- replicate(
    n_sim,
    length(equal_prob_wor(N, n, method = "bernoulli")$sample)
  )

  expect_equal(mean(sizes), n, tolerance = 2)
})

test_that("bernoulli achieves correct inclusion probability", {
  n_sim <- 5000
  N <- 20
  n <- 8 # p = 0.4

  set.seed(42)
  counts <- integer(N)
  for (i in 1:n_sim) {
    idx <- equal_prob_wor(N, n, method = "bernoulli")$sample
    counts[idx] <- counts[idx] + 1
  }

  pi_hat <- counts / n_sim
  expect_equal(pi_hat, rep(n / N, N), tolerance = 0.03)
})

test_that("bernoulli is reproducible with set.seed", {
  set.seed(123)
  idx1 <- equal_prob_wor(100, 50, method = "bernoulli")$sample
  set.seed(123)
  idx2 <- equal_prob_wor(100, 50, method = "bernoulli")$sample
  expect_identical(idx1, idx2)
})

test_that("bernoulli handles boundary cases", {
  idx <- equal_prob_wor(100, 0, method = "bernoulli")$sample
  expect_length(idx, 0)

  idx <- equal_prob_wor(100, 100, method = "bernoulli")$sample
  expect_equal(sort(idx), 1:100)
})

test_that("bernoulli rejects invalid input", {
  expect_error(equal_prob_wor(100, -1, method = "bernoulli"), "non-negative")
  expect_error(
    equal_prob_wor(100, 101, method = "bernoulli"),
    "cannot exceed"
  )
})

test_that("srs rejects Inf inputs", {
  expect_error(equal_prob_wor(10, Inf), "cannot exceed")
  expect_error(equal_prob_wor(Inf, 3), "finite")
  expect_error(equal_prob_wor(10, -Inf), "non-negative")
})

test_that("systematic rejects Inf inputs", {
  expect_error(equal_prob_wor(10, Inf, method = "systematic"), "cannot exceed")
  expect_error(equal_prob_wor(Inf, 3, method = "systematic"), "finite")
})

test_that("srs rejects NA inputs", {
  expect_error(equal_prob_wor(10, NA_real_), "non-negative")
  expect_error(equal_prob_wor(10, NA), "non-negative")
  expect_error(equal_prob_wor(NA_real_, 3), "positive")
  expect_error(equal_prob_wor(NA, 3), "positive")
})

test_that("systematic rejects NA inputs", {
  expect_error(
    equal_prob_wor(10, NA_real_, method = "systematic"),
    "non-negative"
  )
  expect_error(equal_prob_wor(10, NA, method = "systematic"), "non-negative")
  expect_error(equal_prob_wor(NA_real_, 3, method = "systematic"), "positive")
  expect_error(equal_prob_wor(NA, 3, method = "systematic"), "positive")
})

test_that("bernoulli rejects NA inputs", {
  expect_error(
    equal_prob_wor(100, NA_real_, method = "bernoulli"),
    "non-negative"
  )
  expect_error(equal_prob_wor(100, NA, method = "bernoulli"), "non-negative")
  expect_error(equal_prob_wor(NA_real_, 50, method = "bernoulli"), "positive")
  expect_error(equal_prob_wor(NA, 50, method = "bernoulli"), "positive")
})

test_that("srs rejects non-integer n", {
  expect_error(equal_prob_wor(10, 2.9), "not close to an integer")
})

test_that("srs rejects non-integer N", {
  expect_error(equal_prob_wor(10.7, 2), "not close to an integer")
})

test_that("srs silently accepts integer-like doubles", {
  expect_no_error(equal_prob_wor(10, 3.0))
  expect_no_error(equal_prob_wor(10.0, 3))
})

test_that("systematic rejects non-integer n", {
  expect_error(
    equal_prob_wor(10, 2.9, method = "systematic"),
    "not close to an integer"
  )
})

test_that("bernoulli rejects non-integer N", {
  expect_error(
    equal_prob_wor(10.7, 5, method = "bernoulli"),
    "not close to an integer"
  )
})
