# Spatially correlated Poisson sampling (Grafstrom, 2012)

test_that("scps returns a fixed-size spread-only design", {
  set.seed(101)
  pik <- c(0.2, 0.4, 0.6, 0.8)
  spread <- matrix(runif(8), ncol = 2)
  s <- balanced_wor(pik, spread = spread, method = "scps")

  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "balanced")
  expect_s3_class(s, "unequal_prob")
  expect_s3_class(s, "wor")
  expect_equal(s$method, "scps")
  expect_equal(s$n, 2L)
  expect_equal(s$N, 4L)
  expect_equal(s$pik, pik)
  expect_true(s$fixed_size)
  expect_length(s$sample, 2L)
  expect_equal(anyDuplicated(s$sample), 0L)
  expect_false(is.unsorted(s$sample))
})

test_that("scps respects unequal first-order inclusion probabilities", {
  set.seed(102)
  pik <- c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9)
  spread <- matrix(c(pik, runif(20)), ncol = 3)
  nrep <- 10000L
  s <- balanced_wor(
    pik,
    spread = spread,
    method = "scps",
    nrep = nrep
  )
  empirical <- tabulate(as.integer(s$sample), nbins = length(pik)) / nrep

  expect_true(max(abs(empirical - pik)) < 0.02)
  expect_equal(dim(s$sample), c(5L, nrep))
})

test_that("scps handles equal-distance groups without marginal bias", {
  set.seed(103)
  spread <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  storage.mode(spread) <- "double"
  pik <- rep(0.2, nrow(spread))
  nrep <- 5000L
  s <- balanced_wor(
    pik,
    spread = spread,
    method = "scps",
    nrep = nrep
  )
  empirical <- tabulate(as.integer(s$sample), nbins = length(pik)) / nrep

  expect_true(max(abs(empirical - pik)) < 0.03)
})

test_that("scps spreads samples better than srs", {
  set.seed(104)
  N <- 100L
  spread <- matrix(runif(2L * N), ncol = 2)
  pik <- rep(0.2, N)
  mean_nearest_distance <- function(index) {
    distance <- as.matrix(dist(spread[index, , drop = FALSE]))
    diag(distance) <- Inf
    mean(apply(distance, 1L, min))
  }

  scps_spread <- mean(replicate(
    30L,
    mean_nearest_distance(
      balanced_wor(pik, spread = spread, method = "scps")$sample
    )
  ))
  srs_spread <- mean(replicate(
    30L,
    mean_nearest_distance(sample.int(N, 20L))
  ))

  expect_gt(scps_spread, srs_spread)
})

test_that("scps preserves separated coincident-pair sample sizes", {
  set.seed(105)
  groups <- 10L
  centres <- cbind(10 * seq_len(groups), 10 * seq_len(groups))
  spread <- centres[rep(seq_len(groups), each = 2L), ]
  pik <- rep(0.5, 2L * groups)

  for (draw in seq_len(20L)) {
    index <- balanced_wor(pik, spread = spread, method = "scps")$sample
    pair <- (index + 1L) %/% 2L
    expect_equal(sort(pair), seq_len(groups))
  }
})

test_that("scps batch draws are reproducible and pass through certainty units", {
  pik <- c(1, 0, rep(0.5, 8))
  spread <- matrix(seq_len(20) / 20, ncol = 2)
  set.seed(106)
  first <- balanced_wor(
    pik,
    spread = spread,
    method = "scps",
    nrep = 20L
  )$sample
  set.seed(106)
  second <- balanced_wor(
    pik,
    spread = spread,
    method = "scps",
    nrep = 20L
  )$sample

  expect_identical(first, second)
  expect_equal(dim(first), c(5L, 20L))
  expect_true(all(apply(first, 2L, function(x) 1L %in% x)))
  expect_false(any(first == 2L))
})

test_that("scps preserves probability mass over many large-N draws", {
  # Regression: boundary snapping used to accumulate a few ulps per step.
  # Eventually the final two active units could appear to have maximal-weight
  # capacity just below one, even though their probability sum was one.
  set.seed(2600)
  N <- 2500L
  spread <- matrix(runif(2L * N), N, 2L)
  pik <- rep(0.1, N)
  draws <- balanced_wor(
    pik,
    spread = spread,
    method = "scps",
    nrep = 150L
  )$sample

  expect_equal(dim(draws), c(250L, 150L))
  expect_true(all(apply(draws, 2L, anyDuplicated) == 0L))
})

test_that("scps enforces its spread-only capability contract", {
  pik <- rep(0.5, 10)
  spread <- matrix(runif(20), ncol = 2)

  expect_error(
    balanced_wor(pik, aux = matrix(1, 10), spread = spread, method = "scps"),
    "does not use auxiliary balancing variables"
  )
  expect_error(
    balanced_wor(
      pik,
      strata = rep(1:2, 5),
      spread = spread,
      method = "scps"
    ),
    "does not support 'strata'"
  )
  expect_error(
    balanced_wor(pik, method = "scps"),
    "method 'scps' requires 'spread'"
  )
})

test_that("scps metadata and unsupported joint probabilities are explicit", {
  spec <- method_spec("scps")
  expect_equal(spec$type, "balanced")
  expect_true(spec$fixed_size)
  expect_equal(spec$variance_family, "unsupported")
  expect_false(spec$supports_aux)
  expect_false(spec$supports_strata)
  expect_true(spec$supports_spread)

  set.seed(107)
  s <- balanced_wor(
    rep(0.5, 10),
    spread = matrix(runif(20), ncol = 2),
    method = "scps"
  )
  expect_error(
    joint_inclusion_prob(s),
    "joint_inclusion_prob not implemented for method 'scps'"
  )
  expect_error(sampling_cov(s))
  expect_equal(inclusion_prob(s), rep(0.5, 10))
  expect_output(print(s), "Balanced WOR \\[scps\\]")
})

test_that("scps is a reserved built-in method name", {
  expect_error(
    register_method(
      "scps",
      type = "balanced",
      sample_fn = function(pik, n = NULL, aux = NULL, ...) seq_len(n)
    ),
    "'scps' is a built-in method and cannot be overridden"
  )
})
