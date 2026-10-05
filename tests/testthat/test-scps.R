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
  expect_true(spec$supports_prn)

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

test_that("scps tolerates a pik total off an integer by the accepted residue", {
  # Regression: the residue reached the last step amplified by 1 / p, so
  # draws failed as "numerically infeasible" (all of them at 1e-11).
  set.seed(108)
  N <- 30L
  spread <- matrix(runif(2L * N), N)
  for (delta in c(-9e-11, -1e-13, 1e-15, 1e-13, 9e-11)) {
    pik <- inclusion_prob(rexp(N) + 0.01, 6)
    free <- which(pik < 0.9)[1]
    pik[free] <- pik[free] + (6 + delta - sum(pik))
    sizes <- replicate(20L, {
      length(balanced_wor(pik, spread = spread, method = "scps")$sample)
    })
    prn_sizes <- replicate(20L, {
      u <- runif(N)
      length(balanced_wor(pik, spread = spread, method = "scps", prn = u)$sample)
    })
    expect_true(all(c(sizes, prn_sizes) == 6L), info = paste("delta =", delta))
  }
})

## Permanent random numbers

test_that("scps with prn is a function of pik, spread and prn alone", {
  set.seed(109)
  N <- 40L
  pik <- inclusion_prob(rexp(N) + 0.2, 8)
  spread <- matrix(runif(2L * N), N)
  u <- runif(N)

  seed_before <- .Random.seed
  first <- balanced_wor(pik, spread = spread, method = "scps", prn = u)$sample
  expect_identical(.Random.seed, seed_before)
  second <- balanced_wor(pik, spread = spread, method = "scps", prn = u)$sample
  expect_identical(first, second)
  expect_length(first, 8L)

  # Gridded spread has many equal distances, a different quickselect path
  grid <- matrix(as.double(sample(0:2, 2L * N, TRUE)), N)
  expect_identical(
    balanced_wor(pik, spread = grid, method = "scps", prn = u)$sample,
    balanced_wor(pik, spread = grid, method = "scps", prn = u)$sample
  )
})

test_that("scps with prn visits units in row order", {
  # Two separated pairs with pik = 0.5. Unit 1 gives all its weight to
  # unit 2 and unit 3 to unit 4, so u_1 decides the first pair and u_3 the
  # second, and u_2, u_4 are never used.
  pik <- rep(0.5, 4)
  spread <- c(0, 1, 10, 11)
  for (u1 in c(0.1, 0.9)) {
    for (u3 in c(0.2, 0.8)) {
      for (rest in c(0.05, 0.95)) {
        u <- c(u1, rest, u3, 1 - rest)
        s <- balanced_wor(pik, spread = spread, method = "scps", prn = u)
        expected <- c(if (u1 < 0.5) 1L else 2L, if (u3 < 0.5) 3L else 4L)
        expect_identical(s$sample, expected)
      }
    }
  }
})

test_that("scps with prn respects first-order inclusion probabilities", {
  set.seed(110)
  pik <- c(1, 0, 0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9)
  spread <- cbind(c(0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4), c(0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0))
  nrep <- 8000L
  hits <- integer(length(pik))
  sizes <- integer(nrep)
  for (draw in seq_len(nrep)) {
    index <- balanced_wor(
      pik,
      spread = spread,
      method = "scps",
      prn = runif(length(pik))
    )$sample
    sizes[draw] <- length(index)
    hits[index] <- hits[index] + 1L
  }
  empirical <- hits / nrep

  expect_true(all(sizes == 6L))
  free <- pik > 0 & pik < 1
  z <- (empirical[free] - pik[free]) / sqrt(pik[free] * (1 - pik[free]) / nrep)

  expect_equal(empirical[!free], pik[!free])
  expect_lt(max(abs(z)), 4)
})

test_that("scps coordinates samples through prn", {
  set.seed(111)
  N <- 60L
  pik1 <- inclusion_prob(rexp(N) + 0.5, 12)
  pik2 <- inclusion_prob(rexp(N) + 0.5, 15)
  spread <- matrix(runif(2L * N), N)
  draw <- function(pik, u) {
    balanced_wor(pik, spread = spread, method = "scps", prn = u)$sample
  }
  overlap <- function(pair) length(intersect(pair[[1]], pair[[2]]))

  same <- replicate(200L, {
    u <- runif(N)
    overlap(list(draw(pik1, u), draw(pik2, u)))
  })
  opposite <- replicate(200L, {
    u <- runif(N)
    overlap(list(draw(pik1, u), draw(pik2, 1 - u)))
  })
  independent <- sum(pik1 * pik2)

  expect_gt(mean(same), independent + 1)
  expect_lt(mean(opposite), independent - 1)

  # Same probabilities and prn: the same sample
  u <- runif(N)
  expect_identical(draw(pik1, u), draw(pik1, u))
})

test_that("scps without prn keeps its random visiting order", {
  # Without prn the step unit is drawn at random, so row order does not
  # fix the sample. With prn it does.
  pik <- rep(0.5, 4)
  spread <- c(0, 1, 10, 11)
  set.seed(112)
  samples <- replicate(
    200L,
    paste(balanced_wor(pik, spread = spread, method = "scps")$sample, collapse = "-")
  )
  expect_setequal(unique(samples), c("1-3", "1-4", "2-3", "2-4"))
})

test_that("scps validates prn and rejects it where unsupported", {
  pik <- rep(0.5, 10)
  spread <- matrix(runif(20), ncol = 2)

  expect_error(
    balanced_wor(pik, spread = spread, method = "scps", prn = runif(9)),
    "must have length 10"
  )
  expect_error(
    balanced_wor(pik, spread = spread, method = "scps", prn = c(0, runif(9))),
    "open interval"
  )
  expect_error(
    balanced_wor(pik, spread = spread, method = "scps", prn = runif(10), nrep = 2),
    "prn and nrep > 1"
  )
  expect_error(
    balanced_wor(pik, spread = spread, method = "lpm2", prn = runif(10)),
    "method 'lpm2' does not support 'prn'.*'scps'"
  )
  expect_error(
    balanced_wor(pik, prn = runif(10)),
    "method 'cube' does not support 'prn'.*'scps'"
  )
})

test_that("scps shares tied weight equally up to each unit's bound", {
  # All units share one point, so every step's neighbours tie. With prn the
  # thresholds follow by hand, and a prn 1e-9 either side of one flips it.
  draw <- function(pik, u) {
    balanced_wor(pik, spread = rep(0, 4), method = "scps", prn = u)$sample
  }

  # Every tied unit takes an equal third of unit 1's weight, so unit 2 is
  # left at 1/3 once unit 1 is selected.
  pik <- rep(0.5, 4)
  expect_identical(draw(pik, c(0.25, 1 / 3 - 1e-9, 0.9, 0.9)), c(1L, 2L))
  expect_identical(draw(pik, c(0.25, 1 / 3 + 1e-9, 0.4, 0.9)), c(1L, 3L))

  # Unit 4 can take only 0.2, so units 2 and 3 take 0.4 each and unit 2 is
  # left at 0.1.
  pik <- c(0.5, 0.3, 0.3, 0.9)
  expect_identical(draw(pik, c(0.25, 0.1 - 1e-9, 0.5, 0.5)), c(1L, 2L))
  expect_identical(draw(pik, c(0.25, 0.1 + 1e-9, 0.5, 0.5)), c(1L, 4L))
})
