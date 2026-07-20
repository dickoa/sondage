# Local pivotal method 2 (Grafström, Lundström & Schelin, 2012)

test_that("lpm2 returns sondage_sample object", {
  set.seed(1)
  pik <- c(0.2, 0.4, 0.6, 0.8)
  z <- matrix(runif(8), ncol = 2)
  s <- balanced_wor(pik, spread = z, method = "lpm2")
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_s3_class(s, "balanced")
  expect_equal(s$method, "lpm2")
  expect_equal(s$n, 2L)
  expect_equal(s$N, 4L)
  expect_equal(s$pik, pik)
  expect_true(s$fixed_size)
})

test_that("lpm2 draws fixed-size samples of distinct in-range indices", {
  set.seed(2)
  N <- 50
  pik <- rep(0.2, N)
  z <- matrix(runif(2 * N), ncol = 2)
  for (i in 1:20) {
    idx <- balanced_wor(pik, spread = z, method = "lpm2")$sample
    expect_length(idx, 10)
    expect_true(all(idx >= 1 & idx <= N))
    expect_equal(anyDuplicated(idx), 0L)
    expect_false(is.unsorted(idx))
  }
})

test_that("lpm2 is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5, 0.6, 0.4)
  z <- matrix(seq_len(10) / 10, ncol = 2)

  set.seed(999)
  idx1 <- balanced_wor(pik, spread = z, method = "lpm2")$sample
  set.seed(999)
  idx2 <- balanced_wor(pik, spread = z, method = "lpm2")$sample

  expect_identical(idx1, idx2)
})

test_that("lpm2 respects first-order inclusion probabilities", {
  set.seed(42)
  pik <- c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9)
  z <- matrix(c(pik, runif(20)), ncol = 3)
  nrep <- 10000
  s <- balanced_wor(pik, spread = z, method = "lpm2", nrep = nrep)
  emp <- tabulate(as.integer(s$sample), nbins = 10) / nrep
  # 4 x binomial SE at p = 0.5, nrep = 10000
  expect_true(max(abs(emp - pik)) < 0.02)
})

test_that("lpm2 respects marginals on tie-heavy grid coordinates", {
  # Integer grid: exact distance ties everywhere; requires random
  # tie-breaking for E(s) = pik to hold.
  set.seed(43)
  N <- 25
  z <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  storage.mode(z) <- "double"
  pik <- rep(0.2, N)
  nrep <- 5000
  s <- balanced_wor(pik, spread = z, method = "lpm2", nrep = nrep)
  emp <- tabulate(as.integer(s$sample), nbins = N) / nrep
  expect_true(max(abs(emp - pik)) < 0.03)
})

test_that("lpm2 spreads: coincident pairs yield exactly one unit each", {
  # K far-apart locations, two units at each, pik = 0.5: while both
  # units of a pair are undecided each is the other's nearest
  # neighbour, so the first pivotal step touching a pair resolves it
  # to exactly one selected unit.
  set.seed(44)
  K <- 10
  centers <- cbind(10 * (1:K), 10 * (1:K))
  z <- centers[rep(1:K, each = 2), ]
  pik <- rep(0.5, 2 * K)
  for (i in 1:20) {
    idx <- balanced_wor(pik, spread = z, method = "lpm2")$sample
    pair_of <- (idx + 1L) %/% 2L
    expect_length(idx, K)
    expect_equal(sort(pair_of), 1:K)
  }
})

test_that("lpm2 is better spread than srs on random coordinates", {
  set.seed(45)
  N <- 100
  z <- matrix(runif(2 * N), ncol = 2)
  pik <- rep(0.2, N)
  nn_dist <- function(idx) {
    d <- as.matrix(dist(z[idx, ]))
    diag(d) <- Inf
    mean(apply(d, 1, min))
  }
  lpm <- mean(replicate(30, nn_dist(
    balanced_wor(pik, spread = z, method = "lpm2")$sample
  )))
  srs <- mean(replicate(30, nn_dist(sample.int(N, 20))))
  expect_gt(lpm, srs)
})

test_that("lpm2 passes exact 0 and 1 pik through", {
  set.seed(46)
  pik <- c(1, 0, rep(0.5, 8))
  z <- matrix(runif(20), ncol = 2)
  for (i in 1:10) {
    idx <- balanced_wor(pik, spread = z, method = "lpm2")$sample
    expect_true(1L %in% idx)
    expect_false(2L %in% idx)
    expect_length(idx, 5)
  }
})

test_that("lpm2 works when all units are pre-resolved", {
  z <- matrix(1:6 / 6, ncol = 2)
  s <- balanced_wor(c(1, 0, 1), spread = z, method = "lpm2")
  expect_identical(s$sample, c(1L, 3L))
})

test_that("lpm2 works with a single spread column and vector input", {
  set.seed(47)
  pik <- rep(0.5, 6)
  # .check_cube_aux coerces a vector to a one-column matrix
  s <- balanced_wor(pik, spread = as.double(1:6), method = "lpm2")
  expect_length(s$sample, 3)
})

test_that("lpm2 nrep returns matrix with valid fixed-size columns", {
  set.seed(48)
  N <- 20
  pik <- rep(0.25, N)
  z <- matrix(runif(2 * N), ncol = 2)
  s <- balanced_wor(pik, spread = z, method = "lpm2", nrep = 25)
  expect_true(is.matrix(s$sample))
  expect_equal(dim(s$sample), c(5L, 25L))
  for (j in 1:25) {
    col <- s$sample[, j]
    expect_true(all(col >= 1 & col <= N))
    expect_equal(anyDuplicated(col), 0L)
  }
})

test_that("lpm2 nrep is reproducible with set.seed", {
  pik <- rep(0.5, 10)
  z <- matrix(seq_len(20), ncol = 2)
  set.seed(50)
  m1 <- balanced_wor(pik, spread = z, method = "lpm2", nrep = 5)$sample
  set.seed(50)
  m2 <- balanced_wor(pik, spread = z, method = "lpm2", nrep = 5)$sample
  expect_identical(m1, m2)
})

test_that("lpm2 rejects aux, strata, and bounds", {
  pik <- rep(0.5, 10)
  z <- matrix(runif(20), ncol = 2)
  expect_error(
    balanced_wor(pik, aux = matrix(1, 10), spread = z, method = "lpm2"),
    "does not use auxiliary balancing variables"
  )
  expect_error(
    balanced_wor(pik, strata = rep(1:2, 5), spread = z, method = "lpm2"),
    "does not support 'strata'"
  )
  expect_error(
    balanced_wor(
      pik,
      spread = z,
      bounds = list(B = matrix(1, 10), lower = 0, upper = 5),
      method = "lpm2"
    ),
    "'bounds' is only supported by the built-in \"cube\" method"
  )
})

test_that("lpm2 requires spread", {
  expect_error(
    balanced_wor(rep(0.5, 10), method = "lpm2"),
    "method 'lpm2' requires 'spread'"
  )
})

test_that("cube spread error points to lpm2", {
  expect_error(
    balanced_wor(rep(0.5, 10), spread = matrix(runif(20), ncol = 2)),
    "use method = \"lpm2\", method = \"scps\", or a method registered"
  )
})

test_that("lpm2 validates spread", {
  pik <- rep(0.5, 10)
  expect_error(
    balanced_wor(pik, spread = matrix(runif(8), ncol = 2), method = "lpm2"),
    "nrow\\(spread\\) = 4 does not match length\\(pik\\) = 10"
  )
  z_na <- matrix(runif(20), ncol = 2)
  z_na[3, 1] <- NA
  expect_error(
    balanced_wor(pik, spread = z_na, method = "lpm2"),
    "spread must not contain NA, NaN, or Inf"
  )
  expect_error(
    balanced_wor(pik, spread = matrix(numeric(0), nrow = 10), method = "lpm2"),
    "'spread' must have at least one column"
  )
})

test_that("lpm2 rejects pik within eps of the boundaries", {
  z <- matrix(runif(8), ncol = 2)
  expect_error(
    balanced_wor(c(1e-12, 0.5, 0.5, 1 - 1e-12), spread = z, method = "lpm2"),
    "exactly 0 or exactly 1"
  )
  # and accepts them with a smaller eps
  set.seed(51)
  s <- balanced_wor(
    c(1e-12, 0.5, 0.5, 1 - 1e-12),
    spread = z,
    method = "lpm2",
    eps = 1e-14
  )
  expect_length(s$sample, 2)
})

test_that("lpm2 validates eps", {
  pik <- rep(0.5, 4)
  z <- matrix(runif(8), ncol = 2)
  expect_error(
    balanced_wor(pik, spread = z, method = "lpm2", eps = 0.7),
    "eps"
  )
  expect_error(
    balanced_wor(pik, spread = z, method = "lpm2", eps = c(0.1, 0.2)),
    "eps"
  )
})

test_that("print works for lpm2 objects", {
  set.seed(52)
  pik <- rep(0.5, 10)
  z <- matrix(runif(20), ncol = 2)
  s <- balanced_wor(pik, spread = z, method = "lpm2")
  expect_output(print(s), "Balanced WOR \\[lpm2\\]")
  expect_output(print(s), "n=5, N=10")
})

test_that("method_spec reports lpm2 capabilities", {
  spec <- method_spec("lpm2")
  expect_equal(spec$type, "balanced")
  expect_true(spec$fixed_size)
  expect_equal(spec$variance_family, "unsupported")
  expect_false(spec$supports_prn)
  expect_false(spec$supports_aux)
  expect_false(spec$supports_strata)
  expect_true(spec$supports_spread)
})

test_that("lpm2 is a built-in name that cannot be registered", {
  expect_error(
    register_method(
      "lpm2",
      type = "balanced",
      sample_fn = function(pik, n = NULL, aux = NULL, ...) seq_len(n)
    ),
    "'lpm2' is a built-in method and cannot be overridden"
  )
})

test_that("joint_inclusion_prob and sampling_cov error for lpm2", {
  set.seed(53)
  pik <- rep(0.5, 10)
  z <- matrix(runif(20), ncol = 2)
  s <- balanced_wor(pik, spread = z, method = "lpm2")
  expect_error(
    joint_inclusion_prob(s),
    "joint_inclusion_prob not implemented for method 'lpm2'"
  )
  expect_error(
    joint_inclusion_prob(s, sampled_only = TRUE),
    "joint_inclusion_prob not implemented for method 'lpm2'"
  )
  expect_error(sampling_cov(s))
})

test_that("inclusion_prob extracts pik from lpm2 designs", {
  set.seed(54)
  pik <- c(0.2, 0.4, 0.6, 0.8)
  z <- matrix(runif(8), ncol = 2)
  s <- balanced_wor(pik, spread = z, method = "lpm2")
  expect_equal(inclusion_prob(s), pik)
})

test_that("lpm2 handles non-uniform pik summing to 1", {
  set.seed(55)
  pik <- c(0.1, 0.2, 0.3, 0.4)
  z <- matrix(runif(8), ncol = 2)
  nrep <- 5000
  s <- balanced_wor(pik, spread = z, method = "lpm2", nrep = nrep)
  expect_equal(dim(s$sample), c(1L, nrep))
  emp <- tabulate(as.integer(s$sample), nbins = 4) / nrep
  expect_true(max(abs(emp - pik)) < 0.03)
})
