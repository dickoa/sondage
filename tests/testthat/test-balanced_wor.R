test_that("balanced_wor returns sondage_sample object", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))
  s <- balanced_wor(pik, aux = x)
  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "cube")
  expect_equal(s$n, 2L)
  expect_equal(s$N, 4L)
  expect_equal(s$pik, pik)
  expect_true(s$fixed_size)
})

test_that("balanced_wor returns correct number of indices", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))
  expect_length(balanced_wor(pik, aux = x)$sample, 2)
})

test_that("balanced_wor indices are in valid range", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))
  idx <- balanced_wor(pik, aux = x)$sample
  expect_true(all(idx >= 1 & idx <= 4))
})

test_that("balanced_wor has no duplicates", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))
  idx <- balanced_wor(pik, aux = x)$sample
  expect_equal(length(unique(idx)), length(idx))
})

test_that("balanced_wor is reproducible with set.seed", {
  pik <- c(0.2, 0.3, 0.5)
  x <- matrix(c(10, 20, 30))

  set.seed(999)
  idx1 <- balanced_wor(pik, aux = x)$sample
  set.seed(999)
  idx2 <- balanced_wor(pik, aux = x)$sample

  expect_identical(idx1, idx2)
})

test_that("print works for balanced_wor objects", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30, 40)))
  expect_output(print(s), "Unequal prob WOR \\[cube\\]")
  expect_output(print(s), "n=2, N=4")
})


test_that("balanced_wor works with aux = NULL", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik)
  expect_length(s$sample, 2)
  expect_true(all(s$sample >= 1 & s$sample <= 4))
})


test_that("balanced_wor works with multiple aux columns", {
  set.seed(42)
  N <- 30
  pik <- inclusion_prob(runif(N, 1, 10), n = 8)
  x <- matrix(rnorm(N * 3), ncol = 3)

  s <- balanced_wor(pik, aux = x)
  expect_length(s$sample, 8)
  expect_true(all(s$sample >= 1 & s$sample <= N))
})

test_that("balanced_wor accepts vector aux (treated as 1-column)", {
  pik <- c(0.3, 0.3, 0.4)
  x <- c(10, 20, 30)
  s <- balanced_wor(pik, aux = x)
  expect_length(s$sample, 1)
})


test_that("cube achieves correct inclusion probabilities", {
  skip_on_cran()
  pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  N <- length(pik)
  x <- matrix(rnorm(N))

  set.seed(42)
  sim <- balanced_wor(pik, aux = x, nrep = 5000)
  counts <- integer(N)
  for (j in seq_len(ncol(sim$sample))) {
    counts[sim$sample[, j]] <- counts[sim$sample[, j]] + 1L
  }
  pi_hat <- counts / ncol(sim$sample)

  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("cube always produces exact sample size n", {
  skip_on_cran()
  pik <- c(0.2, 0.4, 0.6, 0.8)
  n <- round(sum(pik))
  x <- matrix(c(10, 20, 30, 40))
  nrep <- 500

  set.seed(42)
  sim <- balanced_wor(pik, aux = x, nrep = nrep)
  expect_true(is.matrix(sim$sample))
  expect_equal(nrow(sim$sample), n)
  expect_equal(ncol(sim$sample), nrep)
})


test_that("cube reduces HT estimator variance on balanced variable", {
  skip_on_cran()
  set.seed(1)
  N <- 50
  pik <- inclusion_prob(runif(N, 1, 10), n = 10)
  x <- rnorm(N, mean = 100, sd = 30)
  nrep <- 3000

  sim_cube <- balanced_wor(pik, aux = matrix(x), nrep = nrep)
  ht_cube <- apply(sim_cube$sample, 2, function(s) sum(x[s] / pik[s]))

  sim_cps <- unequal_prob_wor(pik, method = "cps", nrep = nrep)
  ht_cps <- apply(sim_cps$sample, 2, function(s) sum(x[s] / pik[s]))

  # Cube should have substantially lower variance for balanced variable
  expect_true(var(ht_cube) < 0.5 * var(ht_cps))
})

test_that("HT estimator is unbiased for cube", {
  skip_on_cran()
  pik <- c(0.2, 0.4, 0.6, 0.8)
  y <- c(10, 25, 15, 50)
  Y <- sum(y)
  x <- matrix(rnorm(4))

  set.seed(42)
  sim <- balanced_wor(pik, aux = x, nrep = 10000)
  ht_estimates <- apply(sim$sample, 2, function(s) sum(y[s] / pik[s]))
  expect_equal(mean(ht_estimates), Y, tolerance = 0.05 * Y)
})

test_that("stratified cube preserves per-stratum sample sizes", {
  skip_on_cran()
  set.seed(42)
  N <- 40
  strata <- rep(1:4, each = 10)
  pik <- rep(0.3, N)
  x <- matrix(as.double(seq_len(N)), ncol = 1)

  nrep <- 200
  for (i in seq_len(nrep)) {
    s <- balanced_wor(pik, aux = x, strata = strata)
    tab <- tabulate(strata[s$sample], nbins = 4)
    expect_equal(tab, rep(3L, 4), label = paste("rep", i))
  }
})

test_that("stratified cube returns correct class and fields", {
  pik <- rep(0.4, 20)
  x <- matrix(as.double(1:20), ncol = 1)
  strata <- rep(1:4, each = 5)
  s <- balanced_wor(pik, aux = x, strata = strata)
  expect_s3_class(s, "wor")
  expect_equal(s$method, "cube")
  expect_equal(s$N, 20L)
  expect_equal(s$n, 8L)
})

test_that("stratified cube works with unequal stratum sizes", {
  # 3 strata: sizes 5, 10, 5
  strata <- c(rep(1L, 5), rep(2L, 10), rep(3L, 5))
  N <- length(strata)
  pik <- rep(0.4, N)
  x <- matrix(as.double(seq_len(N)), ncol = 1)

  set.seed(42)
  s <- balanced_wor(pik, aux = x, strata = strata)
  tab <- tabulate(strata[s$sample], nbins = 3)
  expect_equal(tab, c(2L, 4L, 2L))
})

test_that("stratified cube achieves correct inclusion probs", {
  skip_on_cran()
  N <- 20
  pik <- rep(0.4, N)
  x <- matrix(as.double(1:N), ncol = 1)
  strata <- rep(1:4, each = 5)

  set.seed(42)
  sim <- balanced_wor(pik, aux = x, strata = strata, nrep = 5000)
  counts <- integer(N)
  for (j in seq_len(ncol(sim$sample))) {
    counts[sim$sample[, j]] <- counts[sim$sample[, j]] + 1L
  }
  pi_hat <- counts / ncol(sim$sample)
  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("inclusion_prob works on cube design", {
  pik <- c(0.2, 0.3, 0.5)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30)))
  expect_equal(inclusion_prob(s), pik)
})

test_that("joint_inclusion_prob works on cube design (HE approx)", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30, 40)))
  pikl <- joint_inclusion_prob(s)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(4, 4))
  # Diagonal = pik
  expect_equal(diag(pikl), pik)
  # Symmetric
  expect_equal(pikl, t(pikl))
  # All entries non-negative
  expect_true(all(pikl >= -1e-10))
  # All joint probs positive for n >= 2
  expect_true(all(pikl > 0))
})

test_that("sampling_cov works on cube design", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30, 40)))
  delta <- sampling_cov(s)

  expect_true(is.matrix(delta))
  expect_equal(dim(delta), c(4, 4))
  # Symmetric
  expect_equal(delta, t(delta))
  # Off-diagonal should be negative for well-behaved WOR
  off_diag <- delta[row(delta) != col(delta)]
  expect_true(all(off_diag < 0))
})

test_that("sampling_cov weighted works on cube design", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30, 40)))
  chk <- sampling_cov(s, weighted = TRUE)

  expect_true(is.matrix(chk))
  # Diagonal: 1 - pik
  expect_equal(diag(chk), 1 - pik, tolerance = 1e-10)
  # Off-diagonal should be non-positive (SYG condition)
  off_diag <- chk[row(chk) != col(chk)]
  expect_true(all(off_diag <= 1e-10))
})

test_that("cube batch returns correct structure", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))
  nrep <- 50

  set.seed(42)
  s <- balanced_wor(pik, aux = x, nrep = nrep)
  expect_s3_class(s, "sondage_sample")
  expect_true(is.matrix(s$sample))
  expect_equal(nrow(s$sample), 2)
  expect_equal(ncol(s$sample), nrep)
  expect_equal(s$pik, pik)
  expect_equal(s$method, "cube")
})

test_that("cube batch is reproducible with set.seed", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))

  set.seed(4321)
  s1 <- balanced_wor(pik, aux = x, nrep = 60)$sample
  set.seed(4321)
  s2 <- balanced_wor(pik, aux = x, nrep = 60)$sample

  expect_identical(s1, s2)
})

test_that("cube batch prints correctly", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- balanced_wor(pik, aux = matrix(c(10, 20, 30, 40)), nrep = 10)
  expect_output(print(s), "10 replicates")
})

test_that("stratified cube batch preserves per-stratum sizes", {
  skip_on_cran()
  N <- 20
  pik <- rep(0.4, N)
  x <- matrix(as.double(1:N), ncol = 1)
  strata <- rep(1:4, each = 5)

  set.seed(42)
  sim <- balanced_wor(pik, aux = x, strata = strata, nrep = 100)
  for (j in seq_len(ncol(sim$sample))) {
    tab <- tabulate(strata[sim$sample[, j]], nbins = 4)
    expect_equal(tab, rep(2L, 4), label = paste("batch rep", j))
  }
})


test_that("stratified cube batch is reproducible when stratum sizes are integer", {
  pik <- rep(0.4, 20)
  x <- matrix(as.double(1:20), ncol = 1)
  strata <- rep(1:4, each = 5)

  set.seed(987)
  s1 <- balanced_wor(pik, aux = x, strata = strata, nrep = 40)$sample
  set.seed(987)
  s2 <- balanced_wor(pik, aux = x, strata = strata, nrep = 40)$sample

  expect_identical(s1, s2)
})

test_that("cube handles certainty units (pik = 1)", {
  pik <- c(1.0, 0.5, 0.5)
  x <- matrix(c(10, 20, 30))
  set.seed(42)
  for (i in 1:50) {
    s <- balanced_wor(pik, aux = x)
    expect_true(1L %in% s$sample)
  }
})

test_that("cube handles near-zero pik", {
  pik <- c(1e-8, 0.5, 0.5)
  x <- matrix(c(10, 20, 30))
  set.seed(42)
  for (i in 1:50) {
    s <- balanced_wor(pik, aux = x)
    expect_false(1L %in% s$sample)
  }
})

test_that("cube handles N=1 census", {
  s <- balanced_wor(1.0)
  expect_equal(s$sample, 1L)
  expect_equal(s$n, 1L)
})

test_that("cube handles all-certainty pik", {
  pik <- c(1.0, 1.0, 1.0)
  s <- balanced_wor(pik)
  expect_equal(sort(s$sample), 1:3)
})

test_that("cube works with equal pik", {
  pik <- rep(0.3, 10)
  x <- matrix(as.double(1:10))
  s <- balanced_wor(pik, aux = x)
  expect_length(s$sample, 3)
})

test_that("cube works with large N", {
  set.seed(123)
  N <- 1000
  n <- 100
  pik <- inclusion_prob(runif(N, 1, 10), n = n)
  x <- matrix(rnorm(N))
  s <- balanced_wor(pik, aux = x)
  expect_length(s$sample, n)
  expect_true(all(s$sample >= 1 & s$sample <= N))
})

test_that("balanced_wor rejects invalid pik", {
  expect_error(balanced_wor(c(0.5, NA, 0.5)), "missing values")
  expect_error(balanced_wor(c("a", "b")), "numeric")
  expect_error(balanced_wor(numeric(0)), "empty")
})

test_that("balanced_wor rejects non-integer sum(pik)", {
  expect_error(balanced_wor(c(0.49, 0.49, 0.49)), "not close to an integer")
})

test_that("balanced_wor rejects invalid nrep", {
  expect_error(balanced_wor(c(0.5, 0.5), nrep = 0), "at least 1")
  expect_error(balanced_wor(c(0.5, 0.5), nrep = Inf), "finite")
  expect_error(balanced_wor(c(0.5, 0.5), nrep = 2.9), "not close to an integer")
})

test_that("balanced_wor rejects mismatched aux dimensions", {
  pik <- c(0.5, 0.5)
  expect_error(balanced_wor(pik, aux = matrix(1:6, ncol = 2)), "does not match")
})

test_that("balanced_wor rejects invalid strata", {
  pik <- rep(0.5, 4)
  x <- matrix(as.double(1:4))
  expect_error(balanced_wor(pik, aux = x, strata = 1:3), "does not match")
  expect_error(
    balanced_wor(pik, aux = x, strata = c(1, 2, NA, 1)),
    "missing values"
  )
  expect_error(
    balanced_wor(pik, aux = x, strata = c(1, 0, 2, 1)),
    "positive integers"
  )
})

test_that("balanced_wor rejects invalid method", {
  expect_error(balanced_wor(c(0.5, 0.5), method = "foo"))
})

test_that("cube handles degenerate/collinear auxiliary variables", {
  N <- 50
  n <- 10
  pik <- rep(n / N, N)

  # Identical columns (rank-1 auxiliary matrix)
  set.seed(42)
  aux_dup <- cbind(1:N, 1:N)
  s <- balanced_wor(pik, aux = aux_dup, method = "cube")
  expect_length(s$sample, n)
  expect_true(all(s$sample %in% 1:N))
  expect_false(anyDuplicated(s$sample) > 0)

  # Constant column (zero variance after A = X/pi)
  aux_const <- cbind(1:N, rep(1, N))
  s2 <- balanced_wor(pik, aux = aux_const, method = "cube")
  expect_length(s2$sample, n)
  expect_true(all(s2$sample %in% 1:N))
  expect_false(anyDuplicated(s2$sample) > 0)

  # All-identical rows
  aux_identical <- matrix(1, nrow = N, ncol = 3)
  s3 <- balanced_wor(pik, aux = aux_identical, method = "cube")
  expect_length(s3$sample, n)
  expect_true(all(s3$sample %in% 1:N))
  expect_false(anyDuplicated(s3$sample) > 0)
})

test_that("stratified cube warns and sets fixed_size=FALSE for non-integer stratum sums", {
  strata <- rep(1:3, each = 4)
  pik <- c(rep(0.35, 4), rep(0.40, 4), rep(0.50, 4))
  # per-stratum sums: 1.4, 1.6, 2.0
  expect_warning(
    balanced_wor(pik, strata = strata),
    "per-stratum sum\\(pik\\) is not close to an integer"
  )
  s <- suppressWarnings(balanced_wor(pik, strata = strata))
  expect_false(s$fixed_size)
})

test_that("stratified cube batch returns list for non-integer stratum sums", {
  strata <- rep(1:3, each = 4)
  pik <- c(rep(0.35, 4), rep(0.40, 4), rep(0.50, 4))
  nrep <- 20
  suppressWarnings(
    s <- balanced_wor(pik, strata = strata, nrep = nrep)
  )
  expect_false(s$fixed_size)
  expect_true(is.list(s$sample))
  expect_length(s$sample, nrep)
  # Each replicate is an integer vector of valid indices
  for (i in seq_len(nrep)) {
    expect_true(all(s$sample[[i]] >= 1L & s$sample[[i]] <= length(pik)))
  }
})

test_that("stratified cube does not warn when per-stratum sums are integer", {
  strata <- rep(1:3, each = 4)
  pik <- rep(0.50, 12)
  # per-stratum sums: 2, 2, 2
  expect_no_warning(balanced_wor(pik, strata = strata))
})

test_that("stratified cube handles sparse stratum labels", {
  set.seed(42)
  N <- 80
  pik <- rep(0.40, N)
  aux <- matrix(rnorm(N), ncol = 1)
  strata <- c(rep(1L, 40), rep(1000000L, 40))

  s <- balanced_wor(pik, aux = aux, strata = strata)
  expect_equal(s$n, 32L)
  expect_length(s$sample, 32)

  # Both strata should be represented
  tab <- table(strata[s$sample])
  expect_equal(length(tab), 2L)
})


# --- condition_aux tests ---

test_that("cube handles collinear aux with condition_aux = TRUE", {
  set.seed(42)
  N <- 60
  n <- 12
  pik <- rep(n / N, N)

  # Build highly collinear aux: 3 independent + 7 dependent
  base <- matrix(rnorm(N * 3), ncol = 3)
  dependent <- base %*% matrix(rnorm(21), nrow = 3, ncol = 7) +
    matrix(rnorm(N * 7, sd = 1e-8), ncol = 7)
  aux <- cbind(base, dependent)

  s <- balanced_wor(pik, aux = aux, condition_aux = TRUE)
  expect_length(s$sample, n)
  expect_true(all(s$sample >= 1 & s$sample <= N))
  expect_false(anyDuplicated(s$sample) > 0)
})

test_that("cube conditioning can be explicitly disabled", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  x <- matrix(c(10, 20, 30, 40))

  set.seed(1)
  s1 <- balanced_wor(pik, aux = x, condition_aux = FALSE)
  set.seed(1)
  s2 <- balanced_wor(pik, aux = x)

  # Default is FALSE, so both should give same result
  expect_identical(s1$sample, s2$sample)
})

test_that("condition_aux preserves correct inclusion probabilities", {
  skip_on_cran()
  set.seed(1)
  N <- 40
  pik <- inclusion_prob(runif(N, 1, 10), n = 10)

  # Collinear aux
  base <- matrix(rnorm(N * 2), ncol = 2)
  aux <- cbind(base, base[, 1] + base[, 2], 2 * base[, 1])

  sim <- balanced_wor(pik, aux = aux, condition_aux = TRUE, nrep = 5000)
  counts <- integer(N)
  for (j in seq_len(ncol(sim$sample))) {
    counts[sim$sample[, j]] <- counts[sim$sample[, j]] + 1L
  }
  pi_hat <- counts / ncol(sim$sample)
  expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("condition_aux works with stratified cube", {
  set.seed(42)
  N <- 40
  pik <- rep(0.4, N)
  strata <- rep(1:4, each = 10)
  base <- matrix(rnorm(N * 2), ncol = 2)
  aux <- cbind(base, base[, 1] + base[, 2])

  s <- balanced_wor(pik, aux = aux, strata = strata, condition_aux = TRUE)
  expect_length(s$sample, 16)
  tab <- tabulate(strata[s$sample], nbins = 4)
  expect_equal(tab, rep(4L, 4))
})

test_that("condition_aux with all-constant aux produces valid sample", {
  pik <- rep(0.4, 10)
  aux <- matrix(1, nrow = 10, ncol = 3)
  s <- balanced_wor(pik, aux = aux, condition_aux = TRUE)
  expect_length(s$sample, 4)
  expect_true(all(s$sample >= 1 & s$sample <= 10))
})
