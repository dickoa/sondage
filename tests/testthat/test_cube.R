compute_ht_estimator <- function(x, prob, sample_idx) {
  colSums(x[sample_idx, , drop = FALSE] / prob[sample_idx])
}

compute_empirical_prob <- function(samples, N) {
  if (is.matrix(samples)) {
    tabulate(as.vector(samples), nbins = N) / ncol(samples)
  } else {
    counts <- integer(N)
    for (s in samples) {
      counts[s] <- counts[s] + 1L
    }
    counts / length(samples)
  }
}

test_that("cube returns integer vector of indices", {
  set.seed(42)
  N <- 100
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  s <- cube(prob, x)

  expect_type(s, "integer")
  expect_true(all(s >= 1 & s <= N))
  expect_true(length(s) == length(unique(s))) # No duplicates
})

test_that("cube respects fixed sample size when prob is balancing variable", {
  set.seed(123)
  N <- 100
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  for (i in 1:20) {
    s <- cube(prob, x)
    expect_equal(length(s), n)
  }
})

test_that("cube handles equal probabilities", {
  set.seed(456)
  N <- 50
  n <- 5
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  s <- cube(prob, x)

  expect_equal(length(s), n)
  expect_true(all(s >= 1 & s <= N))
})

test_that("cube handles unequal probabilities", {
  set.seed(789)
  N <- 100
  prob <- runif(N, 0.05, 0.3)
  prob <- prob / sum(prob) * 20 # Expected sample size ~20
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  s <- cube(prob, x)

  expect_true(length(s) > 0)
  expect_true(length(s) <= N)
  expect_true(all(s >= 1 & s <= N))
})

test_that("cube handles probability bounds correctly", {
  set.seed(101)
  N <- 50
  prob <- c(rep(0, 10), rep(0.2, 30), rep(1, 10)) # Mix of 0, middle, and 1
  x <- cbind(prob, runif(N))

  s <- cube(prob, x)

  # Units with prob=1 must be in sample
  expect_true(all(41:50 %in% s))
  # Units with prob=0 must not be in sample
  expect_true(!any(1:10 %in% s))
})

test_that("cube produces approximately balanced samples", {
  set.seed(202)
  N <- 500
  n <- 50
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 3), ncol = 3))

  s <- cube(prob, x)

  ht_est <- compute_ht_estimator(x, prob, s)
  pop_total <- colSums(x)

  rel_dev <- abs(ht_est - pop_total) / pop_total * 100

  expect_equal(ht_est[1], pop_total[1], tolerance = 1e-9)

  expect_true(
    all(rel_dev[-1] < 5),
    info = sprintf("Max relative deviation: %.2f%%", max(rel_dev[-1]))
  )
})

test_that("cube respects inclusion probabilities on average", {
  skip_on_cran()

  set.seed(404)
  N <- 50
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  n_reps <- 500
  samples <- cube(prob, x, nrep = n_reps)

  emp_prob <- compute_empirical_prob(samples, N)

  # Check that empirical probabilities are close to target
  se <- sqrt(prob[1] * (1 - prob[1]) / n_reps)
  max_dev <- max(abs(emp_prob - prob))

  expect_lt(max_dev, 4 * se)
})

test_that("cube respects unequal inclusion probabilities", {
  skip_on_cran()

  set.seed(505)
  N <- 30
  prob <- seq(0.1, 0.5, length.out = N)
  x <- cbind(prob, runif(N))

  n_reps <- 500
  samples <- cube(prob, x, nrep = n_reps)

  emp_prob <- compute_empirical_prob(samples, N)

  corr <- cor(prob, emp_prob)
  expect_gt(corr, 0.95)

  rmse <- sqrt(mean((emp_prob - prob)^2))
  expect_lt(rmse, 0.05)
})

test_that("cube handles single balancing variable", {
  set.seed(606)
  N <- 100
  n <- 10
  prob <- rep(n / N, N)
  x <- matrix(prob, ncol = 1)

  s <- cube(prob, x)

  expect_equal(length(s), n)
})

test_that("cube handles many balancing variables", {
  set.seed(707)
  N <- 200
  n <- 50
  p <- 15 # Many balancing variables
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * (p - 1)), ncol = p - 1))

  s <- cube(prob, x)

  expect_equal(length(s), n)
  expect_true(all(s >= 1 & s <= N))
})

test_that("cube handles small populations", {
  set.seed(808)
  N <- 10
  n <- 3
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  s <- cube(prob, x)

  expect_equal(length(s), n)
})

test_that("cube handles n close to N", {
  set.seed(909)
  N <- 20
  n <- 18
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  s <- cube(prob, x)

  expect_equal(length(s), n)
})

test_that("cube validates prob input", {
  N <- 10
  x <- matrix(runif(N * 2), ncol = 2)

  expect_error(cube(c(-0.1, rep(0.1, N - 1)), x))
  expect_error(cube(c(1.1, rep(0.1, N - 1)), x))
  expect_error(cube("not numeric", x))
})

test_that("cube validates x input", {
  N <- 10
  prob <- rep(0.2, N)

  expect_error(cube(prob, matrix(runif(5), ncol = 1))) # Wrong N
})

test_that("cube validates eps input", {
  N <- 10
  prob <- rep(0.2, N)
  x <- matrix(runif(N * 2), ncol = 2)

  expect_error(cube(prob, x, eps = -1))
  expect_error(cube(prob, x, eps = 0))
})

test_that("cube with nrep returns matrix of correct dimensions", {
  set.seed(111)
  N <- 50
  n <- 5
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  nrep <- 10
  samples <- cube(prob, x, nrep = nrep)

  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), n)
  expect_equal(ncol(samples), nrep)
})

test_that("cube with nrep produces valid samples", {
  set.seed(222)
  N <- 50
  n <- 5
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  nrep <- 20
  samples <- cube(prob, x, nrep = nrep)

  for (i in seq_len(nrep)) {
    s <- samples[, i]
    expect_type(s, "integer")
    expect_equal(length(s), n)
    expect_true(all(s >= 1 & s <= N))
    expect_true(length(s) == length(unique(s)))
  }
})

test_that("cube with nrep respects inclusion probabilities", {
  skip_on_cran()

  set.seed(333)
  N <- 30
  n <- 6
  prob <- rep(n / N, N)
  x <- cbind(prob, runif(N))

  nrep <- 500
  samples <- cube(prob, x, nrep = nrep)

  emp_prob <- compute_empirical_prob(samples, N)

  se <- sqrt(prob[1] * (1 - prob[1]) / nrep)
  max_dev <- max(abs(emp_prob - prob))

  expect_lt(max_dev, 4 * se)
})

test_that("cube validates nrep input", {
  N <- 10
  prob <- rep(0.2, N)
  x <- cbind(prob, runif(N))

  expect_error(cube(prob, x, nrep = 0))
  expect_error(cube(prob, x, nrep = -1))
})

test_that("stratified cube returns correct indices", {
  set.seed(1001)
  N <- 200
  strata <- rep(1:4, each = 50)
  prob <- rep(10 / 50, N)
  x <- cbind(prob, runif(N))

  s <- cube(prob, x, strata = strata)

  expect_type(s, "integer")
  expect_true(all(s >= 1 & s <= N))
  expect_true(length(s) == length(unique(s)))
})

test_that("stratified cube guarantees exact sample size per stratum", {
  set.seed(1002)
  N <- 200
  strata <- rep(1:4, each = 50)
  n_per_stratum <- 10
  prob <- rep(n_per_stratum / 50, N)
  x <- cbind(prob, runif(N))

  for (i in 1:100) {
    s <- cube(prob, x, strata = strata)
    strata_counts <- table(strata[s])
    expect_true(all(strata_counts == n_per_stratum))
  }
})

test_that("stratified cube works with unequal stratum sizes", {
  set.seed(1003)
  strata <- c(rep(1, 100), rep(2, 50), rep(3, 30), rep(4, 20))
  N <- length(strata)
  # Equal sampling fraction per stratum
  prob <- ifelse(
    strata == 1,
    10 / 100,
    ifelse(strata == 2, 5 / 50, ifelse(strata == 3, 3 / 30, 2 / 20))
  )
  x <- cbind(prob, runif(N))

  s <- cube(prob, x, strata = strata)
  strata_counts <- table(strata[s])

  expect_equal(as.numeric(strata_counts["1"]), 10)
  expect_equal(as.numeric(strata_counts["2"]), 5)
  expect_equal(as.numeric(strata_counts["3"]), 3)
  expect_equal(as.numeric(strata_counts["4"]), 2)
})

test_that("stratified cube with character strata", {
  set.seed(1004)
  N <- 120
  strata <- rep(c("North", "South", "East", "West"), each = 30)
  prob <- rep(5 / 30, N)
  x <- cbind(prob, runif(N))

  s <- cube(prob, x, strata = strata)
  strata_counts <- table(strata[s])

  expect_true(all(strata_counts == 5))
})

test_that("stratified cube balances on auxiliary variables", {
  set.seed(1005)
  N <- 200
  strata <- rep(1:4, each = 50)
  prob <- rep(10 / 50, N)
  x_aux <- runif(N)
  x <- cbind(prob, x_aux)

  s <- cube(prob, x, strata = strata)

  # Check HT estimator is close to population total
  ht_est <- sum(x_aux[s] / prob[s])
  pop_total <- sum(x_aux)

  rel_dev <- abs(ht_est - pop_total) / pop_total
  expect_lt(rel_dev, 0.05) # < 5% deviation
})

test_that("stratified cube respects inclusion probabilities", {
  skip_on_cran()

  set.seed(1006)
  N <- 200
  strata <- rep(1:4, each = 50)
  prob <- rep(10 / 50, N)
  x <- cbind(prob, runif(N))

  nrep <- 500
  samples <- cube(prob, x, strata = strata, nrep = nrep)

  emp_prob <- compute_empirical_prob(samples, N)

  se <- sqrt(prob[1] * (1 - prob[1]) / nrep)
  max_dev <- max(abs(emp_prob - prob))

  expect_lt(max_dev, 4 * se)
})

test_that("stratified cube with nrep returns matrix", {
  set.seed(1007)
  N <- 200
  strata <- rep(1:4, each = 50)
  prob <- rep(10 / 50, N)
  x <- cbind(prob, runif(N))

  nrep <- 10
  samples <- cube(prob, x, strata = strata, nrep = nrep)

  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), 40) # 10 per stratum * 4 strata
  expect_equal(ncol(samples), nrep)
})

test_that("stratified cube validates strata input", {
  N <- 10
  prob <- rep(0.2, N)
  x <- cbind(prob, runif(N))

  expect_error(cube(prob, x, strata = 1:5))
})

test_that("cube is reproducible with set.seed", {
  N <- 50
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  set.seed(12345)
  s1 <- cube(prob, x)

  set.seed(12345)
  s2 <- cube(prob, x)

  expect_equal(s1, s2)
})

test_that("cube with nrep is reproducible with set.seed", {
  N <- 50
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  set.seed(54321)
  samples1 <- cube(prob, x, nrep = 5)

  set.seed(54321)
  samples2 <- cube(prob, x, nrep = 5)

  expect_equal(samples1, samples2)
})

test_that("stratified cube is reproducible with set.seed", {
  N <- 200
  strata <- rep(1:4, each = 50)
  prob <- rep(10 / 50, N)
  x <- cbind(prob, runif(N))

  set.seed(99999)
  s1 <- cube(prob, x, strata = strata)

  set.seed(99999)
  s2 <- cube(prob, x, strata = strata)

  expect_equal(s1, s2)
})

test_that("cube produces similar distribution to BalancedSampling", {
  skip_if_not_installed("BalancedSampling")
  skip_on_cran()

  bs_cube <- getFromNamespace("cube", "BalancedSampling")

  set.seed(999)
  N <- 50
  n <- 10
  prob <- rep(n / N, N)
  x <- cbind(prob, matrix(runif(N * 2), ncol = 2))

  n_reps <- 300

  samples_ours <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    samples_ours[[i]] <- cube(prob, x)
  }
  emp_ours <- compute_empirical_prob(samples_ours, N)

  samples_bs <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    samples_bs[[i]] <- bs_cube(prob, x)
  }
  emp_bs <- compute_empirical_prob(samples_bs, N)

  expect_lt(max(abs(emp_ours - prob)), 0.085)
  expect_lt(max(abs(emp_bs - prob)), 0.085)

  expect_lt(max(abs(emp_ours - emp_bs)), 0.085)
})
