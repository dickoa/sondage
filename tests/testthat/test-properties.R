pik4 <- c(0.2, 0.4, 0.6, 0.8) # n=2, N=4
pik5 <- inclusion_prob(c(10, 20, 30, 40, 50), n = 3) # n=3, N=5

# WR expected hits
hits5 <- expected_hits(c(10, 20, 30, 40, 50), n = 3)

# WOR fixed-size methods (excluding poisson: random size)
wor_methods <- c("cps", "brewer", "systematic")

# Methods with exact joint inclusion probabilities
# (Brewer uses an approximation, so exclude from exact analytical tests)
exact_jip_methods <- c("cps", "systematic")

nrep <- 10000 # simulation reps


# Property 1: HT estimator unbiasedness
test_that("HT estimator is unbiased for fixed-size WOR methods", {
  skip_on_cran()
  y <- c(10, 25, 15, 50)
  Y <- sum(y)

  for (method in wor_methods) {
    set.seed(42)
    sim <- unequal_prob_wor(pik4, method = method, nrep = nrep)
    ht_estimates <- apply(sim$sample, 2, function(s) sum(y[s] / pik4[s]))
    expect_equal(
      mean(ht_estimates),
      Y,
      tolerance = 0.05 * Y,
      label = paste("HT unbiased for", method)
    )
  }
})

test_that("HT estimator is unbiased for Poisson sampling", {
  skip_on_cran()
  y <- c(10, 25, 15, 50)
  Y <- sum(y)

  set.seed(42)
  sim <- unequal_prob_wor(pik4, method = "poisson", nrep = nrep)
  ht_estimates <- vapply(
    sim$sample,
    function(s) {
      if (length(s) == 0) {
        return(0)
      }
      sum(y[s] / pik4[s])
    },
    double(1)
  )
  expect_equal(mean(ht_estimates), Y, tolerance = 0.05 * Y)
})

test_that("HT estimator is unbiased for equal probability WOR", {
  skip_on_cran()
  N <- 20L
  n <- 5L
  y <- seq_len(N)
  Y <- sum(y)

  set.seed(42)
  sim <- equal_prob_wor(N, n, method = "srs", nrep = nrep)
  ht_estimates <- apply(sim$sample, 2, function(s) sum(y[s]) * N / n)
  expect_equal(mean(ht_estimates), Y, tolerance = 0.05 * Y)
})


# Property 2: Generalized HT unbiasedness for WR designs
test_that("Generalized HT estimator is unbiased for WR methods", {
  skip_on_cran()
  y <- c(10, 25, 15, 50, 35)
  Y <- sum(y)
  n <- round(sum(hits5))
  prob <- hits5 / sum(hits5)

  for (method in c("chromy", "multinomial")) {
    set.seed(42)
    sim <- unequal_prob_wr(hits5, method = method, nrep = nrep)
    ht_estimates <- apply(sim$sample, 2, function(s) {
      sum(y[s] / (n * prob[s]))
    })
    expect_equal(
      mean(ht_estimates),
      Y,
      tolerance = 0.05 * Y,
      label = paste("Generalized HT unbiased for", method)
    )
  }
})


# Property 3: Delta matrix row/column sums = 0 for fixed-size designs
test_that("Delta matrix rows sum to zero for fixed-size WOR", {
  for (method in exact_jip_methods) {
    s <- unequal_prob_wor(pik4, method = method)
    delta <- sampling_cov(s)
    row_sums <- rowSums(delta)
    expect_equal(
      row_sums,
      rep(0, length(pik4)),
      tolerance = 1e-6,
      label = paste("Delta row sums for", method)
    )
  }
})

test_that("Delta matrix rows sum to zero for SRS", {
  s <- equal_prob_wor(20, 5)
  delta <- sampling_cov(s)
  expect_equal(rowSums(delta), rep(0, 20), tolerance = 1e-10)
})


# Property 4: Negative covariance for WOR designs
test_that("Off-diagonal covariance is negative for WOR designs", {
  for (method in c("cps", "brewer")) {
    s <- unequal_prob_wor(pik4, method = method)
    delta <- sampling_cov(s)
    off_diag <- delta[row(delta) != col(delta)]
    expect_true(
      all(off_diag < 0),
      label = paste("Negative covariance for", method)
    )
  }
})


# Property 5: SYG check quantities non-positive for WOR
test_that("SYG check quantities are non-positive off-diagonal", {
  for (method in c("cps")) {
    s <- unequal_prob_wor(pik5, method = method)
    chk <- sampling_cov(s, weighted = TRUE)
    off_diag <- chk[row(chk) != col(chk)]
    expect_true(
      all(off_diag <= 1e-10),
      label = paste("SYG non-positive for", method)
    )
  }
})


# Property 6: Joint probability positivity (CPS, Brewer)
test_that("All joint probabilities positive for CPS and Brewer (n >= 2)", {
  for (method in c("cps", "brewer")) {
    s <- unequal_prob_wor(pik4, method = method)
    pikl <- joint_inclusion_prob(s)
    expect_true(all(pikl > 0), label = paste("All pikl > 0 for", method))
  }
})


# Property 7: Variance of sample size
test_that("Sample size is constant for fixed-size methods", {
  skip_on_cran()
  set.seed(42)
  sim_cps <- unequal_prob_wor(pik4, method = "cps", nrep = 500)
  expect_true(nrow(sim_cps$sample) == round(sum(pik4)))

  sim_brewer <- unequal_prob_wor(pik4, method = "brewer", nrep = 500)
  expect_true(nrow(sim_brewer$sample) == round(sum(pik4)))

  sim_chromy <- unequal_prob_wr(hits5, method = "chromy", nrep = 500)
  expect_true(nrow(sim_chromy$sample) == round(sum(hits5)))
})

test_that("Poisson has variable sample size", {
  skip_on_cran()
  set.seed(42)
  sim <- unequal_prob_wor(pik4, method = "poisson", nrep = 1000)
  sizes <- lengths(sim$sample)
  expect_true(var(sizes) > 0)
  expect_equal(mean(sizes), sum(pik4), tolerance = 0.1)
})

test_that("Bernoulli sample size variance matches N*p*(1-p)", {
  skip_on_cran()
  N <- 100L
  n <- 30L
  p <- n / N

  set.seed(42)
  sim <- equal_prob_wor(N, n, method = "bernoulli", nrep = 5000)
  sizes <- lengths(sim$sample)
  expect_equal(var(sizes), N * p * (1 - p), tolerance = 2)
  expect_equal(mean(sizes), n, tolerance = 1)
})


# Property 8: Poisson conditional independence
test_that("Poisson has independent selections (pi_ij = pi_i * pi_j)", {
  s <- unequal_prob_wor(pik4, method = "poisson")
  pikl <- joint_inclusion_prob(s)
  expected <- outer(pik4, pik4)
  diag(expected) <- pik4
  expect_equal(pikl, expected, tolerance = 1e-14)
})


# Property 9: Marginal constraints on joint inclusion probabilities
test_that("Joint prob marginal constraint for fixed-size WOR", {
  n <- round(sum(pik4))
  for (method in exact_jip_methods) {
    s <- unequal_prob_wor(pik4, method = method)
    pikl <- joint_inclusion_prob(s)
    for (i in seq_along(pik4)) {
      row_sum_off_diag <- sum(pikl[i, -i])
      expect_equal(
        row_sum_off_diag,
        (n - 1) * pik4[i],
        tolerance = 1e-5,
        label = paste(method, "marginal constraint, unit", i)
      )
    }
  }
})


# Property 10: WR covariance matrix properties
test_that("Multinomial covariance matches known formula", {
  s <- unequal_prob_wr(hits5, method = "multinomial")
  cov_mat <- sampling_cov(s)
  prob <- hits5 / sum(hits5)
  n <- round(sum(hits5))
  ehits <- n * prob

  # Off-diagonal: Cov(n_i, n_j) = -n * p_i * p_j
  for (i in seq_along(prob)) {
    for (j in seq_along(prob)) {
      if (i != j) {
        expect_equal(cov_mat[i, j], -n * prob[i] * prob[j], tolerance = 1e-10)
      }
    }
  }

  # Diagonal: Var(n_i) = n * p_i * (1 - p_i)
  for (i in seq_along(prob)) {
    expect_equal(cov_mat[i, i], n * prob[i] * (1 - prob[i]), tolerance = 1e-10)
  }
})


# Property 11: SRS known variance formula
test_that("SRS HT variance matches formula", {
  skip_on_cran()
  N <- 20L
  n <- 5L
  y <- seq_len(N)
  S2 <- var(y) # uses N-1 denominator
  theoretical_var <- N^2 * (1 - n / N) * S2 / n

  set.seed(42)
  ht_estimates <- replicate(nrep, {
    s <- equal_prob_wor(N, n, method = "srs")
    sum(y[s$sample]) * N / n
  })

  expect_equal(
    var(ht_estimates),
    theoretical_var,
    tolerance = 0.1 * theoretical_var
  )
})


# Property 12: WR sampling covariance has non-negative diagonal (variance)
test_that("WR sampling covariance has non-negative diagonal (variance)", {
  # Multinomial with E(n_i) > 1
  prob <- c(0.05, 0.15, 0.30, 0.50)
  hits <- 10 * prob
  s <- unequal_prob_wr(hits, method = "multinomial")
  cov_mat <- sampling_cov(s)
  expect_true(all(diag(cov_mat) >= 0))

  # SRS WR with E(n_i) > 1
  s2 <- equal_prob_wr(5, 10)
  cov_mat2 <- sampling_cov(s2)
  expect_true(all(diag(cov_mat2) >= 0))
})


# Property 13: print works for non-integer expected n
test_that("print works for Poisson/Bernoulli with non-integer n", {
  s1 <- unequal_prob_wor(c(0.2, 0.3, 0.4), method = "poisson")
  expect_output(print(s1), "n=0.9")

  s2 <- equal_prob_wor(10, 3.5, method = "bernoulli")
  expect_output(print(s2), "n=3.5")

  # Integer n still prints cleanly
  s3 <- equal_prob_wor(10, 3, method = "srs")
  expect_output(print(s3), "n=3,")
})


# Property 14: sampling_cov.wr weighted handles zero-prob units
test_that("sampling_cov.wr weighted returns NA (not NaN) for zero-prob units", {
  s <- unequal_prob_wr(c(0, 1, 1), method = "multinomial")
  m <- sampling_cov(s, weighted = TRUE)
  # NA for zero-prob unit, not NaN
  expect_true(all(is.na(m[1, ])))
  expect_true(all(is.na(m[, 1])))
  expect_false(any(is.nan(m)))
  # Non-zero prob entries should be finite
  expect_true(all(is.finite(m[2:3, 2:3])))
})


# Property 15: sampling_cov.wor weighted handles zero-prob units
test_that("sampling_cov.wor weighted returns NA (not NaN) for zero-prob units", {
  s <- unequal_prob_wor(c(0, 0.5, 0.5, 0.5, 0.5), method = "cps")
  m <- sampling_cov(s, weighted = TRUE)
  # No NaN anywhere
  expect_false(any(is.nan(m)))
  # Off-diagonals involving zero-prob unit should be NA
  expect_true(all(is.na(m[1, -1])))
  expect_true(all(is.na(m[-1, 1])))
  # Diagonal for zero-prob unit: 1 - 0 = 1
  expect_equal(m[1, 1], 1)
  # Non-zero-prob entries should be finite
  expect_true(all(is.finite(m[2:5, 2:5])))
})


# Property 16: Joint probability functions reject large N
test_that("joint probability functions reject large N", {
  pik <- rep(0.5, 10001)
  expect_error(
    joint_inclusion_prob(unequal_prob_wor(pik, method = "poisson")),
    "sampled_only = TRUE"
  )
})

# Property 17: sampled_only batch error
test_that("sampled_only errors for batch designs with subsetting instructions", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "cps", nrep = 3)
  expect_error(
    joint_inclusion_prob(s, sampled_only = TRUE),
    "pikl\\[s\\$sample\\[, i\\]"
  )
})
