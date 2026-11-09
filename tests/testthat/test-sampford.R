sampford_brute_jip <- function(pik) {
  N <- length(pik)
  certain <- which(pik == 1)
  active <- which(pik > 0 & pik < 1)
  n_active <- as.integer(round(sum(pik[active])))
  J <- matrix(0, N, N)
  diag(J) <- pik

  for (i in certain) {
    J[i, ] <- pik
    J[, i] <- pik
  }

  if (n_active > 0L && length(active) >= n_active) {
    samples <- combn(active, n_active)
    odds <- pik / (1 - pik)
    prob <- apply(samples, 2L, function(s) {
      prod(odds[s]) * (n_active - sum(pik[s]))
    })
    prob <- prob / sum(prob)
    for (k in seq_len(ncol(samples))) {
      s <- samples[, k]
      J[s, s] <- J[s, s] + (1 - diag(n_active)) * prob[k]
    }
  }
  diag(J) <- pik
  J
}

test_that("sampford returns a fixed-size WOR design", {
  pik <- c(0.2, 0.7, 0.8, 0.5, 0.4, 0.4)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "sampford")

  expect_s3_class(s, c("unequal_prob", "wor", "sondage_sample"))
  expect_identical(s$method, "sampford")
  expect_true(s$fixed_size)
  expect_equal(s$n, 3L)
  expect_length(s$sample, 3L)
  expect_identical(s$sample, sort(unique(s$sample)))
  expect_equal(inclusion_prob(s), pik)
})

test_that("sampford handles equal probabilities and n = 1", {
  set.seed(2)
  s_equal <- unequal_prob_wor(rep(0.5, 8), method = "sampford")
  expect_length(s_equal$sample, 4L)
  expect_identical(s_equal$sample, sort(unique(s_equal$sample)))

  pik_one <- c(0.1, 0.2, 0.3, 0.4)
  s_one <- unequal_prob_wor(pik_one, method = "sampford")
  expect_length(s_one$sample, 1L)
})

test_that("sampford handles complement, zero, and certainty units", {
  set.seed(3)
  pik <- c(1, 0, 0.6, 0.7, 0.8, 0.9)
  s <- unequal_prob_wor(pik, method = "sampford")

  expect_length(s$sample, 4L)
  expect_true(1L %in% s$sample)
  expect_false(2L %in% s$sample)
  expect_identical(s$sample, sort(unique(s$sample)))
})

test_that("sampford exact joint probabilities match enumeration", {
  cases <- list(
    c(0.2, 0.7, 0.8, 0.5, 0.4, 0.4),
    c(0.6, 0.7, 0.8, 0.9),
    c(1, 0, 0.2, 0.3, 0.5),
    c(rep(0.2, 5), rep(0.8, 5))
  )

  for (pik in cases) {
    s <- unequal_prob_wor(pik, method = "sampford")
    J <- joint_inclusion_prob(s)
    expect_equal(J, sampford_brute_jip(pik), tolerance = 2e-12)
    expect_equal(rowSums(J), sum(pik) * pik, tolerance = 2e-11)
    expect_true(isSymmetric(J))
  }
})

test_that("sampford sampled-only joints equal the full submatrix", {
  pik <- c(1, 0, 0.01, 0.19, 0.3, 0.6, 0.9)
  set.seed(4)
  s <- unequal_prob_wor(pik, method = "sampford")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  expect_equal(
    unname(sub),
    full[s$sample, s$sample, drop = FALSE],
    tolerance = 1e-12
  )
  expect_equal(unname(diag(sub)), pik[s$sample])
})

test_that("sampford batch mode returns valid samples", {
  pik <- c(0.2, 0.7, 0.8, 0.5, 0.4, 0.4)
  set.seed(5)
  s <- unequal_prob_wor(pik, method = "sampford", nrep = 100L)

  expect_equal(dim(s$sample), c(3L, 100L))
  expect_true(all(apply(s$sample, 2L, function(x) {
    identical(x, sort(unique(x)))
  })))
})

test_that("sampford empirical inclusion probabilities equal targets", {
  skip_on_cran()
  pik <- c(0.2, 0.7, 0.8, 0.5, 0.4, 0.4)
  set.seed(6)
  nrep <- 20000L
  s <- unequal_prob_wor(pik, method = "sampford", nrep = nrep)
  empirical <- tabulate(as.vector(s$sample), nbins = length(pik)) / nrep
  se <- sqrt(pik * (1 - pik) / nrep)
  expect_true(all(abs(empirical - pik) < 5 * se))
})

test_that("sampford method metadata is available", {
  spec <- method_spec("sampford")
  expect_identical(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
  expect_identical(spec$variance_family, "pps_brewer")
})

## Joint probabilities at moderate population sizes

# Enumeration is limited to small populations, where the polynomial deletion
# could not lose accuracy. These designs broke it while every value stayed
# inside its valid bounds. The tolerances allow for platforms without an
# extended long double.

test_that("sampford joints satisfy the fixed-size identities at moderate N", {
  designs <- list(c(100, 30), c(100, 40), c(100, 50), c(100, 70), c(200, 100))
  for (design in designs) {
    N <- design[1]
    n <- design[2]
    pik <- inclusion_prob(exp(seq(-2, 0, length.out = N)), n)
    s <- unequal_prob_wor(pik, method = "sampford")
    J <- joint_inclusion_prob(s)
    delta <- sampling_cov(s)

    expect_lt(max(abs(rowSums(J) - n * pik)), 1e-8)
    expect_lt(abs(sum(delta)), 1e-7)
    expect_gt(
      min(eigen(delta, symmetric = TRUE, only.values = TRUE)$values),
      -1e-8
    )
  }
})

test_that("sampford joints match a high-precision reference", {
  # Reference values evaluate the design definition with 500 significant
  # digits, P(s) proportional to prod(q[s]) * sum(1 - pik[s]).
  pik <- inclusion_prob(exp(seq(-2, 0, length.out = 100)), 50)
  J <- joint_inclusion_prob(unequal_prob_wor(pik, method = "sampford"))
  expect_equal(J[88, 84], 0.776503241134, tolerance = 1e-9)
  expect_equal(J[47, 88], 0.366945710949758, tolerance = 1e-9)

  pik <- inclusion_prob(exp(seq(-2, 0, length.out = 100)), 70)
  J <- joint_inclusion_prob(unequal_prob_wor(pik, method = "sampford"))
  expect_equal(J[34, 63], 0.470028854100290, tolerance = 1e-9)
})

test_that("sampford sampled-only joints match the full matrix at moderate N", {
  pik <- inclusion_prob(exp(seq(-2, 0, length.out = 100)), 50)
  set.seed(11)
  s <- unequal_prob_wor(pik, method = "sampford")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  expect_lt(max(abs(unname(sub) - unname(full[s$sample, s$sample]))), 1e-10)
})

test_that("sampford sampled-only joints are exact for a dominant unit", {
  # Forward deletion cannot remove unit 1 here, and building the polynomial to
  # full degree would cost more than computing these pairs directly.
  pik <- c(0.9, rep((10 - 0.9) / 24999, 24999))
  set.seed(12)
  s <- unequal_prob_wor(pik, method = "sampford")
  s$sample <- 1:3
  J <- joint_inclusion_prob(s, sampled_only = TRUE)

  expect_equal(J[1, 2], 0.000324012960518421, tolerance = 1e-9)
  expect_equal(J[2, 3], 1.18094170922609e-07, tolerance = 1e-9)
})
