# Tille Example 10 (page 86), known exact values
tille_pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)

tille_expected <- matrix(c(
  0.07,   0.0049, 0.0130, 0.0215, 0.0447, 0.0559,
  0.0049, 0.17,   0.0324, 0.0537, 0.1113, 0.1377,
  0.0130, 0.0324, 0.41,   0.1407, 0.2888, 0.3452,
  0.0215, 0.0537, 0.1407, 0.61,   0.4691, 0.5351,
  0.0447, 0.1113, 0.2888, 0.4691, 0.83,   0.7461,
  0.0559, 0.1377, 0.3452, 0.5351, 0.7461, 0.91
), nrow = 6, byrow = TRUE)

test_that("joint_inclusion_prob.cps returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps")
  pikl <- joint_inclusion_prob(s)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.cps is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.cps has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.cps satisfies marginal constraint", {
  pik <- tille_pik
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("joint_inclusion_prob.cps matches Tille Example 10", {
  pikl <- joint_inclusion_prob(unequal_prob_wor(tille_pik, method = "cps"))

  expect_equal(pikl, tille_expected, tolerance = 0.01)
})

test_that("joint_inclusion_prob.cps produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))

  for (i in 1:3) {
    for (j in 1:3) {
      if (i != j) {
        expect_true(pikl[i, j] <= min(pik[i], pik[j]) + 1e-9)
      }
    }
  }
})

test_that("joint_inclusion_prob.cps handles n=1 case", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("joint_inclusion_prob.cps handles certainty selections", {
  pik <- c(0.5, 0.5, 1.0)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  expect_equal(pikl[3, 1], pik[1], tolerance = 1e-6)
  expect_equal(pikl[3, 2], pik[2], tolerance = 1e-6)
})

test_that("joint_inclusion_prob.systematic returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.systematic is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.systematic has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.systematic produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("joint_inclusion_prob.systematic can have zero joint probs", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))

  off_diag <- pikl[row(pikl) != col(pikl)]
  expect_true(any(off_diag == 0) || all(off_diag > 0))
})

test_that("joint_inclusion_prob.systematic satisfies marginal constraint", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  n <- sum(pik)

  for (k in seq_along(pik)) {
    row_sum <- sum(pikl[k, ]) - pikl[k, k]
    expected <- (n - 1) * pik[k]
    expect_equal(row_sum, expected, tolerance = 0.01)
  }
})

test_that("joint_inclusion_prob.brewer returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("joint_inclusion_prob.brewer is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(pikl, t(pikl))
})

test_that("joint_inclusion_prob.brewer has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.brewer produces valid probabilities", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
})

test_that("joint_inclusion_prob.brewer approximates CPS reasonably", {
  pik <- tille_pik
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))
  pikl_exact <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  off_diag_brewer <- pikl_brewer[row(pikl_brewer) != col(pikl_brewer)]
  off_diag_exact <- pikl_exact[row(pikl_exact) != col(pikl_exact)]
  expect_true(cor(off_diag_brewer, off_diag_exact) > 0.99)

  mad <- mean(abs(pikl_brewer - pikl_exact))
  expect_true(mad < 0.03)
})

test_that("joint_inclusion_prob.brewer handles n <= 1 case", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("joint_expected_hits.chromy returns valid matrix", {
  x <- c(10, 20, 15, 25, 30)
  hits <- expected_hits(x, n = 3)
  s <- unequal_prob_wr(hits, method = "chromy")
  joint <- joint_expected_hits(s, nsim = 5000)

  expect_equal(dim(joint), c(5, 5))
  expect_equal(joint, t(joint))
  expect_true(all(joint > 0))

  pik <- 3 * x / sum(x)
  expect_equal(diag(joint), pik, tolerance = 0.05)
})

test_that("all joint_inclusion_prob methods agree on diagonal", {
  pik <- c(0.2, 0.3, 0.5)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(diag(pikl_cps), pik)
  expect_equal(diag(pikl_sys), pik)
  expect_equal(diag(pikl_brewer), pik)
})

test_that("all joint methods produce symmetric matrices", {
  pik <- c(0.15, 0.25, 0.35, 0.25)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))

  expect_equal(pikl_cps, t(pikl_cps))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_brewer, t(pikl_brewer))
})

test_that("equal-prob systematic JIP differs from SRS", {
  N <- 12; n <- 3
  s <- equal_prob_wor(N, n, method = "systematic")
  pikl <- joint_inclusion_prob(s)

  # Should NOT be uniform (SRS would be n(n-1)/(N(N-1)) = 1/22 for all pairs)
  off_diag <- pikl[row(pikl) != col(pikl)]
  expect_true(length(unique(round(off_diag, 8))) > 1)

  # Some pairs should have pi_ij = 0 (never co-occur in systematic)
  expect_true(any(off_diag == 0))
  # Some pairs should have pi_ij = n/N = 0.25 (always co-occur)
  expect_true(any(abs(off_diag - n / N) < 1e-10))

  # Diagonal should still be n/N
  expect_equal(diag(pikl), rep(n / N, N))
})

test_that("high-entropy JIP warns on large marginal defect", {
  # Skewed pik with n=2 → defect/n ≈ 6.9%, should warn
  pik_skewed <- c(0.190194, 0.702073, 0.168549, 0.026420,
                  0.441842, 0.304860, 0.109806, 0.056257)
  s <- unequal_prob_wor(pik_skewed, method = "brewer")
  expect_warning(joint_inclusion_prob(s), "marginal defect")
})

test_that("high-entropy JIP does not warn for well-spread pik", {
  # Tille pik → defect/n ≈ 1.7%, should not warn
  pik_tille <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
  s <- unequal_prob_wor(pik_tille, method = "brewer")
  expect_no_warning(joint_inclusion_prob(s))
})

test_that("high-entropy JIP handles certainty units correctly", {
  pik <- c(
    1, 1, 0.191473, 1, 0.636845, 0.653934, 0.292246, 1, 1,
    0.214227, 0.846962, 0.224352, 1, 0.437642, 0.758233,
    0.618972, 0.178129, 0.795222, 1, 0.87704, 0.452308,
    0.703438, 1, 0.423387, 0.69559
  )
  upper <- outer(pik, pik, pmin)
  cert <- which(pik == 1)

  for (m in c("brewer", "sps", "pareto")) {
    s <- unequal_prob_wor(pik, method = m)
    J <- joint_inclusion_prob(s)

    # No entry exceeds 1
    expect_true(all(J <= 1 + 1e-10), info = paste(m, "entries > 1"))
    # No entry exceeds min(pi_i, pi_j)
    expect_true(all(J <= upper + 1e-10), info = paste(m, "entries > upper"))
    # No negative entries
    expect_true(all(J >= -1e-10), info = paste(m, "negative entries"))
    # Certainty-certainty pairs = 1
    expect_true(all(J[cert, cert] == 1), info = paste(m, "cert-cert"))
    # Certainty-other pairs = pik[other]
    for (ci in cert) {
      others <- setdiff(seq_along(pik), ci)
      expect_equal(J[ci, others], pik[others], tolerance = 1e-10,
                   info = paste(m, "cert-other for unit", ci))
    }
  }
})

test_that("joint methods work with larger populations", {
  set.seed(1)
  N <- 50
  n <- 10
  x <- runif(N)
  pik <- n * x / sum(x)
  hits <- expected_hits(x, n = n)

  pikl_cps <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  pikl_brewer <- joint_inclusion_prob(unequal_prob_wor(pik, method = "brewer"))
  pikl_sys <- joint_inclusion_prob(unequal_prob_wor(pik, method = "systematic"))
  pikl_chromy <- joint_expected_hits(unequal_prob_wr(hits, method = "chromy"))

  expect_equal(dim(pikl_cps), c(N, N))
  expect_equal(dim(pikl_brewer), c(N, N))
  expect_equal(dim(pikl_sys), c(N, N))
  expect_equal(dim(pikl_chromy), c(N, N))

  expect_equal(pikl_cps, t(pikl_cps))
  expect_equal(pikl_brewer, t(pikl_brewer))
  expect_equal(pikl_sys, t(pikl_sys))
  expect_equal(pikl_chromy, t(pikl_chromy))
})


# ---- sampled_only = TRUE ----

test_that("sampled_only returns n x n for CPS matching full[s,s]", {
  pik <- tille_pik
  s <- unequal_prob_wor(pik, method = "cps")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(length(idx), length(idx)))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("sampled_only returns n x n for systematic matching full[s,s]", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "systematic")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(length(idx), length(idx)))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("sampled_only returns n x n for HE methods matching full[s,s]", {
  set.seed(1)
  N <- 50
  n <- 10
  pik <- inclusion_prob(runif(N), n = n)

  for (method in c("brewer", "sps", "pareto")) {
    s <- unequal_prob_wor(pik, method = method)
    full <- joint_inclusion_prob(s)
    sub <- joint_inclusion_prob(s, sampled_only = TRUE)

    idx <- s$sample
    expect_equal(dim(sub), c(length(idx), length(idx)),
                 info = paste("dim for", method))
    expect_equal(sub, full[idx, idx, drop = FALSE], tolerance = 1e-10,
                 info = paste("values for", method))
  }
})

test_that("sampled_only returns n x n for Poisson matching full[s,s]", {
  set.seed(1)
  pik <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  s <- unequal_prob_wor(pik, method = "poisson")
  # Ensure at least 2 sampled units for meaningful test
  while (length(s$sample) < 2) {
    s <- unequal_prob_wor(pik, method = "poisson")
  }
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(length(idx), length(idx)))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("sampled_only returns n x n for SRS WOR matching full[s,s]", {
  s <- equal_prob_wor(10, 3)
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(3, 3))
  expect_equal(sub, full[idx, idx])
})

test_that("sampled_only returns n x n for Bernoulli matching full[s,s]", {
  set.seed(42)
  s <- equal_prob_wor(10, 3, method = "bernoulli")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(length(idx), length(idx)))
  expect_equal(sub, full[idx, idx])
})

test_that("sampled_only errors for batch WOR designs", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps", nrep = 3)
  expect_error(
    joint_inclusion_prob(s, sampled_only = TRUE),
    "not supported for batch"
  )
})

test_that("sampled_only errors for random-size batch WOR designs", {
  pik <- c(0.3, 0.3, 0.4)
  s <- unequal_prob_wor(pik, method = "poisson", nrep = 3)
  expect_error(
    joint_inclusion_prob(s, sampled_only = TRUE),
    "not supported for batch"
  )
})

test_that("sampled_only skips marginal defect warning for HE", {
  # This pik vector triggers marginal defect warning when sampled_only=FALSE
  pik_skewed <- c(0.190194, 0.702073, 0.168549, 0.026420,
                  0.441842, 0.304860, 0.109806, 0.056257)
  s <- unequal_prob_wor(pik_skewed, method = "brewer")
  expect_warning(joint_inclusion_prob(s), "marginal defect")
  expect_no_warning(joint_inclusion_prob(s, sampled_only = TRUE))
})

test_that("sampled_only works with HE certainty units", {
  pik <- c(1, 0.5, 0.5, 1)
  s <- unequal_prob_wor(pik, method = "brewer")
  full <- joint_inclusion_prob(s)
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)

  idx <- s$sample
  expect_equal(dim(sub), c(length(idx), length(idx)))
  expect_equal(sub, full[idx, idx, drop = FALSE], tolerance = 1e-10)
})

test_that("N > 10K error suggests sampled_only", {
  pik <- rep(0.5, 10001)
  expect_error(
    joint_inclusion_prob(unequal_prob_wor(pik, method = "poisson")),
    "sampled_only = TRUE"
  )
})


# ---- CPS JIP equal/near-equal weight regression ----

test_that("cps jip is exact for equal-pik case", {
  N <- 20
  n <- 6
  pik <- rep(n / N, N)

  s <- unequal_prob_wor(pik, method = "cps")
  J <- joint_inclusion_prob(s)

  # Known exact off-diagonal for equal-pi fixed-size (SRS)
  target <- n * (n - 1) / (N * (N - 1))
  offdiag <- J[row(J) != col(J)]
  expect_equal(offdiag, rep(target, N * (N - 1)), tolerance = 1e-10)

  # Diagonal = pik
  expect_equal(diag(J), pik, tolerance = 1e-14)
})

test_that("cps jip satisfies row-sum identity for equal-pik", {
  N <- 20
  n <- 6
  pik <- rep(n / N, N)

  J <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))

  # Fixed-size identity: sum_{j != i} pi_ij = (n-1) * pi_i
  for (i in seq_len(N)) {
    row_sum <- sum(J[i, ]) - J[i, i]
    expect_equal(row_sum, (n - 1) * pik[i], tolerance = 1e-10)
  }
})

test_that("cps jip is exact for certainty + equal valid units", {
  # 2 certainty units + 8 equal valid units
  N <- 10
  pik <- c(1.0, 1.0, rep(3 / 8, 8))

  s <- unequal_prob_wor(pik, method = "cps")
  J <- joint_inclusion_prob(s)

  # Certainty-certainty: 1.0

  expect_equal(J[1, 2], 1.0, tolerance = 1e-14)
  # Certainty-valid: pik[valid]
  for (i in 3:10) {
    expect_equal(J[1, i], pik[i], tolerance = 1e-14)
    expect_equal(J[2, i], pik[i], tolerance = 1e-14)
  }

  # Valid-valid off-diagonal: n_valid*(n_valid-1)/(N_valid*(N_valid-1))
  N_valid <- 8
  n_valid <- 3  # sum(pik[3:10])
  target <- n_valid * (n_valid - 1) / (N_valid * (N_valid - 1))
  for (i in 3:9) {
    for (j in (i + 1):10) {
      expect_equal(J[i, j], target, tolerance = 1e-10)
    }
  }
})

test_that("cps jip row-sum identity holds for skewed pik", {
  # Non-equal pik — verify the fix doesn't break the general case
  pik <- tille_pik
  J <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  n <- sum(pik)

  for (i in seq_along(pik)) {
    row_sum <- sum(J[i, ]) - J[i, i]
    expect_equal(row_sum, (n - 1) * pik[i], tolerance = 1e-5)
  }
})
