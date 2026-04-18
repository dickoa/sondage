# Tille Example 10 (page 86), known exact values
tille_pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)

test_that("he_jip returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- he_jip(pik)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("he_jip is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- he_jip(pik)

  expect_equal(pikl, t(pikl))
})

test_that("he_jip has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- he_jip(pik)

  expect_equal(diag(pikl), pik)
})

test_that("he_jip produces valid probabilities", {
  pik <- tille_pik
  pikl <- he_jip(pik)

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
  upper <- outer(pik, pik, pmin)
  expect_true(all(pikl <= upper + 1e-10))
})

test_that("he_jip matches internal dispatch for brewer", {
  pik <- tille_pik
  s <- unequal_prob_wor(pik, method = "brewer")
  from_dispatch <- joint_inclusion_prob(s)
  from_he <- he_jip(pik)

  expect_equal(from_he, from_dispatch)
})

test_that("he_jip approximates CPS reasonably", {
  pik <- tille_pik
  exact <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  approx <- he_jip(pik)

  off_exact <- exact[row(exact) != col(exact)]
  off_approx <- approx[row(approx) != col(approx)]
  expect_true(cor(off_exact, off_approx) > 0.99)
  expect_true(mean(abs(exact - approx)) < 0.03)
})

test_that("he_jip sample_idx returns submatrix matching full[s,s]", {
  pik <- tille_pik
  full <- he_jip(pik)
  idx <- c(1L, 3L, 5L)
  sub <- he_jip(pik, sample_idx = idx)

  expect_equal(dim(sub), c(3, 3))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("he_jip handles certainty units", {
  pik <- c(1, 0.5, 0.5, 1)
  pikl <- he_jip(pik)
  cert <- c(1, 4)

  # Certainty-certainty = 1
  expect_equal(pikl[1, 4], 1)
  # Certainty-other = pik[other]
  expect_equal(pikl[1, 2], pik[2])
  expect_equal(pikl[1, 3], pik[3])
  expect_equal(pikl[4, 2], pik[2])
  expect_equal(pikl[4, 3], pik[3])
})

test_that("he_jip handles n <= 1", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- he_jip(pik)

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("he_jip validates inputs", {
  expect_error(he_jip(c(-0.1, 0.5)), "between 0 and 1")
  expect_error(he_jip(c(0.5, NA)), "missing values")
  expect_error(he_jip("abc"), "numeric")
  expect_error(he_jip(c(0.2, 0.3), eps = -1), "open interval")
})

test_that("he_jip works as joint_fn in register_method", {
  sampford_sample <- function(pik, n = NULL, prn = NULL, ...) {
    N <- length(pik)
    if (n == 1L) {
      return(sample.int(N, 1L, prob = pik))
    }
    q <- pik / (1 - pik)
    repeat {
      first <- sample.int(N, 1L, prob = pik)
      rest <- sample.int(N, n - 1L, replace = TRUE, prob = q)
      s <- c(first, rest)
      if (anyDuplicated(s) == 0L) return(sort(s))
    }
  }

  register_method(
    "test_he_sampford",
    type = "wor",
    sample_fn = sampford_sample,
    joint_fn = he_jip,
    fixed_size = TRUE
  )
  on.exit(unregister_method("test_he_sampford"))

  pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
  s <- unequal_prob_wor(pik, method = "test_he_sampford")

  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(8, 8))
  expect_equal(pikl, t(pikl))
  expect_equal(diag(pikl), pik)

  # sampled_only also works
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)
  idx <- s$sample
  expect_equal(sub, pikl[idx, idx, drop = FALSE])

  # sampling_cov chain works
  cov_mat <- sampling_cov(s)
  expect_equal(dim(cov_mat), c(8, 8))
})

test_that("he_jip with larger population", {
  set.seed(1)
  N <- 100
  n <- 15
  pik <- inclusion_prob(runif(N), n = n)
  pikl <- he_jip(pik)

  expect_equal(dim(pikl), c(N, N))
  expect_equal(pikl, t(pikl))
  expect_equal(diag(pikl), pik)
  expect_true(all(pikl >= -1e-10))
  expect_true(all(pikl <= outer(pik, pik, pmin) + 1e-10))
})

test_that("he_jip with sample_idx and certainty units", {
  pik <- c(1, 0.2, 0.3, 0.5, 1)
  full <- he_jip(pik)
  idx <- c(1L, 2L, 4L, 5L)
  sub <- he_jip(pik, sample_idx = idx)

  expect_equal(dim(sub), c(4, 4))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that(".he_jip_sampled handles sample with no valid units", {
  # all certainty units => valid_s has length 0, hits the early return branch
  pik <- c(1, 1, 1, 0.5, 0.5)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "brewer")
  J <- joint_inclusion_prob(s, sampled_only = TRUE)
  expect_equal(dim(J), c(sum(pik), sum(pik)))
  # certainty rows are filled with pik of sampled units
  expect_true(all(is.finite(J)))
})

test_that("hajek_jip returns correct dimensions", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- hajek_jip(pik)

  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
})

test_that("hajek_jip is symmetric", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- hajek_jip(pik)

  expect_equal(pikl, t(pikl))
})

test_that("hajek_jip has correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- hajek_jip(pik)

  expect_equal(diag(pikl), pik)
})

test_that("hajek_jip produces valid probabilities", {
  pik <- tille_pik
  pikl <- hajek_jip(pik)

  expect_true(all(pikl >= 0))
  expect_true(all(pikl <= 1))
  upper <- outer(pik, pik, pmin)
  expect_true(all(pikl <= upper + 1e-10))
})

test_that("hajek_jip approximates CPS reasonably", {
  pik <- tille_pik
  exact <- joint_inclusion_prob(unequal_prob_wor(pik, method = "cps"))
  approx <- hajek_jip(pik)

  off_exact <- exact[row(exact) != col(exact)]
  off_approx <- approx[row(approx) != col(approx)]
  expect_true(cor(off_exact, off_approx) > 0.99)
  expect_true(mean(abs(exact - approx)) < 0.05)
})

test_that("hajek_jip is close to he_jip", {
  set.seed(1)
  N <- 50
  n <- 10
  pik <- inclusion_prob(runif(N), n = n)
  he <- he_jip(pik)
  hj <- hajek_jip(pik)

  off_he <- he[row(he) != col(he)]
  off_hj <- hj[row(hj) != col(hj)]
  expect_true(cor(off_he, off_hj) > 0.999)
  expect_true(max(abs(he - hj)) < 0.01)
})

test_that("hajek_jip sample_idx returns submatrix matching full[s,s]", {
  pik <- tille_pik
  full <- hajek_jip(pik)
  idx <- c(1L, 3L, 5L)
  sub <- hajek_jip(pik, sample_idx = idx)

  expect_equal(dim(sub), c(3, 3))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("hajek_jip handles certainty units", {
  pik <- c(1, 0.5, 0.5, 1)
  pikl <- hajek_jip(pik)

  # Certainty-certainty = 1
  expect_equal(pikl[1, 4], 1)
  # Certainty-other = pik[other]
  expect_equal(pikl[1, 2], pik[2])
  expect_equal(pikl[1, 3], pik[3])
  expect_equal(pikl[4, 2], pik[2])
  expect_equal(pikl[4, 3], pik[3])
})

test_that("hajek_jip handles n <= 1", {
  pik <- c(0.3, 0.3, 0.4)
  pikl <- hajek_jip(pik)

  expect_equal(diag(pikl), pik)
  off_diag <- pikl
  diag(off_diag) <- 0
  expect_equal(sum(off_diag), 0, tolerance = 1e-9)
})

test_that("hajek_jip validates inputs", {
  expect_error(hajek_jip(c(-0.1, 0.5)), "between 0 and 1")
  expect_error(hajek_jip(c(0.5, NA)), "missing values")
  expect_error(hajek_jip("abc"), "numeric")
  expect_error(hajek_jip(c(0.2, 0.3), eps = -1), "open interval")
})

test_that("hajek_jip works as joint_fn in register_method", {
  simple_pps <- function(pik, n = NULL, prn = NULL, ...) {
    sort(sample.int(length(pik), n, prob = pik))
  }

  register_method(
    "test_hajek_pps",
    type = "wor",
    sample_fn = simple_pps,
    joint_fn = hajek_jip,
    fixed_size = TRUE
  )
  on.exit(unregister_method("test_hajek_pps"))

  pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
  s <- unequal_prob_wor(pik, method = "test_hajek_pps")

  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(8, 8))
  expect_equal(pikl, t(pikl))
  expect_equal(diag(pikl), pik)

  # sampled_only also works
  sub <- joint_inclusion_prob(s, sampled_only = TRUE)
  idx <- s$sample
  expect_equal(sub, pikl[idx, idx, drop = FALSE])

  # Full chain: sampling_cov
  cov_mat <- sampling_cov(s)
  expect_equal(dim(cov_mat), c(8, 8))
})

test_that("hajek_jip with larger population", {
  set.seed(1)
  N <- 100
  n <- 15
  pik <- inclusion_prob(runif(N), n = n)
  pikl <- hajek_jip(pik)

  expect_equal(dim(pikl), c(N, N))
  expect_equal(pikl, t(pikl))
  expect_equal(diag(pikl), pik)
  expect_true(all(pikl >= -1e-10))
  expect_true(all(pikl <= outer(pik, pik, pmin) + 1e-10))
})

test_that("hajek_jip with sample_idx and certainty units", {
  pik <- c(1, 0.2, 0.3, 0.5, 1)
  full <- hajek_jip(pik)
  idx <- c(1L, 2L, 4L, 5L)
  sub <- hajek_jip(pik, sample_idx = idx)

  expect_equal(dim(sub), c(4, 4))
  expect_equal(sub, full[idx, idx, drop = FALSE])
})

test_that("hajek_jip formula is correct analytically", {
  # For equal pi, hajek should give pi_ij = pi^2 * [1 - (1-pi)^2 / (N*pi*(1-pi))]
  N <- 10
  n <- 4
  pik <- rep(n / N, N)
  pikl <- hajek_jip(pik)

  pi <- n / N
  D <- N * pi * (1 - pi)
  expected_ij <- pi^2 * (1 - (1 - pi)^2 / D)

  off_diag <- pikl[1, 2]
  expect_equal(off_diag, expected_ij, tolerance = 1e-12)
})


# ---- both: eps passthrough via ... ----

test_that("eps passes through registered joint_fn dispatch", {
  simple_pps <- function(pik, n = NULL, prn = NULL, ...) {
    sort(sample.int(length(pik), n, prob = pik))
  }

  register_method(
    "test_eps_pass",
    type = "wor",
    sample_fn = simple_pps,
    joint_fn = he_jip,
    fixed_size = TRUE
  )
  on.exit(unregister_method("test_eps_pass"))

  pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
  s <- unequal_prob_wor(pik, method = "test_eps_pass")

  # Default eps = 1e-6
  pikl_default <- joint_inclusion_prob(s)
  # Custom eps passed through ...
  pikl_custom <- joint_inclusion_prob(s, eps = 1e-4)

  # Both should work and be valid matrices
  expect_equal(dim(pikl_default), c(8, 8))
  expect_equal(dim(pikl_custom), c(8, 8))
})
