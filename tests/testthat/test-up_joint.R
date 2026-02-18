test_that("joint_inclusion_prob.poisson returns correct structure", {
  pik <- c(0.2, 0.5, 0.8)
  s <- unequal_prob_wor(pik, method = "poisson")
  pikl <- joint_inclusion_prob(s)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl))
  expect_equal(diag(pikl), pik)
})

test_that("joint_inclusion_prob.poisson computes pi_i * pi_j off-diagonal", {
  pik <- c(0.2, 0.5, 0.8)
  s <- unequal_prob_wor(pik, method = "poisson")
  pikl <- joint_inclusion_prob(s)
  expect_equal(pikl[1, 2], 0.2 * 0.5)
  expect_equal(pikl[1, 3], 0.2 * 0.8)
  expect_equal(pikl[2, 3], 0.5 * 0.8)
})

test_that("joint_inclusion_prob.poisson does not require integer sum(pik)", {
  pik <- c(0.2, 0.3, 0.4)
  s <- unequal_prob_wor(pik, method = "poisson")
  expect_no_warning(joint_inclusion_prob(s))
})

test_that("joint_inclusion_prob.cps returns symmetric matrix with correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps")
  pikl <- joint_inclusion_prob(s)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})

test_that("joint_inclusion_prob.brewer returns symmetric matrix with correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "brewer")
  pikl <- joint_inclusion_prob(s)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})

# ---- joint_inclusion_prob for systematic ----

test_that("joint_inclusion_prob.systematic returns symmetric matrix", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "systematic")
  pikl <- joint_inclusion_prob(s)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})

# ---- joint_inclusion_prob for equal probability designs ----

test_that("joint_inclusion_prob.srs returns known formula", {
  s <- equal_prob_wor(10, 3)
  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(10, 10))
  expect_equal(diag(pikl), rep(3 / 10, 10))
  expect_equal(pikl[1, 2], 3 * 2 / (10 * 9))
})

test_that("joint_inclusion_prob.bernoulli returns p^2 off-diagonal", {
  s <- equal_prob_wor(5, 1.5, method = "bernoulli")
  pikl <- joint_inclusion_prob(s)
  p <- 1.5 / 5
  expect_equal(dim(pikl), c(5, 5))
  expect_equal(diag(pikl), rep(p, 5))
  expect_equal(pikl[1, 2], p * p)
})

# ---- joint_expected_hits for chromy ----

test_that("joint_expected_hits.chromy returns matrix of correct size", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 2)
  s <- unequal_prob_wr(hits, method = "chromy")
  pikl <- joint_expected_hits(s, nsim = 1000)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 0.05)
})

test_that("joint_expected_hits.chromy rejects non-integer nsim", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 1)
  s <- unequal_prob_wr(hits, method = "chromy")
  expect_error(joint_expected_hits(s, nsim = 100.7), "not close to an integer")
})

test_that("joint_expected_hits.chromy rejects Inf nsim", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 1)
  s <- unequal_prob_wr(hits, method = "chromy")
  expect_error(joint_expected_hits(s, nsim = Inf), "finite")
})

test_that("joint_expected_hits.chromy accepts integer-like nsim", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 2)
  s <- unequal_prob_wr(hits, method = "chromy")
  expect_no_error(joint_expected_hits(s, nsim = 100.0))
})

# ---- joint_expected_hits for multinomial ----

test_that("joint_expected_hits.multinomial returns correct values", {
  x <- c(10, 20, 30, 40)
  hits <- expected_hits(x, n = 5)
  s <- unequal_prob_wr(hits, method = "multinomial")
  pikl <- joint_expected_hits(s)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(4, 4))

  p <- x / sum(x)
  n <- 5
  # Off-diagonal: n*(n-1)*p_i*p_j
  expect_equal(pikl[1, 2], n * (n - 1) * p[1] * p[2])
  # Diagonal: n*p_k
  expect_equal(diag(pikl), n * p)
})

# ---- joint_expected_hits for srs_wr ----

test_that("joint_expected_hits.srs_wr returns correct values", {
  s <- equal_prob_wr(10, 5)
  pikl <- joint_expected_hits(s)
  expect_equal(dim(pikl), c(10, 10))
  p <- 1 / 10
  expect_equal(diag(pikl), rep(5 * p, 10))
  expect_equal(pikl[1, 2], 5 * 4 * p * p)
})

# ---- sampling_cov ----

test_that("sampling_cov.wor returns correct matrix", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps")
  delta <- sampling_cov(s)
  pikl <- joint_inclusion_prob(s)

  expect_equal(delta, pikl - outer(pik, pik), tolerance = 1e-10)
  # Diagonal should be pik * (1 - pik)
  expect_equal(diag(delta), pik * (1 - pik), tolerance = 1e-10)
})

test_that("sampling_cov works for all WOR design types", {
  expect_no_error(sampling_cov(equal_prob_wor(10, 3)))
  expect_no_error(sampling_cov(unequal_prob_wor(
    c(0.2, 0.3, 0.5),
    method = "brewer"
  )))
  expect_no_error(sampling_cov(unequal_prob_wor(
    c(0.2, 0.3, 0.5),
    method = "cps"
  )))
  expect_no_error(sampling_cov(unequal_prob_wor(
    c(0.2, 0.3, 0.5),
    method = "systematic"
  )))
  expect_no_error(sampling_cov(unequal_prob_wor(
    c(0.2, 0.3, 0.5),
    method = "poisson"
  )))
})

test_that("sampling_cov works for WR designs", {
  x <- c(10, 20, 30)
  hits <- expected_hits(x, n = 2)
  expect_no_error(sampling_cov(
    unequal_prob_wr(hits, method = "chromy"),
    nsim = 500
  ))
  expect_no_error(sampling_cov(unequal_prob_wr(hits, method = "multinomial")))
  expect_no_error(sampling_cov(equal_prob_wr(10, 3)))
})

# ---- sampling_cov weighted (SYG check) ----

test_that("sampling_cov weighted returns correct matrix for WOR", {
  pik <- c(0.2, 0.4, 0.6, 0.8) # n = 2, all pi_ij > 0
  s <- unequal_prob_wor(pik, method = "cps")
  chk <- sampling_cov(s, weighted = TRUE)
  pikl <- joint_inclusion_prob(s)

  expected <- 1 - outer(pik, pik) / pikl
  diag(expected) <- 1 - pik
  expect_equal(chk, expected, tolerance = 1e-10)
})

test_that("sampling_cov weighted warns and returns NA when pi_ij = 0", {
  pik <- c(0.2, 0.3, 0.5) # n = 1, all off-diagonal pi_ij = 0
  s <- unequal_prob_wor(pik, method = "cps")

  expect_warning(
    chk <- sampling_cov(s, weighted = TRUE),
    "Sen-Yates-Grundy"
  )

  # diagonal is 1 - pik
  expect_equal(diag(chk), 1 - pik)
  # off-diagonal entries are NA
  off_diag <- chk[row(chk) != col(chk)]
  expect_true(all(is.na(off_diag)))
})

test_that("sampling_cov weighted off-diagonal is non-positive for high-entropy designs", {
  pik <- c(0.2, 0.4, 0.6, 0.8) # n = 2

  chk_cps <- sampling_cov(unequal_prob_wor(pik, method = "cps"), weighted = TRUE)
  off_diag <- chk_cps[row(chk_cps) != col(chk_cps)]
  expect_true(all(off_diag <= 1e-10))
  expect_equal(diag(chk_cps), 1 - pik, tolerance = 1e-10)

  chk_brewer <- sampling_cov(
    unequal_prob_wor(pik, method = "brewer"),
    weighted = TRUE
  )
  off_diag_b <- chk_brewer[row(chk_brewer) != col(chk_brewer)]
  expect_true(all(off_diag_b <= 1e-10))
  expect_equal(diag(chk_brewer), 1 - pik, tolerance = 1e-10)
})
