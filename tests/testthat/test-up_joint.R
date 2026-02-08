test_that("up_poisson_jip returns correct structure", {
  pik <- c(0.2, 0.5, 0.8)
  pikl <- up_poisson_jip(pik)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl))  # symmetric
  expect_equal(diag(pikl), pik)
})

test_that("up_poisson_jip computes pi_i * pi_j off-diagonal", {
  pik <- c(0.2, 0.5, 0.8)
  pikl <- up_poisson_jip(pik)
  expect_equal(pikl[1, 2], 0.2 * 0.5)
  expect_equal(pikl[1, 3], 0.2 * 0.8)
  expect_equal(pikl[2, 3], 0.5 * 0.8)
})

test_that("up_poisson_jip does not require integer sum(pik)", {
  pik <- c(0.2, 0.3, 0.4)  # sum = 0.9
  expect_no_warning(up_poisson_jip(pik))
})

# Fixed-size joint wrappers should validate sum(pik)

test_that("up_maxent_jip rejects non-integer sum(pik)", {
  expect_error(up_maxent_jip(c(0.49, 0.49, 0.49)), "not close to an integer")
})

test_that("up_maxent_jip silently accepts integer sum(pik)", {
  expect_no_error(up_maxent_jip(c(0.2, 0.3, 0.5)))
})

test_that("up_systematic_jip rejects non-integer sum(pik)", {
  expect_error(up_systematic_jip(c(0.49, 0.49, 0.49)), "not close to an integer")
})

test_that("up_systematic_jip silently accepts integer sum(pik)", {
  expect_no_error(up_systematic_jip(c(0.2, 0.3, 0.5)))
})

test_that("up_brewer_jip rejects non-integer sum(pik)", {
  expect_error(up_brewer_jip(c(0.49, 0.49, 0.49)), "not close to an integer")
})

test_that("up_brewer_jip silently accepts integer sum(pik)", {
  expect_no_error(up_brewer_jip(c(0.2, 0.3, 0.5)))
})

# Basic validation shared across all joint wrappers

test_that("joint wrappers reject invalid pik", {
  expect_error(up_poisson_jip(c(0.5, NA)), "missing values")
  expect_error(up_maxent_jip(c(0.5, NA)), "missing values")
  expect_error(up_systematic_jip(c(0.5, NA)), "missing values")
  expect_error(up_brewer_jip(c(0.5, NA)), "missing values")

  expect_error(up_poisson_jip(c("a", "b")), "numeric")
  expect_error(up_maxent_jip(c("a", "b")), "numeric")
  expect_error(up_systematic_jip(c("a", "b")), "numeric")
  expect_error(up_brewer_jip(c("a", "b")), "numeric")
})

# Structural checks for fixed-size joint wrappers

test_that("up_maxent_jip returns symmetric matrix with correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_maxent_jip(pik)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})

test_that("up_systematic_jip returns symmetric matrix with correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_systematic_jip(pik)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})

# up_chromy_pairexp integer checks

test_that("up_chromy_pairexp rejects non-integer n", {
  expect_error(up_chromy_pairexp(c(10, 20, 30), n = 1.7), "not close to an integer")
})

test_that("up_chromy_pairexp rejects Inf n", {
  expect_error(up_chromy_pairexp(c(10, 20, 30), n = Inf), "finite")
})

test_that("up_chromy_pairexp rejects non-integer nsim", {
  expect_error(up_chromy_pairexp(c(10, 20, 30), n = 1, nsim = 100.7), "not close to an integer")
})

test_that("up_chromy_pairexp rejects Inf nsim", {
  expect_error(up_chromy_pairexp(c(10, 20, 30), n = 1, nsim = Inf), "finite")
})

test_that("up_chromy_pairexp accepts integer-like doubles", {
  expect_no_error(up_chromy_pairexp(c(10, 20, 30), n = 2.0, nsim = 100.0))
})

test_that("up_brewer_jip returns symmetric matrix with correct diagonal", {
  pik <- c(0.2, 0.3, 0.5)
  pikl <- up_brewer_jip(pik)
  expect_true(is.matrix(pikl))
  expect_equal(dim(pikl), c(3, 3))
  expect_equal(pikl, t(pikl), tolerance = 1e-10)
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
})
