test_that("inclusion_prob returns correct length", {
  pik <- inclusion_prob(c(10, 20, 30, 40), n = 2)
  expect_length(pik, 4)
})

test_that("inclusion_prob sums to n", {
  pik <- inclusion_prob(c(10, 20, 30, 40), n = 2)
  expect_equal(sum(pik), 2)

  pik <- inclusion_prob(1:100, n = 25)
  expect_equal(sum(pik), 25)
})

test_that("inclusion_prob values are proportional to size", {
  size <- c(10, 20, 30, 40)
  pik <- inclusion_prob(size, n = 2)

  expect_equal(pik[2] / pik[1], 2, tolerance = 1e-6)
  expect_equal(pik[3] / pik[1], 3, tolerance = 1e-6)
  expect_equal(pik[4] / pik[1], 4, tolerance = 1e-6)
})

test_that("inclusion_prob handles certainty selections", {
  # One very large unit
  size <- c(1, 1, 1, 100)
  pik <- inclusion_prob(size, n = 2)

  expect_equal(pik[4], 1) # Large unit gets probability 1
  expect_equal(sum(pik[1:3]), 1) # Remaining units share 1
})

test_that("inclusion_prob handles multiple certainty selections", {
  size <- c(1, 100, 200, 300) # Three large units
  pik <- inclusion_prob(size, n = 3)

  expect_true(pik[2] == 1 || pik[3] == 1 || pik[4] == 1)
})

test_that("inclusion_prob returns values in [0, 1]", {
  set.seed(42)
  size <- runif(100, 1, 1000)
  pik <- inclusion_prob(size, n = 30)

  expect_true(all(pik >= 0))
  expect_true(all(pik <= 1))
})

test_that("inclusion_prob handles equal sizes", {
  size <- rep(1, 10)
  pik <- inclusion_prob(size, n = 5)

  expect_equal(pik, rep(0.5, 10))
})

test_that("inclusion_prob warns about negative values", {
  expect_warning(
    pik <- inclusion_prob(c(-1, 2, 3, 4), n = 2),
    "negative"
  )
})

test_that("inclusion_prob rejects invalid input", {
  expect_error(inclusion_prob(1:10, n = 15), "cannot exceed")
  expect_error(inclusion_prob(1:10, n = -1), "non-negative")
  expect_error(inclusion_prob(1:10, n = NA), "single numeric")
  expect_error(inclusion_prob(1:10, n = NA_real_), "not NA")
})

test_that("inclusion_prob rejects non-numeric n", {
  expect_error(inclusion_prob(1:10, n = "5"), "single numeric")
  expect_error(inclusion_prob(1:10, n = TRUE), "single numeric")
})

test_that("inclusion_prob rejects vector n", {
  expect_error(inclusion_prob(1:10, n = c(1, 2)), "single numeric")
})

test_that("inclusion_prob rejects non-numeric x", {
  expect_error(inclusion_prob(c("a", "b", "c"), n = 1), "numeric vector")
})

test_that("inclusion_prob rejects non-integer n", {
  expect_error(inclusion_prob(1:10, n = 2.9), "not close to an integer")
})

test_that("inclusion_prob rejects Inf n", {
  expect_error(inclusion_prob(1:10, n = Inf), "finite")
})

test_that("inclusion_prob silently accepts integer-like n", {
  expect_no_error(inclusion_prob(1:10, n = 3.0))
})

test_that("inclusion_prob works with unequal_prob_wor", {
  size <- c(500, 1200, 800, 3000, 600)
  pik <- inclusion_prob(size, n = 3)

  set.seed(42)
  s <- unequal_prob_wor(pik, method = "cps")
  expect_length(s$sample, 3)
})

test_that("inclusion_prob.wor extracts pik from design", {
  pik <- c(0.2, 0.3, 0.5)
  s <- unequal_prob_wor(pik, method = "cps")
  expect_equal(inclusion_prob(s), pik)
})

test_that("expected_hits.default computes n * x / sum(x)", {
  x <- c(10, 20, 30, 40)
  expect_equal(expected_hits(x, n = 3), 3 * x / sum(x))
})

test_that("expected_hits.default allows n > length(x) (WR context)", {
  x <- c(10, 20, 30, 40)
  hits <- expected_hits(x, n = 10)
  expect_equal(sum(hits), 10)
  expect_equal(hits, 10 * x / sum(x))
})

test_that("expected_hits.wr extracts from design", {
  x <- c(10, 20, 30, 40)
  hits <- expected_hits(x, n = 3)
  s <- unequal_prob_wr(hits, method = "chromy")
  expect_equal(expected_hits(s), hits)
})

# expected_hits input validation

test_that("expected_hits errors when n is missing", {
  expect_error(expected_hits(c(10, 20, 30)), "required")
})

test_that("expected_hits errors when n is not a single numeric", {
  expect_error(expected_hits(c(10, 20), n = "a"), "single numeric")
  expect_error(expected_hits(c(10, 20), n = c(1, 2)), "single numeric")
})

test_that("expected_hits errors when n is NA or negative", {
  expect_error(expected_hits(c(10, 20), n = NA_real_), "non-negative")
  expect_error(expected_hits(c(10, 20), n = -1), "non-negative")
})

test_that("expected_hits errors when x is not numeric", {
  expect_error(expected_hits("a", n = 2), "numeric vector")
})

test_that("expected_hits rejects NA in x", {
  expect_error(expected_hits(c(10, NA, 30), n = 2), "missing values")
})

test_that("expected_hits rejects Inf in x", {
  expect_error(expected_hits(c(10, Inf, 30), n = 2), "finite")
})

test_that("expected_hits rejects negative x", {
  expect_error(expected_hits(c(10, -5, 30), n = 2), "non-negative")
})

test_that("expected_hits errors when sum(x) is zero", {
  expect_error(expected_hits(c(0, 0, 0), n = 2), "sum")
})

# Non-finite and NA input validation

test_that("inclusion_prob rejects Inf in x", {
  expect_error(inclusion_prob(c(1, Inf, 2), 2), "finite")
})

test_that("inclusion_prob rejects -Inf in x", {
  expect_error(inclusion_prob(c(1, -Inf, 2), 2), "finite")
})

test_that("inclusion_prob rejects NaN in x", {
  expect_error(inclusion_prob(c(1, NaN, 2), 2), "missing values")
})

test_that("inclusion_prob rejects NA in x", {
  expect_error(inclusion_prob(c(1, NA, 2), 2), "missing values")
})

# n exceeds achievable positive units

test_that("inclusion_prob errors when n exceeds positive units", {
  expect_error(inclusion_prob(c(0, 0, 0, 1), 3), "exceeds")
})

test_that("inclusion_prob errors when all zeros with n > 0", {
  expect_error(inclusion_prob(c(0, 0, 0), 1), "exceeds")
})

test_that("inclusion_prob handles n = 0 with all-zero x", {
  pik <- inclusion_prob(c(0, 0, 0), n = 0)
  expect_equal(pik, c(0, 0, 0))
})

test_that("inclusion_prob correctly sums to n after certainty capping", {
  # 2 very large units, 2 small, n = 3
  # After capping the 2 large ones, remaining n = 1 is achievable
  size <- c(1, 1, 100, 200)
  pik <- inclusion_prob(size, n = 3)
  expect_equal(sum(pik), 3)
  expect_equal(pik[3], 1)
  expect_equal(pik[4], 1)
})
