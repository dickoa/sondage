test_that("inclusion_prob returns correct length", {
  pik <- inclusion_prob(c(10, 20, 30, 40), n = 2)
  expect_length(pik, 4)
})

test_that("inclusion_prob preserves unit names", {
  size <- c(north = 10, south = 20, east = 30, west = 40)
  pik <- inclusion_prob(size, n = 2)

  expect_identical(names(pik), names(size))
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
  expect_error(inclusion_prob(1:10, n = NA), "not be NA")
  expect_error(inclusion_prob(1:10, n = NA_real_), "not be NA")
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

test_that("inclusion_prob accepts a fractional expected size", {
  expect_equal(inclusion_prob(rep(1, 10), n = 2.6), rep(0.26, 10))
  expect_equal(inclusion_prob(c(1, 2, 3), n = 1.2), c(0.2, 0.4, 0.6))
  # Capping redistributes and the sum is still n.
  expect_equal(inclusion_prob(c(100, 1, 1), n = 1.6), c(1, 0.3, 0.3))
  expect_equal(
    inclusion_prob(c(1, 1, 1, 100, 100), n = 2.6),
    c(0.2, 0.2, 0.2, 1, 1)
  )
  # Several rounds of capping with a fractional remainder.
  expect_equal(
    inclusion_prob(c(1, 1, 1, 1, 10, 50, 100), n = 3.5),
    c(rep(0.125, 4), 1, 1, 1)
  )
  expect_equal(inclusion_prob(c(0, 0, 0, 5, 5), n = 1.6), c(0, 0, 0, 0.8, 0.8))
})

test_that("inclusion_prob keeps small and near-integer targets", {
  expect_equal(sum(inclusion_prob(1:5, n = 0.2)), 0.2)
  expect_equal(sum(inclusion_prob(1:5, n = 1e-8)), 1e-8)
  expect_equal(sum(inclusion_prob(1:5, n = 1e-12)), 1e-12)
  expect_equal(sum(inclusion_prob(1:10, n = 2.5)), 2.5)
  # No longer snapped to 7, and a fixed-size sampler refuses it.
  p <- inclusion_prob(1:10, n = 6.99995)
  expect_equal(sum(p), 6.99995)
  expect_error(
    unequal_prob_wor(p, method = "brewer"),
    "sum\\(pik\\) = 6.99995 is not close to an integer"
  )
})

test_that("inclusion_prob handles the ends of its domain", {
  expect_identical(inclusion_prob(c(0, 0, 0), n = 0), c(0, 0, 0))
  expect_identical(inclusion_prob(numeric(0), n = 0), numeric(0))
  expect_equal(inclusion_prob(c(0, 3, 0, 5, 5), n = 3), c(0, 1, 0, 1, 1))
  x <- c(3, 7, 1, 20, 9)
  expect_equal(inclusion_prob(x * 1e-200, n = 2.6), inclusion_prob(x, 2.6))
  expect_equal(inclusion_prob(x * 1e200, n = 2.6), inclusion_prob(x, 2.6))
})

test_that("a target above the positive sizes is refused", {
  msg <- "exceeds the number of units with positive size"
  expect_error(inclusion_prob(c(0, 0, 0, 5, 5), n = 2.4), paste(msg, "\\(2\\)"))
  expect_error(
    inclusion_prob(c(0, 0, 0, 5, 5), n = 2 + 1e-9),
    paste(msg, "\\(2\\)")
  )
  expect_error(inclusion_prob(c(0, 0, 0), n = 0.3), paste(msg, "\\(0\\)"))
  expect_error(inclusion_prob(c(0, 0, 0), n = 1e-12), paste(msg, "\\(0\\)"))
  # The message reports n as given.
  expect_error(inclusion_prob(c(0, 0, 0, 5, 5), n = 2.6), "'n' \\(2.6\\)")
  expect_error(inclusion_prob(numeric(0), n = 1), "cannot exceed length")
})

test_that("inclusion_prob survives extreme size ratios", {
  # Once the large units are certain, the rest are subnormal or 0 against
  # the largest size. A zero-size unit used to get 0 * Inf = NaN, and the
  # small units were all capped or all lost.
  expect_identical(inclusion_prob(c(0, 1e-300, 1e10), n = 2), c(0, 1, 1))
  expect_equal(inclusion_prob(c(0, 1e-300, 1e10), n = 1.6), c(0, 0.6, 1))
  expect_equal(inclusion_prob(c(1e-300, 1e-300, 1e10), n = 2), c(0.5, 0.5, 1))
  expect_equal(inclusion_prob(c(1e-300, 3e-300, 1e10), n = 2), c(0.25, 0.75, 1))
  expect_equal(
    inclusion_prob(c(1e-200, 1e-200, 1e200), n = 1.5),
    c(0.25, 0.25, 1)
  )
  expect_equal(
    inclusion_prob(c(1e-300, 1e-300, 1e-300, 1e10, 1e10), n = 2.7),
    c(rep(0.7 / 3, 3), 1, 1)
  )
  # A size that underflows to 0 against the largest is still positive.
  expect_equal(inclusion_prob(c(1e-320, 1e10), n = 1.5), c(0.5, 1))
  expect_equal(inclusion_prob(c(1e-320, 1e10), n = 2), c(1, 1))
  expect_equal(
    inclusion_prob(c(1e-320, 3e-320, 1e10, 1e300), n = 2.5),
    c(0.125, 0.375, 1, 1)
  )
  # The message reports n to 15 digits.
  expect_error(inclusion_prob(c(0, 5, 5), n = 2 + 1e-9), "'n' \\(2.000000001\\)")
})

test_that("inclusion_prob matches an exact solution across extreme ratios", {
  # The k largest units are certain and lambda = (n - k) / sum(rest),
  # solved in log space so that no size ratio overflows.
  ref_exact <- function(x, n) {
    out <- numeric(length(x))
    pos <- which(x > 0)
    o <- pos[order(x[pos], decreasing = TRUE)]
    lx <- log(x[o])
    P <- length(o)
    for (k in 0:(P - 1)) {
      if (n - k <= 0) break
      rest <- lx[(k + 1):P]
      ll <- log(n - k) - (max(rest) + log(sum(exp(rest - max(rest)))))
      if (ll + lx[k + 1] <= 0) {
        out[o] <- c(rep(1, k), exp(ll + rest))
        return(out)
      }
    }
    out[o] <- 1
    out
  }
  set.seed(20261011)
  for (k in 1:300) {
    N <- sample(2:30, 1)
    x <- 10^runif(N, sample(c(-320, -250, -100), 1), 300)
    x[sample(N, sample(0:(N - 1), 1))] <- 0
    P <- sum(x > 0)
    n <- switch(sample(3, 1), runif(1, 0, P), P, P - runif(1) * 1e-9)
    expect_lt(max(abs(inclusion_prob(x, n) - ref_exact(x, n))), 1e-10)
  }
})

test_that("inclusion_prob matches a bisection reference on random inputs", {
  ref_pik <- function(x, n) {
    f <- function(l) sum(pmin(1, l * x)) - n
    hi <- 2 / min(x[x > 0])
    pmin(1, stats::uniroot(f, c(0, hi), tol = 1e-15)$root * x)
  }
  set.seed(20261010)
  for (k in 1:300) {
    N <- sample(2:40, 1)
    x <- rexp(N)^sample(c(1, 4), 1)
    x[sample(N, sample(0:(N - 1), 1))] <- 0
    n <- runif(1, 0.01, sum(x > 0) - 0.01)
    expect_equal(inclusion_prob(x, n), ref_pik(x, n), tolerance = 1e-8)
  }
})

test_that("fixed-size samplers refuse a fractional sum", {
  p <- inclusion_prob(1:10, n = 2.6)
  for (m in c("brewer", "sampford", "cps", "systematic")) {
    expect_error(unequal_prob_wor(p, method = m), "not close to an integer")
  }
  s <- unequal_prob_wor(p, method = "poisson")
  expect_equal(s$pik, p)
  expect_equal(s$n, 2.6)
})

test_that("n near an integer is reported with six digits", {
  expect_error(equal_prob_wor(10, 7.0002), "n \\(7.0002\\)")
})

test_that("inclusion_prob rejects Inf n", {
  expect_error(inclusion_prob(1:10, n = Inf), "finite")
})

test_that("inclusion_prob rejects -Inf n", {
  expect_error(inclusion_prob(1:10, n = -Inf), "finite")
})

test_that("inclusion_prob rejects NaN n", {
  expect_error(inclusion_prob(1:10, n = NaN), "NA")
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

test_that("expected_hits errors when n is missing", {
  expect_error(expected_hits(c(10, 20, 30)), "required")
})

test_that("expected_hits errors when n is not a single numeric", {
  expect_error(expected_hits(c(10, 20), n = "a"), "single numeric")
  expect_error(expected_hits(c(10, 20), n = c(1, 2)), "single numeric")
})

test_that("expected_hits errors when n is NA or negative", {
  expect_error(expected_hits(c(10, 20), n = NA_real_), "NA")
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

test_that("expected_hits rejects Inf n", {
  expect_error(expected_hits(c(1, 2, 3), n = Inf), "finite")
})

test_that("expected_hits rejects -Inf n", {
  expect_error(expected_hits(c(1, 2, 3), n = -Inf), "finite")
})

test_that("expected_hits rejects NaN n", {
  expect_error(expected_hits(c(1, 2, 3), n = NaN), "NA")
})

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

# Behavior lock tests (pre-flight for issue #5 dead-branch removal).
# These exercise the C paths that survived, and also assert the R-layer
# refuses non-finite inputs so the C-level NA handling is provably
# unreachable via the public API.

test_that("inclusion_prob handles all-zero input without NaN", {
  expect_equal(inclusion_prob(c(0, 0, 0, 0), n = 0), c(0, 0, 0, 0))
  expect_error(
    inclusion_prob(c(0, 0, 0), n = 1),
    "exceeds the number of units with positive size"
  )
})

test_that("inclusion_prob mixed zero and positive", {
  pik <- inclusion_prob(c(0, 10, 20, 30), n = 2)
  expect_equal(pik[1], 0)
  expect_equal(sum(pik), 2)
})

test_that("inclusion_prob iterative capping with several certainty selections", {
  # Multiple capping rounds, the rescale loop must handle the cascade.
  pik <- inclusion_prob(c(1, 1, 1, 1000, 1000, 1000), n = 4)
  expect_equal(sum(pik[4:6]), 3)
  expect_equal(sum(pik), 4)
  expect_equal(sum(pik[1:3]), 1)
})

test_that("inclusion_prob rejects NA, NaN, Inf at the R layer", {
  # These lock in the R-layer gate that keeps non-finite out of C.
  expect_error(inclusion_prob(c(NA, 1, 2, 3), n = 2), "missing values")
  expect_error(inclusion_prob(c(NaN, 1, 2, 3), n = 2), "missing values|finite")
  expect_error(inclusion_prob(c(Inf, 1, 2, 3), n = 2), "finite")
  expect_error(inclusion_prob(c(-Inf, 1, 2, 3), n = 2), "finite")
})
