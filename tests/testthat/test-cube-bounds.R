# Cube method with inequality constraints (Tripet & Tillé 2026)

# Indicator matrix for the levels of a grouping vector
.ind <- function(f) {
  sapply(unique(f), function(g) as.double(f == g))
}

test_that("bounds must be a complete, named list", {
  pik <- rep(0.5, 4)
  B <- matrix(1, 4, 1)
  expect_error(
    balanced_wor(pik, bounds = TRUE),
    "must be a list"
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = 0)),
    "missing element"
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = 0, upper = 3, extra = 1)),
    "unknown element"
  )
})

test_that("bounds$B is validated", {
  pik <- rep(0.5, 4)
  expect_error(
    balanced_wor(
      pik,
      bounds = list(B = matrix(1, 3, 1), lower = 0, upper = 3)
    ),
    "does not match length\\(pik\\)"
  )
  expect_error(
    balanced_wor(
      pik,
      bounds = list(B = matrix(NA_real_, 4, 1), lower = 0, upper = 3)
    ),
    "NA, NaN, or Inf"
  )
  expect_error(
    balanced_wor(
      pik,
      bounds = list(
        B = matrix(numeric(0), 4, 0),
        lower = double(0),
        upper = double(0)
      )
    ),
    "at least one column"
  )
})

test_that("bound vectors are validated", {
  pik <- rep(0.5, 4)
  B <- matrix(1, 4, 1)
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = c(0, 0), upper = 3)),
    "length ncol"
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = NA_real_, upper = 3)),
    "missing values"
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = 3, upper = 1)),
    "'bounds\\$lower' must be <="
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = -Inf, upper = Inf)),
    "no finite bound"
  )
})

test_that("infeasible starting probabilities error", {
  pik <- rep(0.5, 4) # colSums(B * pik) = 2
  B <- matrix(1, 4, 1)
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = 3, upper = 4)),
    "violate 'bounds'"
  )
  expect_error(
    balanced_wor(pik, bounds = list(B = B, lower = 0, upper = 1)),
    "violate 'bounds'"
  )
})

test_that("bounds cannot be combined with strata", {
  pik <- rep(0.5, 8)
  B <- matrix(1, 8, 1)
  expect_error(
    balanced_wor(
      pik,
      strata = rep(1:2, each = 4),
      bounds = list(B = B, lower = 3, upper = 5)
    ),
    "cannot be combined with 'strata'"
  )
})

test_that("bounds are rejected for registered methods", {
  register_method(
    "bounds_toy",
    type = "balanced",
    sample_fn = function(pik, n = NULL, aux = NULL, ...) {
      sample.int(length(pik), size = n)
    }
  )
  on.exit(unregister_method("bounds_toy"), add = TRUE)
  pik <- rep(0.5, 4)
  expect_error(
    balanced_wor(
      pik,
      bounds = list(B = matrix(1, 4, 1), lower = 1, upper = 3),
      method = "bounds_toy"
    ),
    "registered balanced methods do not support 'bounds'"
  )
})

test_that("bounds are stored on the design object", {
  pik <- rep(0.5, 8)
  bounds <- list(
    B = .ind(rep(1:2, each = 4)),
    lower = c(1, 1),
    upper = c(3, 3)
  )
  s <- balanced_wor(pik, bounds = bounds)
  expect_s3_class(s, "balanced")
  expect_identical(s$bounds, bounds)
  expect_true(s$fixed_size)
  expect_equal(s$n, 4L)
  expect_equal(s$method, "cube")

  s0 <- balanced_wor(pik)
  expect_null(s0$bounds)
})

test_that("bounds draws are reproducible with set.seed", {
  pik <- rep(0.5, 8)
  bounds <- list(B = .ind(rep(1:2, each = 4)), lower = c(2, 2), upper = c(2, 2))
  set.seed(42)
  s1 <- balanced_wor(pik, bounds = bounds)$sample
  set.seed(42)
  s2 <- balanced_wor(pik, bounds = bounds)$sample
  expect_identical(s1, s2)
})

# Category bounding (controlled selection)
test_that("single-factor category counts stay within floor/ceil", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  groups <- rep(c("a", "b", "c"), each = 4)
  B <- .ind(groups)
  S <- colSums(B * pik) # 1.9, 2.1, 2.0

  set.seed(1)
  for (i in 1:200) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_length(s$sample, 6)
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
  }
})

test_that("two-way margins stay within floor/ceil (Goodman-Kish control)", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  rowf <- rep(c("a", "b", "c"), each = 4)
  colf <- rep(c("x", "y"), times = 6)
  B <- cbind(.ind(rowf), .ind(colf))
  S <- colSums(B * pik)

  set.seed(2)
  for (i in 1:200) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_length(s$sample, 6)
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
  }
})

test_that("equality bounds (lower == upper) give exact counts", {
  pik <- rep(0.5, 12)
  groups <- rep(c("a", "b", "c"), each = 4)
  B <- .ind(groups)

  set.seed(3)
  for (i in 1:100) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = c(2, 2, 2), upper = c(2, 2, 2))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_equal(unname(cnt), c(2, 2, 2))
  }
})

test_that("one-sided bounds work (-Inf / Inf sides)", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  groups <- rep(c("a", "b", "c"), each = 4)
  B <- .ind(groups)
  S <- colSums(B * pik)

  set.seed(4)
  # Minimum group sizes (paper section 5): count >= floor(S)
  for (i in 1:100) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = rep(Inf, 3))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_true(all(cnt >= floor(S)))
  }
  # Maximum group sizes: count <= ceiling(S)
  for (i in 1:100) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = rep(-Inf, 3), upper = ceiling(S))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_true(all(cnt <= ceiling(S)))
  }
})

test_that("overlapping (non-nested) group bounds are all respected", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  g1 <- rep(c("a", "b", "c"), each = 4) # blocks
  g2 <- rep(c("u", "v", "w"), times = 4) # interleaved
  B <- cbind(.ind(g1), .ind(g2))
  S <- colSums(B * pik)

  set.seed(5)
  for (i in 1:200) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
  }
})

test_that("bounds combine with auxiliary balancing", {
  set.seed(6)
  N <- 24
  pik <- rep(0.25, N)
  x <- matrix(as.double(1:N) + rnorm(N), ncol = 1)
  groups <- rep(c("a", "b", "c"), each = 8)
  B <- .ind(groups)
  S <- colSums(B * pik) # 2, 2, 2

  tx <- sum(x)
  for (i in 1:100) {
    s <- balanced_wor(
      pik,
      aux = x,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
    # HT estimate of the x total stays in a sane band around the truth
    ht <- sum(x[s$sample, ] / pik[s$sample])
    expect_lt(abs(ht - tx) / tx, 0.5)
  }
})

test_that("non-indicator constraint columns work (HT estimator bound)", {
  # Continuous-valued constraints that end the flight tight against a
  # boundary generally admit no exact 0/1 completion (the boundary
  # value is not integer-representable), so the landing may relax them
  # with a warning; the paper's own supplementary code inflates the
  # bounds in this situation. The contract: bounds hold whenever no
  # relaxation warning was raised, and E(s) = pik holds regardless.
  set.seed(7)
  N <- 20
  pik <- rep(0.3, N)
  x <- as.double(1:N)
  tx <- sum(x)
  # Bound the HT estimator of the x total within +/- 15% of the truth
  B <- matrix(x / pik, ncol = 1)
  n_clean <- 0
  for (i in 1:100) {
    relaxed <- FALSE
    s <- withCallingHandlers(
      balanced_wor(
        pik,
        bounds = list(B = B, lower = 0.85 * tx, upper = 1.15 * tx)
      ),
      warning = function(w) {
        relaxed <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    if (!relaxed) {
      n_clean <- n_clean + 1
      ht <- sum(x[s$sample] / pik[s$sample])
      expect_true(ht >= 0.85 * tx && ht <= 1.15 * tx)
    }
  }
  expect_gt(n_clean, 0)
})

# Paper applications
test_that("controlled matrix rounding: Cochran example (paper section 6)", {
  # Cochran (1977, p. 124): 165 schools, 5 x 4 table, divided by 16.5
  # with integer parts removed. sum(pik) = 165/16.5 - 2 = 8 exactly.
  counts <- matrix(
    c(15, 21, 17, 9, 10, 8, 13, 7, 6, 9, 5, 8, 4, 3, 6, 6, 3, 2, 5, 8),
    nrow = 5,
    byrow = TRUE
  )
  M <- counts / 16.5
  M <- M - floor(M)
  pik <- as.vector(M) # column-major cells
  rowi <- rep(1:5, times = 4)
  coli <- rep(1:4, each = 5)
  B <- cbind(.ind(rowi), .ind(coli))
  S <- colSums(B * pik)

  expect_equal(sum(pik), 8)

  set.seed(8)
  for (i in 1:200) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    expect_length(s$sample, 8)
    cnt <- colSums(B[s$sample, , drop = FALSE])
    # Every row and column margin rounds to an adjacent integer
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
  }
})

test_that("MU284 category bounding (paper example 1)", {
  skip_if_not_installed("sampling")
  data("MU284", package = "sampling", envir = environment())

  # Cross REG (8 regions) with P75 quartile classes (4), n = 50
  q <- quantile(MU284$P75, (0:4) / 4, type = 3)
  q[1] <- 0
  q[5] <- 700
  z <- cut(MU284$P75, q)
  N <- nrow(MU284)
  pik <- rep(50 / N, N)
  B <- cbind(.ind(MU284$REG), .ind(as.integer(z)))
  S <- colSums(B * pik)

  set.seed(9)
  for (i in 1:20) {
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = floor(S), upper = ceiling(S))
    )
    expect_length(s$sample, 50)
    cnt <- colSums(B[s$sample, , drop = FALSE])
    expect_true(all(cnt >= floor(S) & cnt <= ceiling(S)))
  }
})

test_that("inclusion probabilities are respected under bounds", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  rowf <- rep(c("a", "b", "c"), each = 4)
  colf <- rep(c("x", "y"), times = 6)
  B <- cbind(.ind(rowf), .ind(colf))
  S <- colSums(B * pik)

  set.seed(10)
  M <- 4000
  s <- balanced_wor(
    pik,
    bounds = list(B = B, lower = floor(S), upper = ceiling(S)),
    nrep = M
  )
  emp <- tabulate(as.vector(s$sample), nbins = 12) / M
  se <- sqrt(pik * (1 - pik) / M)
  expect_true(all(abs(emp - pik) < 5 * se))
})

test_that("infeasible integer structures relax with a warning, unbiasedly", {
  # Groups {1,2}, {2,3}, {1,3} with exact count 1 each imply
  # s1 + s2 + s3 = 1.5: no integer solution exists alongside n = 2.
  pik <- rep(0.5, 4)
  B <- cbind(
    c(1, 1, 0, 0),
    c(0, 1, 1, 0),
    c(1, 0, 1, 0)
  )
  set.seed(11)
  expect_warning(
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = c(1, 1, 1), upper = c(1, 1, 1))
    ),
    "relaxed during landing"
  )
  expect_length(s$sample, 2)

  # E(s) = pik still holds across replicates
  M <- 3000
  hits <- integer(4)
  for (i in 1:M) {
    suppressWarnings(
      s <- balanced_wor(
        pik,
        bounds = list(B = B, lower = c(1, 1, 1), upper = c(1, 1, 1))
      )
    )
    hits[s$sample] <- hits[s$sample] + 1L
  }
  emp <- hits / M
  se <- sqrt(pik * (1 - pik) / M)
  expect_true(all(abs(emp - pik) < 5 * se))
})

test_that("nrep > 1 respects bounds in every replicate", {
  pik <- c(0.3, 0.5, 0.7, 0.4, 0.6, 0.5, 0.2, 0.8, 0.5, 0.45, 0.55, 0.5)
  groups <- rep(c("a", "b", "c"), each = 4)
  B <- .ind(groups)
  S <- colSums(B * pik)
  bounds <- list(B = B, lower = floor(S), upper = ceiling(S))

  set.seed(12)
  s <- balanced_wor(pik, bounds = bounds, nrep = 500)
  expect_true(is.matrix(s$sample))
  expect_equal(dim(s$sample), c(6L, 500L))
  expect_identical(s$bounds, bounds)
  cnt <- apply(s$sample, 2, function(ix) colSums(B[ix, , drop = FALSE]))
  expect_true(all(cnt >= floor(S)) && all(cnt <= ceiling(S)))
})

test_that("batch relaxation warning reports affected replicates", {
  pik <- rep(0.5, 4)
  B <- cbind(
    c(1, 1, 0, 0),
    c(0, 1, 1, 0),
    c(1, 0, 1, 0)
  )
  set.seed(13)
  expect_warning(
    s <- balanced_wor(
      pik,
      bounds = list(B = B, lower = c(1, 1, 1), upper = c(1, 1, 1)),
      nrep = 20
    ),
    "relaxed during landing in \\d+ of 20 replicate"
  )
  expect_equal(dim(s$sample), c(2L, 20L))
})

test_that("no relaxation attribute leaks into the design object", {
  pik <- rep(0.5, 4)
  B <- cbind(
    c(1, 1, 0, 0),
    c(0, 1, 1, 0),
    c(1, 0, 1, 0)
  )
  set.seed(14)
  suppressWarnings({
    s1 <- balanced_wor(
      pik,
      bounds = list(B = B, lower = c(1, 1, 1), upper = c(1, 1, 1))
    )
    s2 <- balanced_wor(
      pik,
      bounds = list(B = B, lower = c(1, 1, 1), upper = c(1, 1, 1)),
      nrep = 5
    )
  })
  expect_null(attr(s1$sample, "relaxed", exact = TRUE))
  expect_null(attr(s2$sample, "relaxed_reps", exact = TRUE))
})
