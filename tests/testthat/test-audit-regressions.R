# Regression tests for the 2026-07-10 algorithm audit
# (dev/ALGORITHM_AUDIT.md, fix plan dev/ALGORITHM_FIX_PLAN.md).
#
# These encode the post-fix contracts:
#   - eps never modifies the design (only exact 0/1 are special),
#   - fixed-size methods always return exactly the declared size,
#   - CPS is unbiased for large ordinary designs,
#   - CPS joint probabilities are finite and satisfy design identities,
#   - size conversion is scale invariant,
#   - cube auxiliary conditioning detects rank.

## Audit #2: eps must not silently change the design

test_that("eps is rejected by fixed-size unequal probability methods", {
  p <- rep(0.1, 10)
  for (m in c("cps", "sampford", "brewer", "systematic", "sps", "pareto")) {
    expect_error(
      unequal_prob_wor(p, m, eps = 0.2),
      "no longer modifies the design",
      info = m
    )
  }
})

test_that("balanced_wor rejects pik that eps would reclassify", {
  expect_error(
    balanced_wor(rep(0.1, 10), eps = 0.2),
    "exactly 0 or exactly 1"
  )
  # default eps = 1e-10 does not reclassify ordinary pik
  s <- balanced_wor(rep(0.1, 10))
  expect_identical(length(s$sample), s$n)
})

test_that("fixed-size validation rejects non-integer pik sums strictly", {
  # off by 1e-3: statistically impossible fixed-size design
  p <- rep(0.1, 10)
  p[1] <- p[1] + 1e-3
  for (m in c("cps", "sampford", "brewer", "systematic", "sps", "pareto")) {
    expect_error(
      unequal_prob_wor(p, m),
      "not close to an integer",
      info = m
    )
  }
  # off by 1e-5: previously accepted by tol = 1e-4, now rejected
  p <- rep(0.1, 10)
  p[1] <- p[1] + 1e-5
  expect_error(unequal_prob_wor(p, "systematic"), "not close to an integer")
})

test_that("fixed-size validation accepts floating-point-level sum error", {
  # inclusion_prob output at large N carries ~1e-9 accumulation error
  set.seed(42)
  x <- rlnorm(1e5, 5, 2)
  pik <- inclusion_prob(x, 500)
  s <- unequal_prob_wor(pik, "systematic")
  expect_identical(length(s$sample), 500L)
})

## Audit #1 + #2: fixed-size methods return exactly n

test_that("fixed-size WOR methods return exactly n with exact 0/1 units", {
  set.seed(101)
  # exact zeros and ones mixed with interior values summing to an integer
  base <- c(0, 1, 1, 0.3, 0.7, 0.25, 0.75, 0.5, 0.5, 0)
  n_declared <- as.integer(round(sum(base)))
  for (m in c("cps", "sampford", "brewer", "systematic", "sps", "pareto")) {
    for (rep in 1:20) {
      s <- suppressWarnings(unequal_prob_wor(base, m))
      expect_identical(length(s$sample), n_declared, info = m)
      expect_true(all(s$sample >= 1L & s$sample <= 10L), info = m)
      expect_false(anyDuplicated(s$sample) > 0, info = m)
      # exact-1 units always selected, exact-0 never
      expect_true(all(c(2L, 3L) %in% s$sample), info = m)
      expect_false(any(c(1L, 10L) %in% s$sample), info = m)
    }
  }
  s <- balanced_wor(base)
  expect_identical(length(s$sample), n_declared)
})

test_that("systematic PPS batch never writes out of bounds (audit crash)", {
  # The original crash used eps = .1 to create a non-integer residual.
  # That input is now rejected; the sampler itself must return exactly n
  # for every replicate of a valid design.
  set.seed(4)
  p <- c(0.95, 0.04, 0.20, 0.81)  # sums to 2 exactly
  s <- unequal_prob_wor(p, "systematic", nrep = 5000)
  expect_identical(dim(s$sample), c(2L, 5000L))
  expect_true(all(s$sample >= 1L & s$sample <= 4L))
})

test_that("balanced batch size mismatches fail before writing past a column", {
  pik <- rep(5 / 11, 11)
  spread <- matrix(as.double(seq_along(pik)), ncol = 1L)
  aux <- spread

  set.seed(1)
  expect_error(
    balanced_wor(pik, spread = spread, method = "lpm2", nrep = 3, eps = 0.45),
    "lpm2 batch draw 2 produced size 6, expected 5"
  )

  # With this deliberately extreme eps, Windows can hit SCPS's defensive
  # feasibility check before extraction; other platforms reach the bounded
  # extractor and report the intended size mismatch. Both paths fail safely.
  set.seed(1)
  expect_error(
    balanced_wor(pik, spread = spread, method = "scps", nrep = 3, eps = 0.45),
    paste0(
      "SCPS (batch draw 1 produced size 4, expected 5|",
      "maximal weights are numerically infeasible in draw 1)"
    )
  )

  set.seed(1)
  expect_error(
    balanced_wor(pik, aux = aux, method = "cube", nrep = 3, eps = 0.45),
    "cube batch draw 1 produced size 7, expected 5"
  )

  set.seed(1)
  expect_error(
    balanced_wor(
      pik,
      aux = aux,
      strata = rep(1L, length(pik)),
      method = "cube",
      nrep = 3,
      eps = 0.45
    ),
    "stratified cube batch draw 1 produced size 7, expected 5"
  )
})

test_that("CPS zero-work design initializes its weight buffer", {
  pik <- 1 - c(1, 2, 3, 4) * 1e-13
  design <- .Call(sondage:::C_cps_design, as.double(pik))

  expect_identical(design$n_work, 0L)
  expect_identical(design$w, numeric(length(pik)))

  sample <- unequal_prob_wor(pik, method = "cps", nrep = 2)
  expect_identical(dim(sample$sample), c(4L, 2L))
})

test_that("near-boundary interior pik are sampled, not snapped", {
  # pik close to (but not exactly) 0/1 stay in the design
  p <- c(1e-9, 1 - 1e-9, 0.5, 0.5)  # sums to 2 within fp tolerance
  s <- suppressWarnings(unequal_prob_wor(p, "systematic"))
  expect_identical(length(s$sample), 2L)
  # unit 2 selected with prob 1 - 1e-9: must be present in practice
  expect_true(2L %in% s$sample)
})

## Audit #3: CPS unbiased for large ordinary designs

test_that("CPS equal-probability large design has no positional bias", {
  # audit repro: N = 2000, pik = 0.5 gave first half 0.235, second half 0.765
  set.seed(2026)
  N <- 2000L
  z <- unequal_prob_wor(rep(0.5, N), "cps", nrep = 200)$sample
  counts <- tabulate(as.vector(z), N) / ncol(z)
  expect_lt(abs(mean(counts[1:(N / 2)]) - 0.5), 0.02)
  expect_lt(abs(mean(counts[(N / 2 + 1):N]) - 0.5), 0.02)
  expect_gt(min(counts), 0.25)
  expect_lt(max(counts), 0.75)
})

test_that("CPS unequal large design achieves target marginals", {
  # unequal pik so the conditional-table path (not an SRS fast path) is used
  set.seed(7)
  N <- 800L
  x <- exp(rnorm(N, 0, 0.8))
  pik <- inclusion_prob(x, 300)
  nrep <- 400L
  z <- unequal_prob_wor(pik, "cps", nrep = nrep)$sample
  emp <- tabulate(as.vector(z), N) / nrep
  se <- sqrt(pik * (1 - pik) / nrep)
  zscore <- (emp - pik) / pmax(se, 1e-6)
  # no unit may deviate wildly; average standardized error near 0
  expect_lt(max(abs(zscore)), 6)
  expect_lt(abs(mean(zscore)), 0.5)
})

test_that("CPS single draws also avoid positional bias at large N", {
  set.seed(11)
  N <- 1500L
  counts <- integer(N)
  nrep <- 60L
  for (i in seq_len(nrep)) {
    idx <- unequal_prob_wor(rep(0.5, N), "cps")$sample
    counts[idx] <- counts[idx] + 1L
  }
  # crude but catches the audit failure mode (halves at 0.235 vs 0.765)
  expect_lt(
    abs(mean(counts[1:(N / 2)]) - mean(counts[(N / 2 + 1):N])) / nrep,
    0.15
  )
})

## Audit #4: CPS joint inclusion probabilities

test_that("CPS joint probabilities are stable for repeated groups (N = 200)", {
  N <- 200L
  p <- c(rep(0.2, N / 2), rep(0.8, N / 2))
  s <- unequal_prob_wor(p, "cps")
  J <- joint_inclusion_prob(s)
  n <- sum(p)

  expect_true(all(is.finite(J)))
  expect_identical(J, t(J))
  expect_equal(diag(J), p, tolerance = 1e-12)
  # fixed-size row identity: sum_j pi_ij = n * pi_i
  expect_lt(max(abs(rowSums(J) - n * p)), 1e-6)
  # off-diagonal joint < both marginals (strict for non-degenerate design)
  offdiag <- J[upper.tri(J)]
  expect_true(all(offdiag > 0))
  expect_lt(J[N - 1L, N], 0.8 - 0.05)
})

test_that("CPS joint probabilities contain no NaN at N = 800", {
  N <- 800L
  p <- c(rep(0.2, N / 2), rep(0.8, N / 2))
  s <- unequal_prob_wor(p, "cps")
  J <- joint_inclusion_prob(s)
  expect_true(all(is.finite(J)))
  expect_lt(max(abs(rowSums(J) - sum(p) * p)), 1e-6)
})

test_that("CPS joint probabilities match Monte Carlo co-occurrence", {
  set.seed(99)
  N <- 10L
  pik <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.6, 0.6, 0.5, 0.5)
  s <- unequal_prob_wor(pik, "cps")
  J <- joint_inclusion_prob(s)

  nrep <- 40000L
  z <- unequal_prob_wor(pik, "cps", nrep = nrep)$sample
  co <- matrix(0, N, N)
  for (r in seq_len(ncol(z))) {
    idx <- z[, r]
    co[idx, idx] <- co[idx, idx] + 1
  }
  co <- co / nrep
  # 5 sigma on each pairwise estimate
  tol <- 5 * sqrt(pmax(J * (1 - J), 0.05) / nrep)
  expect_true(all(abs(co - J)[upper.tri(J)] <= tol[upper.tri(J)]))
})

test_that("CPS full and sampled_only joint probabilities agree", {
  set.seed(5)
  N <- 60L
  pik <- inclusion_prob(exp(rnorm(N)), 20)
  s <- unequal_prob_wor(pik, "cps")
  J_full <- joint_inclusion_prob(s)
  J_sub <- joint_inclusion_prob(s, sampled_only = TRUE)
  idx <- s$sample
  expect_equal(J_sub, J_full[idx, idx, drop = FALSE], tolerance = 1e-10,
               ignore_attr = TRUE)
})

## Audit #5: scale-invariant size conversion

test_that("inclusion_prob is invariant to positive rescaling", {
  x <- c(1, 2, 3, 4, 10)
  base <- inclusion_prob(x, 2)
  for (c_scale in c(1e-300, 1e-30, 1, 1e30, 1e300)) {
    expect_equal(inclusion_prob(x * c_scale, 2), base, tolerance = 1e-12,
                 info = paste("scale", c_scale))
  }
})

test_that("inclusion_prob handles extreme magnitudes (audit repro)", {
  x <- rep(.Machine$double.xmax, 2)
  expect_equal(inclusion_prob(x, 1), c(0.5, 0.5))
  x <- rep(.Machine$double.xmin * .Machine$double.eps, 2)
  expect_equal(inclusion_prob(x, 1), c(0.5, 0.5))
})

test_that("expected_hits is invariant to positive rescaling", {
  x <- c(1, 2, 3, 4)
  base <- expected_hits(x, 3)
  expect_equal(sum(base), 3)
  for (c_scale in c(1e-300, 1e300)) {
    expect_equal(expected_hits(x * c_scale, 3), base, tolerance = 1e-12)
  }
  x <- rep(.Machine$double.xmax, 2)
  expect_equal(expected_hits(x, 1), c(0.5, 0.5))
})

## Audit #6: cube auxiliary conditioning rank detection

test_that("cube conditioning removes exactly dependent columns", {
  x <- as.double(1:10)
  pik <- rep(0.2, 10)

  # scalar multiples and duplicates reduce to one column
  z <- sondage:::.condition_cube_aux(cbind(x, 2 * x, x), pik,
                                     qr_tol = 1e-8)
  expect_identical(ncol(z), 1L)
  expect_identical(qr(z)$rank, 1L)

  # independent columns are kept
  z2 <- sondage:::.condition_cube_aux(cbind(x, x^2), pik, qr_tol = 1e-8)
  expect_identical(ncol(z2), 2L)

  # constant column is dropped (zero weighted variance after centering)
  z3 <- sondage:::.condition_cube_aux(cbind(rep(1, 10), x), pik,
                                      qr_tol = 1e-8)
  expect_identical(ncol(z3), 1L)
})

test_that("cube conditioning respects qr_tol for near-dependence", {
  set.seed(3)
  x <- as.double(1:20)
  pik <- rep(0.25, 20)
  noise <- rnorm(20)
  # second column = first + tiny noise: dependent at loose tol,
  # independent at strict tol
  aux <- cbind(x, x + 1e-10 * noise)
  z_loose <- sondage:::.condition_cube_aux(aux, pik, qr_tol = 1e-4)
  z_strict <- sondage:::.condition_cube_aux(aux, pik, qr_tol = 1e-14)
  expect_identical(ncol(z_loose), 1L)
  expect_identical(ncol(z_strict), 2L)
})

## Audit #7: systematic joint probabilities (reference check)

test_that("systematic JIP matches direct interval reference", {
  # R reference: pi_ij = length of overlap of the two circular arcs
  sys_jip_ref <- function(pik) {
    N <- length(pik)
    V <- cumsum(pik)
    starts <- c(0, V[-N]) %% 1
    J <- outer(pik, pik)
    diag(J) <- pik
    overlap <- function(a1, l1, a2, l2) {
      # overlap of circular arcs [a1, a1+l1), [a2, a2+l2) on [0,1):
      # sum the linear overlaps over integer shifts of the second arc
      tot <- 0
      for (k in -1:1) {
        lo <- max(a1, a2 + k)
        hi <- min(a1 + l1, a2 + k + l2)
        tot <- tot + max(0, hi - lo)
      }
      tot
    }
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        J[i, j] <- J[j, i] <- overlap(starts[i], pik[i], starts[j], pik[j])
      }
    }
    J
  }

  set.seed(21)
  for (trial in 1:5) {
    N <- 30L
    pik <- inclusion_prob(exp(rnorm(N)), 8)
    s <- unequal_prob_wor(pik, "systematic")
    J <- joint_inclusion_prob(s)
    expect_equal(unclass(J), sys_jip_ref(pik), tolerance = 1e-10,
                 ignore_attr = TRUE)
    expect_lt(max(abs(rowSums(J) - sum(pik) * pik)), 1e-8)
  }
})

test_that("systematic JIP preserves exact structural zeros", {
  s <- equal_prob_wor(10, 3, method = "systematic")
  J <- joint_inclusion_prob(s)

  expect_identical(J[8, 9], 0)
  expect_identical(J[9, 10], 0)

  sampled <- s
  sampled$sample <- c(8L, 9L)
  J_sub <- joint_inclusion_prob(sampled, sampled_only = TRUE)
  expect_identical(J_sub[1, 2], 0)

  expect_warning(
    syg <- sampling_cov(s, weighted = TRUE),
    "joint inclusion probabilities are zero"
  )
  expect_true(is.na(syg[8, 9]))
  expect_true(is.na(syg[9, 10]))
})

test_that("equal-probability batch draws normalize near-integer n", {
  near_three <- 3 + 5e-5

  for (method in c("srs", "systematic")) {
    one <- equal_prob_wor(10, near_three, method = method)
    many <- equal_prob_wor(10, near_three, method = method, nrep = 2)

    expect_identical(one$n, 3L, info = method)
    expect_identical(many$n, 3L, info = method)
    expect_identical(dim(many$sample), c(3L, 2L), info = method)
  }
})

test_that("systematic PPS returns sorted indices with certainty units", {
  set.seed(17)
  for (i in seq_len(20)) {
    s <- unequal_prob_wor(c(0.5, 1, 0.5), method = "systematic")
    expect_false(is.unsorted(s$sample))
  }
})

test_that("optimized batch paths preserve sequential single-draw streams", {
  compare_fixed <- function(batch, draw_one) {
    singles <- lapply(seq_len(ncol(batch)), function(i) draw_one())
    expect_identical(batch, do.call(cbind, singles))
  }

  for (method in c("srs", "systematic")) {
    set.seed(12)
    batch <- equal_prob_wor(20, 5, method, nrep = 4)$sample
    set.seed(12)
    compare_fixed(batch, function() equal_prob_wor(20, 5, method)$sample)
  }

  set.seed(12)
  batch <- equal_prob_wr(20, 5, nrep = 4)
  set.seed(12)
  singles <- lapply(seq_len(4), function(i) equal_prob_wr(20, 5))
  expect_identical(batch$sample, do.call(cbind, lapply(singles, `[[`, "sample")))
  expect_identical(batch$hits, do.call(cbind, lapply(singles, `[[`, "hits")))

  hits <- c(0.2, 0.5, 0.8, 0.5)
  for (method in c("chromy", "multinomial")) {
    set.seed(12)
    batch <- unequal_prob_wr(hits, method, nrep = 4)
    set.seed(12)
    singles <- lapply(seq_len(4), function(i) unequal_prob_wr(hits, method))
    expect_identical(
      batch$sample, do.call(cbind, lapply(singles, `[[`, "sample"))
    )
    expect_identical(batch$hits, do.call(cbind, lapply(singles, `[[`, "hits")))
  }
})

test_that("random-size WOR batching uses the selected raw hook", {
  local_mocked_bindings(
    .get_builtin_spec = function(method, context) {
      list(draw = function(pik, prn = NULL) which(pik > 0.4))
    },
    .method_is_fixed_size = function(method, context) FALSE,
    .package = "sondage"
  )

  s <- sondage:::.batch_wor(c(0.2, 0.6, 0.2), "synthetic_random", 2L)
  expect_identical(s$sample, list(2L, 2L))
})

test_that("shared spatial C validation remains defensive", {
  pik <- rep(0.5, 4)
  symbols <- list(sondage:::C_lpm2, sondage:::C_scps)
  for (symbol in symbols) {
    expect_error(.Call(symbol, as.double(pik), 1:4, 1e-10), "numeric matrix")
    expect_error(
      .Call(symbol, as.double(pik), matrix(as.double(1:6), nrow = 3), 1e-10),
      "does not match"
    )
    expect_error(
      .Call(symbol, as.double(pik), matrix(numeric(), nrow = 4), 1e-10),
      "at least one column"
    )
  }
})

test_that("cube option validation agrees between single and batch paths", {
  pik <- rep(0.5, 6)
  aux <- matrix(as.double(seq_along(pik)), ncol = 1)
  for (nrep in c(1L, 2L)) {
    expect_no_error(
      balanced_wor(
        pik, aux = aux, nrep = nrep, condition_aux = TRUE, qr_tol = 0
      )
    )
    expect_error(
      balanced_wor(pik, aux = aux, nrep = nrep, condition_aux = 1),
      "condition_aux"
    )
    expect_error(
      balanced_wor(pik, aux = aux, nrep = nrep, qr_tol = -1),
      "non-negative"
    )
  }
})

## Audit #9: tied order-sampling keys stay correct

test_that("SPS and Pareto select exactly n with fully tied keys", {
  N <- 5000L
  p <- rep(0.5, N)
  u <- rep(0.5, N)
  for (m in c("sps", "pareto")) {
    s <- unequal_prob_wor(p, m, prn = u)
    expect_identical(length(s$sample), 2500L, info = m)
    expect_false(anyDuplicated(s$sample) > 0, info = m)
  }
})
