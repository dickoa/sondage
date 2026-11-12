## Without-replacement systematic draws stay the same

# Recorded from C_up_systematic before it gained the rescaling guard and the
# end pin it shares with "systematic_mr". The guard and the pin only act
# within rounding error of a boundary, so ordinary inputs must keep drawing
# the same units and consume the same random numbers.

test_that("systematic WOR samples match the values recorded before the guard", {
  cases <- list(
    plain = inclusion_prob(c(5, 12, 30, 60, 8, 40, 45, 12, 70, 25), 3),
    certain = c(1, 0.4, 0, 0.35, 1, 0.25),
    drift = {
      p <- inclusion_prob(c(3, 7, 11, 13, 17, 19, 23), 4)
      p[1] <- p[1] + 3e-11
      p
    },
    tiny = c(1e-12, 0.5, 0.5 - 1e-12, 1)
  )
  expected <- list(
    plain = list(c(3L, 6L, 9L), c(4L, 7L, 10L), c(4L, 7L, 9L)),
    certain = list(c(1L, 2L, 5L), c(1L, 5L, 6L), c(1L, 4L, 5L)),
    drift = list(c(2L, 4L, 6L, 7L), 4:7, c(3L, 5L, 6L, 7L)),
    tiny = list(c(2L, 4L), 3:4, 3:4)
  )
  seeds <- c(1L, 42L, 2026L)
  for (nm in names(cases)) {
    for (i in seq_along(seeds)) {
      set.seed(seeds[[i]])
      s <- unequal_prob_wor(cases[[nm]], method = "systematic")$sample
      expect_identical(s, expected[[nm]][[i]], info = paste(nm, seeds[[i]]))
    }
  }

  set.seed(9)
  s <- unequal_prob_wor(cases$plain, method = "systematic", nrep = 3)$sample
  expect_identical(s, matrix(c(3L, 6L, 9L, 1L, 4L, 8L, 3L, 6L, 9L), 3L, 3L))
})

test_that("a systematic WOR draw uses one uniform, and none when all are certain", {
  pik <- inclusion_prob(c(5, 12, 30, 60, 8, 40, 45, 12, 70, 25), 3)
  set.seed(3)
  unequal_prob_wor(pik, method = "systematic")
  after_draw <- runif(1)
  set.seed(3)
  runif(1)
  expect_identical(after_draw, runif(1))

  set.seed(3)
  unequal_prob_wor(c(1, 0, 1), method = "systematic")
  after_draw <- runif(1)
  set.seed(3)
  expect_identical(after_draw, runif(1))
})

## Systematic PPS with minimum replacement

# Expected hits with whole parts 0, 1 and 2 and no fractional part near 0
# or 1, so that every hit count can occur.
mr_hits <- expected_hits(c(5, 120, 30, 60, 8, 300, 45, 12, 70, 25, 450), 7)

# Fractional part of each expected hit, before the walk rescales it.
mr_frac <- function(e) e - floor(e)

test_that("systematic_mr returns a WR design object", {
  set.seed(1)
  s <- unequal_prob_wr(mr_hits, method = "systematic_mr")
  expect_s3_class(s, c("unequal_prob", "wr", "sondage_sample"), exact = TRUE)
  expect_identical(s$method, "systematic_mr")
  expect_identical(s$n, 7L)
  expect_true(s$fixed_size)
  expect_identical(s$sample, rep(seq_along(s$hits), s$hits))
  expect_equal(expected_hits(s), mr_hits)
  expect_output(print(s), "Unequal prob WR \\[systematic_mr\\]")
})

test_that("each unit gets floor(e) or ceiling(e) hits and the hits sum to n", {
  e <- c(0, 1, 2.5, 0.3, 0.7, 0.25, 0.75, 1.5)
  set.seed(4)
  s <- unequal_prob_wr(e, method = "systematic_mr", nrep = 2000)
  expect_true(all(s$hits >= floor(e) & s$hits <= ceiling(e)))
  expect_true(all(colSums(s$hits) == 7L))
  expect_true(all(s$hits[1, ] == 0L))
  expect_true(all(s$hits[2, ] == 1L))

  set.seed(4)
  one <- unequal_prob_wr(e, method = "systematic_mr")
  expect_identical(one$hits, s$hits[, 1])
  expect_identical(one$sample, s$sample[, 1])
})

test_that("sums a few ulps off a whole number keep the hits in range", {
  e <- mr_hits
  for (drift in c(-1e-11, 1e-11, 64 * .Machine$double.eps)) {
    d <- e
    d[3] <- d[3] + drift
    set.seed(5)
    s <- unequal_prob_wr(d, method = "systematic_mr", nrep = 500)
    expect_true(all(s$hits >= floor(d) & s$hits <= ceiling(d)))
    expect_true(all(colSums(s$hits) == 7L))
  }
})

test_that("whole expected hits give a deterministic draw that uses no uniform", {
  e <- c(2, 0, 1, 3, 1)
  set.seed(6)
  s <- unequal_prob_wr(e, method = "systematic_mr")
  after_draw <- runif(1)
  set.seed(6)
  expect_identical(after_draw, runif(1))
  expect_identical(s$hits, as.integer(e))
  expect_identical(s$sample, c(1L, 1L, 3L, 4L, 4L, 4L, 5L))
  expect_identical(
    joint_expected_hits(s),
    outer(e, e)
  )
})

test_that("mean hits and cross-products match the exact moments", {
  skip_on_cran()
  B <- 40000L
  set.seed(7)
  H <- unequal_prob_wr(mr_hits, method = "systematic_mr", nrep = B)$hits
  J <- joint_expected_hits(unequal_prob_wr(mr_hits, method = "systematic_mr"))

  f <- mr_frac(mr_hits)
  se_mean <- sqrt(f * (1 - f) / B)
  expect_true(all(abs(rowMeans(H) - mr_hits) <= 5 * se_mean + 1e-12))

  cross <- tcrossprod(H) / B
  sd_cross <- sqrt(pmax((tcrossprod(H^2) / B) - cross^2, 0))
  expect_true(all(abs(cross - J) <= 5 * sd_cross / sqrt(B) + 1e-12))
})

test_that("the joint expected hits have the exact structure of the design", {
  set.seed(8)
  s <- unequal_prob_wr(mr_hits, method = "systematic_mr")
  J <- joint_expected_hits(s)
  f <- mr_frac(mr_hits)

  expect_true(isSymmetric(J))
  expect_equal(rowSums(J), 7 * mr_hits, tolerance = 1e-12)
  expect_equal(diag(J), mr_hits^2 + f * (1 - f), tolerance = 1e-12)
  expect_true(all(J >= 0))

  # A unit with a whole part is hit in every sample.
  big <- which(mr_hits >= 1)
  expect_true(all(J[big, -big] >= rep(mr_hits[-big], each = length(big)) -
    1e-12))

  idx <- which(s$hits > 0)
  expect_identical(
    unname(joint_expected_hits(s, sampled_only = TRUE)),
    unname(J[idx, idx])
  )
  expect_identical(joint_expected_hits(s, nsim = 5L), J)
})

test_that("pairs that never occur together are exactly zero", {
  set.seed(9)
  s <- unequal_prob_wr(mr_hits, method = "systematic_mr")
  J <- joint_expected_hits(s)

  # Brute-force the support of the design over a fine grid of starts.
  e <- mr_hits
  f <- (e - floor(e)) * (7 - sum(floor(e))) / sum(e - floor(e))
  hi <- cumsum(f)
  lo <- hi - f
  u <- (seq_len(20000) - 0.5) / 20000
  I <- vapply(u, function(v) {
    as.integer(floor(hi - v) > floor(lo - v))
  }, integer(length(e)))
  ever <- tcrossprod(I + floor(e)) > 0
  expect_true(all(J[!ever] == 0))
  expect_true(all(J[ever] > 0))
  expect_true(any(!ever))

  expect_warning(
    m <- sampling_cov(s, weighted = TRUE),
    "pairwise expected hits are zero"
  )
  expect_true(all(is.na(m[!ever])))
  expect_false(anyNA(m[ever & row(m) != col(m)]))
})

test_that("with every expected hit at least 1, all pairs are positive", {
  e <- c(1.4, 2.6, 1, 3)
  set.seed(10)
  s <- unequal_prob_wr(e, method = "systematic_mr")
  expect_true(all(joint_expected_hits(s) > 0))
  expect_no_warning(sampling_cov(s, weighted = TRUE))
})

test_that("a whole-number hit stays whole through n * prob", {
  # sum(h) is 2 + 4e-16, so n * prob returns the first unit as
  # 0.99999999999999978. The draw and the joint expected hits both treat
  # it as one certain hit.
  h <- c(1, 0.1, 0.2, 0.7 + 4e-16)
  set.seed(11)
  s <- unequal_prob_wr(h, method = "systematic_mr", nrep = 2000)
  expect_true(all(s$hits[1, ] == 1L))
  s1 <- unequal_prob_wr(h, method = "systematic_mr")
  expect_lt(expected_hits(s1)[1], 1)
  J <- joint_expected_hits(s1)
  expect_identical(J[1, 1], 1)
  expect_equal(J[1, -1], h[-1], tolerance = 1e-12)
})

test_that("with every expected hit below 1 it is the WOR systematic design", {
  pik <- inclusion_prob(c(5, 12, 30, 60, 8, 40, 45, 12, 70, 25), 3)
  for (seed in 1:200) {
    set.seed(seed)
    wor <- unequal_prob_wor(pik, method = "systematic")$sample
    set.seed(seed)
    wr <- unequal_prob_wr(pik, method = "systematic_mr")$sample
    expect_identical(wr, wor, info = seed)
  }
  s_wor <- unequal_prob_wor(pik, method = "systematic")
  s_wr <- unequal_prob_wr(pik, method = "systematic_mr")
  expect_equal(
    joint_expected_hits(s_wr),
    joint_inclusion_prob(s_wor),
    tolerance = 1e-12
  )
})

test_that("a batch of draws matches the same number of single draws", {
  set.seed(12)
  batch <- unequal_prob_wr(mr_hits, method = "systematic_mr", nrep = 4)
  set.seed(12)
  singles <- lapply(1:4, function(i) {
    unequal_prob_wr(mr_hits, method = "systematic_mr")
  })
  expect_identical(batch$sample, do.call(cbind, lapply(singles, `[[`, "sample")))
  expect_identical(batch$hits, do.call(cbind, lapply(singles, `[[`, "hits")))
  expect_identical(dim(batch$sample), c(7L, 4L))
  expect_identical(dim(batch$hits), c(length(mr_hits), 4L))
})

test_that("systematic_mr refuses what the WR methods refuse", {
  expect_error(
    unequal_prob_wr(c(0.5, 0.5 + 5e-5), method = "systematic_mr"),
    "not close to an integer"
  )
  expect_error(
    unequal_prob_wr(mr_hits, method = "systematic_mr", prn = runif(11)),
    "prn"
  )
  expect_error(
    unequal_prob_wr(mr_hits, method = "systematic_mr", n = 7),
    "sum\\(hits\\)"
  )
  expect_error(
    unequal_prob_wr(mr_hits, method = "systematic_mr", foo = 1),
    "unused argument"
  )
  expect_error(
    unequal_prob_wr(mr_hits, method = "systematic_mr", nrep = 2, foo = 1),
    "unused argument"
  )
})

test_that("method_spec() describes systematic_mr", {
  spec <- method_spec("systematic_mr")
  expect_identical(spec$dispatcher, "unequal_prob_wr")
  expect_identical(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_identical(spec$variance_family, "wr")
  expect_false(spec$supports_prn)
  expect_identical(spec$probabilities, "exact")
  expect_identical(
    method_spec("systematic_mr", dispatcher = "unequal_prob_wr"),
    spec
  )
  expect_error(
    register_method(
      "systematic_mr",
      "wr",
      sample_fn = function(hits, n, prn, ...) seq_len(n)
    ),
    "built-in"
  )
})
