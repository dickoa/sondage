test_that("sps returns correct object structure", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "sps")

  expect_s3_class(s, c("unequal_prob", "wor", "sondage_sample"))
  expect_equal(s$method, "sps")
  expect_true(s$fixed_size)
  expect_equal(s$N, 4L)
  expect_equal(s$n, 2L)
  expect_equal(s$pik, pik)
})

test_that("sps produces valid samples", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "sps")

  expect_equal(length(s$sample), 2L)
  expect_true(all(s$sample >= 1L & s$sample <= 4L))
  expect_equal(length(unique(s$sample)), length(s$sample))
  expect_equal(s$sample, sort(s$sample))
})

test_that("sps inclusion probabilities match nominal pik", {
  set.seed(42)
  pik <- c(0.15, 0.25, 0.35, 0.45, 0.80)
  nrep <- 5000L
  sim <- unequal_prob_wor(pik, method = "sps", nrep = nrep)
  emp <- rowMeans(apply(sim, 2, function(col) {
    tabulate(col, nbins = length(pik))
  }))
  expect_equal(emp, pik, tolerance = 0.05)
})

test_that("sps with prn is deterministic", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.3, 0.7, 0.1, 0.5)
  s1 <- unequal_prob_wor(pik, method = "sps", prn = prn)
  s2 <- unequal_prob_wor(pik, method = "sps", prn = prn)
  expect_identical(s1$sample, s2$sample)
})

test_that("sps top-up coordination: n subset of n+1", {
  set.seed(1)
  prn <- runif(10)
  pik_small <- inclusion_prob(1:10, n = 3)
  pik_large <- inclusion_prob(1:10, n = 4)

  s_small <- unequal_prob_wor(pik_small, method = "sps", prn = prn)
  s_large <- unequal_prob_wor(pik_large, method = "sps", prn = prn)

  expect_true(all(s_small$sample %in% s_large$sample))
})

test_that("sps handles certainty units", {
  pik <- c(0.999999, 0.999999, 0.3, 0.4, 0.3)
  s <- unequal_prob_wor(pik, method = "sps")
  expect_true(all(c(1L, 2L) %in% s$sample))
})

test_that("sps works with n = 1", {
  pik <- c(0.1, 0.2, 0.3, 0.4)
  s <- unequal_prob_wor(pik, method = "sps")
  expect_equal(length(s$sample), 1L)
})

test_that("pareto returns correct object structure", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  s <- unequal_prob_wor(pik, method = "pareto")

  expect_s3_class(s, c("unequal_prob", "wor", "sondage_sample"))
  expect_equal(s$method, "pareto")
  expect_true(s$fixed_size)
  expect_equal(s$N, 4L)
  expect_equal(s$n, 2L)
})

test_that("pareto produces valid samples", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "pareto")

  expect_equal(length(s$sample), 2L)
  expect_true(all(s$sample >= 1L & s$sample <= 4L))
  expect_equal(length(unique(s$sample)), length(s$sample))
  expect_equal(s$sample, sort(s$sample))
})

test_that("pareto inclusion probabilities match nominal pik", {
  set.seed(42)
  pik <- c(0.15, 0.25, 0.35, 0.45, 0.80)
  nrep <- 10000L
  sim <- unequal_prob_wor(pik, method = "pareto", nrep = nrep)
  emp <- rowMeans(apply(sim, 2, function(col) {
    tabulate(col, nbins = length(pik))
  }))
  expect_equal(emp, pik, tolerance = 0.03)
})

test_that("pareto with prn is deterministic", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.3, 0.7, 0.1, 0.5)
  s1 <- unequal_prob_wor(pik, method = "pareto", prn = prn)
  s2 <- unequal_prob_wor(pik, method = "pareto", prn = prn)
  expect_identical(s1$sample, s2$sample)
})

test_that("pareto top-up coordination: n subset of n+1", {
  set.seed(1)
  prn <- runif(10)
  pik_small <- inclusion_prob(1:10, n = 3)
  pik_large <- inclusion_prob(1:10, n = 4)

  s_small <- unequal_prob_wor(pik_small, method = "pareto", prn = prn)
  s_large <- unequal_prob_wor(pik_large, method = "pareto", prn = prn)

  expect_true(all(s_small$sample %in% s_large$sample))
})

test_that("pareto handles certainty units", {
  pik <- c(0.999999, 0.999999, 0.3, 0.4, 0.3)
  s <- unequal_prob_wor(pik, method = "pareto")
  expect_true(all(c(1L, 2L) %in% s$sample))
})

test_that("pareto works with n = 1", {
  pik <- c(0.1, 0.2, 0.3, 0.4)
  s <- unequal_prob_wor(pik, method = "pareto")
  expect_equal(length(s$sample), 1L)
})

test_that("sps batch mode returns correct dimensions", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  sim <- unequal_prob_wor(pik, method = "sps", nrep = 50)
  expect_equal(dim(sim), c(2L, 50L))
})

test_that("pareto batch mode returns correct dimensions", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  sim <- unequal_prob_wor(pik, method = "pareto", nrep = 50)
  expect_equal(dim(sim), c(2L, 50L))
})

test_that("batch with prn produces identical columns", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.3, 0.7, 0.1, 0.5)

  sim_sps <- unequal_prob_wor(pik, method = "sps", nrep = 5, prn = prn)
  expect_true(all(apply(sim_sps, 2, identical, sim_sps[, 1])))

  sim_par <- unequal_prob_wor(pik, method = "pareto", nrep = 5, prn = prn)
  expect_true(all(apply(sim_par, 2, identical, sim_par[, 1])))
})

test_that("joint_inclusion_prob works for sps and pareto", {
  pik <- c(0.2, 0.3, 0.5)
  for (m in c("sps", "pareto")) {
    s <- unequal_prob_wor(pik, method = m)
    pikl <- joint_inclusion_prob(s)
    expect_equal(dim(pikl), c(3L, 3L))
    expect_equal(diag(pikl), pik)
    expect_true(isSymmetric(pikl))
  }
})

test_that("poisson with prn is deterministic", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.3, 0.7, 0.1, 0.5)
  s1 <- unequal_prob_wor(pik, method = "poisson", prn = prn)
  s2 <- unequal_prob_wor(pik, method = "poisson", prn = prn)
  expect_identical(s1$sample, s2$sample)
})

test_that("poisson prn selects units where prn < pik", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.1, 0.5, 0.3, 0.9)
  s <- unequal_prob_wor(pik, method = "poisson", prn = prn)
  expected <- which(prn < pik)
  expect_equal(s$sample, expected)
})

test_that("poisson batch with prn produces identical lists", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.1, 0.5, 0.3, 0.9)
  sim <- unequal_prob_wor(pik, method = "poisson", nrep = 5, prn = prn)
  for (i in 2:5) {
    expect_identical(sim[[i]], sim[[1]])
  }
})

test_that("bernoulli with prn is deterministic", {
  prn <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.2, 0.4, 0.6, 0.8, 0.05)
  s1 <- equal_prob_wor(10, 3, method = "bernoulli", prn = prn)
  s2 <- equal_prob_wor(10, 3, method = "bernoulli", prn = prn)
  expect_identical(s1$sample, s2$sample)
})

test_that("bernoulli prn selects units where prn < p", {
  prn <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.2, 0.4, 0.6, 0.8, 0.05)
  s <- equal_prob_wor(10, 3, method = "bernoulli", prn = prn)
  p <- 3 / 10
  expected <- which(prn < p)
  expect_equal(s$sample, expected)
})

test_that("bernoulli batch with prn produces identical lists", {
  prn <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.2, 0.4, 0.6, 0.8, 0.05)
  sim <- equal_prob_wor(10, 3, method = "bernoulli", nrep = 5, prn = prn)
  for (i in 2:5) {
    expect_identical(sim[[i]], sim[[1]])
  }
})

test_that("warning when prn used with unsupported WOR methods", {
  pik <- c(0.2, 0.4, 0.6, 0.8)
  prn <- c(0.3, 0.7, 0.1, 0.5)
  expect_warning(
    unequal_prob_wor(pik, method = "cps", prn = prn),
    "prn is not used by method 'cps'"
  )
  expect_warning(
    unequal_prob_wor(pik, method = "brewer", prn = prn),
    "prn is not used by method 'brewer'"
  )
  expect_warning(
    unequal_prob_wor(pik, method = "systematic", prn = prn),
    "prn is not used by method 'systematic'"
  )
})

test_that("warning when prn used with WR methods", {
  hits <- c(0.4, 0.8, 0.5, 0.6, 0.7)
  prn <- runif(5)
  expect_warning(
    unequal_prob_wr(hits, method = "chromy", prn = prn),
    "prn is not used by method 'chromy'"
  )
  expect_warning(
    unequal_prob_wr(hits, method = "multinomial", prn = prn),
    "prn is not used by method 'multinomial'"
  )
})

test_that("warning when prn used with unsupported equal_prob_wor methods", {
  prn <- runif(10)
  expect_warning(
    equal_prob_wor(10, 3, method = "srs", prn = prn),
    "prn is not used by method 'srs'"
  )
  expect_warning(
    equal_prob_wor(10, 3, method = "systematic", prn = prn),
    "prn is not used by method 'systematic'"
  )
})

test_that("warning when prn used with equal_prob_wr", {
  prn <- runif(10)
  expect_warning(
    equal_prob_wr(10, 3, prn = prn),
    "prn is not used by method 'srs'"
  )
})

test_that("check_prn rejects invalid inputs", {
  expect_error(
    unequal_prob_wor(c(0.5, 0.5), method = "sps", prn = "a"),
    "prn must be a numeric"
  )
  expect_error(
    unequal_prob_wor(c(0.5, 0.5), method = "sps", prn = c(0.3)),
    "prn must have length 2"
  )
  expect_error(
    unequal_prob_wor(c(0.5, 0.5), method = "sps", prn = c(0.3, NA)),
    "missing values in prn"
  )
  expect_error(
    unequal_prob_wor(c(0.5, 0.5), method = "sps", prn = c(0.0, 0.5)),
    "open interval"
  )
  expect_error(
    unequal_prob_wor(c(0.5, 0.5), method = "sps", prn = c(1.0, 0.5)),
    "open interval"
  )
})
