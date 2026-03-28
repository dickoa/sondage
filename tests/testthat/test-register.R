# --- Registration API ---

test_that("register_method validates inputs", {
  expect_error(register_method(123, "wor", sample_fn = identity), "non-empty character")
  expect_error(register_method("", "wor", sample_fn = identity), "non-empty character")
  expect_error(register_method("x", "wor", sample_fn = "not_fn"), "must be a function")
  expect_error(register_method("x", "wor", sample_fn = identity, joint_fn = 42),
               "must be a function or NULL")
  expect_error(register_method("x", "wor", sample_fn = identity, fixed_size = NA),
               "TRUE or FALSE")
  expect_error(register_method("x", "wor", sample_fn = identity, supports_prn = "yes"),
               "TRUE or FALSE")
  expect_error(register_method("x", "badtype", sample_fn = identity),
               "should be one of")
})

test_that("register_method rejects built-in names", {
  expect_error(register_method("cps", "wor", sample_fn = identity), "built-in")
  expect_error(register_method("brewer", "wor", sample_fn = identity), "built-in")
  expect_error(register_method("chromy", "wr", sample_fn = identity), "built-in")
  expect_error(register_method("srs", "wor", sample_fn = identity), "built-in")
  expect_error(register_method("systematic", "wor", sample_fn = identity), "built-in")
})

test_that("register/unregister lifecycle works", {
  on.exit(unregister_method("test_method"), add = TRUE)

  expect_false(is_registered_method("test_method"))
  expect_equal(length(registered_methods()), 0L)

  register_method("test_method", "wor", sample_fn = identity)
  expect_true(is_registered_method("test_method"))
  expect_equal(registered_methods(), "test_method")

  expect_true(unregister_method("test_method"))
  expect_false(is_registered_method("test_method"))
  expect_false(unregister_method("test_method"))
})

# --- Dispatch with a toy WOR method ---

toy_wor_sample <- function(pik, n = NULL, prn = NULL, ...) {
  # Deterministic: always pick the n units with largest pik
  order(pik, decreasing = TRUE)[seq_len(n)]
}

test_that("registered WOR method dispatches through unequal_prob_wor", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")

  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_equal(s$method, "toy_wor")
  expect_equal(s$N, 6L)
  expect_equal(s$n, 2L)
  expect_true(s$fixed_size)
  expect_equal(sort(s$sample), c(5L, 6L))
  expect_equal(s$pik, pik)
})

test_that("registered WOR method works in batch mode", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor", nrep = 5)

  expect_equal(dim(s$sample), c(2L, 5L))
  # All replicates identical (deterministic sampler)
  for (i in seq_len(5)) {
    expect_equal(sort(s$sample[, i]), c(5L, 6L))
  }
})

test_that("registered method without joint_fn errors on joint_inclusion_prob", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")
  expect_error(joint_inclusion_prob(s), "not implemented for method 'toy_wor'")
})

test_that("registered method with joint_fn works", {
  on.exit(unregister_method("toy_wor"), add = TRUE)

  toy_joint <- function(pik, sample_idx = NULL, ...) {
    # Fake: just outer product
    if (!is.null(sample_idx)) pik <- pik[sample_idx]
    J <- outer(pik, pik)
    diag(J) <- pik
    J
  }

  register_method("toy_wor", "wor",
                   sample_fn = toy_wor_sample, joint_fn = toy_joint)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")

  pikl <- joint_inclusion_prob(s)
  expect_equal(nrow(pikl), 6L)
  expect_equal(ncol(pikl), 6L)
  expect_equal(diag(pikl), pik)

  pikl_s <- joint_inclusion_prob(s, sampled_only = TRUE)
  expect_equal(nrow(pikl_s), 2L)
})

test_that("sampling_cov works for registered WOR method with joint_fn", {
  on.exit(unregister_method("toy_wor"), add = TRUE)

  toy_joint <- function(pik, sample_idx = NULL, ...) {
    if (!is.null(sample_idx)) pik <- pik[sample_idx]
    J <- outer(pik, pik)
    diag(J) <- pik
    J
  }

  register_method("toy_wor", "wor",
                   sample_fn = toy_wor_sample, joint_fn = toy_joint)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")

  delta <- sampling_cov(s)
  expect_equal(nrow(delta), 6L)
})

test_that("PRN warning for registered method without PRN support", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample,
                   supports_prn = FALSE)

  pik <- c(0.5, 0.5)
  expect_warning(
    unequal_prob_wor(pik, method = "toy_wor", prn = c(0.3, 0.7)),
    "prn is not used"
  )
})

test_that("print works for registered method", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")
  out <- capture.output(print(s))
  expect_true(any(grepl("toy_wor", out)))
})

# --- Dispatch with a toy WR method ---

toy_wr_sample <- function(hits, n = NULL, prn = NULL, ...) {
  sample.int(length(hits), size = n, replace = TRUE, prob = hits)
}

test_that("registered WR method dispatches through unequal_prob_wr", {
  on.exit(unregister_method("toy_wr"), add = TRUE)
  register_method("toy_wr", "wr", sample_fn = toy_wr_sample)

  hits <- c(0.5, 1.0, 0.5)
  s <- unequal_prob_wr(hits, method = "toy_wr")

  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wr")
  expect_equal(s$method, "toy_wr")
  expect_equal(s$N, 3L)
  expect_equal(s$n, 2L)
})

test_that("registered WR method works in batch mode", {
  on.exit(unregister_method("toy_wr"), add = TRUE)
  register_method("toy_wr", "wr", sample_fn = toy_wr_sample)

  hits <- c(0.5, 1.0, 0.5)
  s <- unequal_prob_wr(hits, method = "toy_wr", nrep = 10)

  expect_equal(dim(s$sample), c(2L, 10L))
  expect_equal(dim(s$hits), c(3L, 10L))
})

# --- Integration with sampling::UPtille ---

test_that("UPtille from sampling package works end-to-end", {
  skip_if_not_installed("sampling")
  on.exit(unregister_method("tille"), add = TRUE)

  tille_sample <- function(pik, n = NULL, prn = NULL, ...) {
    which(sampling::UPtille(pik) == 1L)
  }

  tille_joint <- function(pik, sample_idx = NULL, ...) {
    pikl <- sampling::UPtillepi2(pik)
    if (!is.null(sample_idx)) {
      pikl <- pikl[sample_idx, sample_idx, drop = FALSE]
    }
    pikl
  }

  register_method("tille", "wor",
                   sample_fn = tille_sample,
                   joint_fn = tille_joint,
                   fixed_size = TRUE)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  set.seed(42)
  s <- unequal_prob_wor(pik, method = "tille")

  # Basic structure

  expect_s3_class(s, c("unequal_prob", "wor", "sondage_sample"))
  expect_equal(s$method, "tille")
  expect_equal(length(s$sample), 2L)
  expect_true(all(s$sample >= 1L & s$sample <= 6L))

  # Joint inclusion probs
  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(6L, 6L))
  expect_equal(diag(pikl), pik, tolerance = 1e-10)
  expect_true(isSymmetric(pikl))

  # sampled_only submatrix
  pikl_s <- joint_inclusion_prob(s, sampled_only = TRUE)
  expect_equal(dim(pikl_s), c(2L, 2L))
  expect_equal(pikl_s, pikl[s$sample, s$sample])

  # Sampling covariance
  delta <- sampling_cov(s)
  expect_equal(dim(delta), c(6L, 6L))

  # inclusion_prob generic
  expect_equal(inclusion_prob(s), pik)
})

test_that("UPtille batch mode produces valid replicates", {
  skip_if_not_installed("sampling")
  on.exit(unregister_method("tille"), add = TRUE)

  tille_sample <- function(pik, n = NULL, prn = NULL, ...) {
    which(sampling::UPtille(pik) == 1L)
  }

  register_method("tille", "wor", sample_fn = tille_sample, fixed_size = TRUE)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "tille", nrep = 100)

  expect_equal(dim(s$sample), c(2L, 100L))

  # Empirical inclusion probs should be roughly close to pik
  freq <- tabulate(s$sample, nbins = 6) / 100
  expect_equal(freq, pik, tolerance = 0.2)
})

test_that("UPtille sampling_cov(weighted=TRUE) works", {
  skip_if_not_installed("sampling")
  on.exit(unregister_method("tille"), add = TRUE)

  tille_sample <- function(pik, n = NULL, prn = NULL, ...) {
    which(sampling::UPtille(pik) == 1L)
  }
  tille_joint <- function(pik, sample_idx = NULL, ...) {
    pikl <- sampling::UPtillepi2(pik)
    if (!is.null(sample_idx)) pikl <- pikl[sample_idx, sample_idx, drop = FALSE]
    pikl
  }

  register_method("tille", "wor",
                   sample_fn = tille_sample, joint_fn = tille_joint)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  set.seed(42)
  s <- unequal_prob_wor(pik, method = "tille")

  syg <- sampling_cov(s, weighted = TRUE)
  expect_equal(dim(syg), c(6L, 6L))
  # Off-diagonal should be non-positive for well-behaved WOR designs
  off_diag <- syg[lower.tri(syg)]
  off_diag <- off_diag[!is.na(off_diag)]
  expect_true(all(off_diag <= 1e-10))
})

# --- Custom Poisson with PRN (random-size + sample coordination) ---

poisson_prn_sample <- function(pik, n = NULL, prn = NULL, ...) {
  N <- length(pik)
  if (is.null(prn)) prn <- runif(N)
  which(prn < pik)
}

poisson_prn_joint <- function(pik, sample_idx = NULL, ...) {
  # Independent draws: pi_ij = pi_i * pi_j
  if (!is.null(sample_idx)) pik <- pik[sample_idx]
  J <- outer(pik, pik)
  diag(J) <- pik
  J
}

test_that("registered Poisson with PRN produces correct sample", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   joint_fn = poisson_prn_joint,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.2, 0.8, 0.5, 0.9, 0.1)
  prn <- c(0.3, 0.1, 0.4, 0.2, 0.9)

  s <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn)

  # prn < pik selects units 2, 4 (0.1<0.8, 0.2<0.9)
  # unit 3: prn=0.4 < pik=0.5 -> selected too
  expect_equal(s$sample, c(2L, 3L, 4L))
  expect_false(s$fixed_size)
  expect_equal(s$n, sum(pik))
  expect_equal(s$method, "poisson_prn")
})

test_that("PRN coordination: same PRN gives same sample", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   joint_fn = poisson_prn_joint,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.2, 0.8, 0.5, 0.9, 0.1)
  prn <- c(0.15, 0.5, 0.3, 0.85, 0.05)

  s1 <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn)
  s2 <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn)
  expect_equal(s1$sample, s2$sample)
})

test_that("PRN coordination: different PRN can change sample", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   joint_fn = poisson_prn_joint,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.5, 0.5, 0.5, 0.5)
  prn_a <- c(0.1, 0.6, 0.1, 0.6)  # selects units 1, 3
  prn_b <- c(0.6, 0.1, 0.6, 0.1)  # selects units 2, 4

  s_a <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn_a)
  s_b <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn_b)
  expect_equal(s_a$sample, c(1L, 3L))
  expect_equal(s_b$sample, c(2L, 4L))
})

test_that("PRN + nrep > 1 errors for registered method", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.3, 0.7)
  expect_error(
    unequal_prob_wor(pik, method = "poisson_prn", prn = c(0.5, 0.5), nrep = 2),
    "prn and nrep > 1 cannot be used together"
  )
})

test_that("random-size batch mode returns list of varying-length samples", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.3, 0.5, 0.7, 0.2, 0.8)
  set.seed(123)
  s <- unequal_prob_wor(pik, method = "poisson_prn", nrep = 50)

  expect_true(is.list(s$sample))
  expect_equal(length(s$sample), 50L)
  expect_false(s$fixed_size)

  # Sample sizes should vary (Poisson is random-size)
  sizes <- lengths(s$sample)
  expect_true(length(unique(sizes)) > 1L)

  # All indices in valid range
  all_idx <- unlist(s$sample)
  expect_true(all(all_idx >= 1L & all_idx <= 5L))
})

test_that("joint probs work for registered random-size method", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   joint_fn = poisson_prn_joint,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.3, 0.5, 0.7, 0.2, 0.8)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "poisson_prn")

  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(5L, 5L))
  expect_equal(diag(pikl), pik)
  # Independent: pi_ij = pi_i * pi_j
  expect_equal(pikl[1, 2], pik[1] * pik[2])

  # sampled_only
  pikl_s <- joint_inclusion_prob(s, sampled_only = TRUE)
  ns <- length(s$sample)
  expect_equal(dim(pikl_s), c(ns, ns))
})

test_that("sampling_cov works for registered random-size method", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method("poisson_prn", "wor",
                   sample_fn = poisson_prn_sample,
                   joint_fn = poisson_prn_joint,
                   fixed_size = FALSE,
                   supports_prn = TRUE)

  pik <- c(0.3, 0.5, 0.7, 0.2, 0.8)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "poisson_prn")

  delta <- sampling_cov(s)
  expect_equal(dim(delta), c(5L, 5L))
  # Independent design: Cov(I_i, I_j) = 0 for i != j
  expect_equal(delta[1, 2], 0)
  expect_equal(delta[3, 5], 0)
  # Diagonal: Var(I_i) = pi_i(1 - pi_i)
  expect_equal(diag(delta), pik * (1 - pik), tolerance = 1e-10)
})

# --- Unknown method still errors ---

test_that("unknown method (not built-in, not registered) still errors", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "nonexistent"))
  expect_error(unequal_prob_wr(c(1, 1), method = "nonexistent"))
})
