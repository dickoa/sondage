test_that("register_method validates inputs", {
  expect_error(
    register_method(123, "wor", sample_fn = identity),
    "non-empty character"
  )
  expect_error(
    register_method("", "wor", sample_fn = identity),
    "non-empty character"
  )
  expect_error(
    register_method("x", "wor", sample_fn = "not_fn"),
    "must be a function"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, joint_fn = 42),
    "must be a function or NULL"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, fixed_size = NA),
    "TRUE or FALSE"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_prn = "yes"),
    "TRUE or FALSE"
  )
  expect_error(
    register_method("x", "badtype", sample_fn = identity),
    "should be one of"
  )
})

test_that("register_method rejects built-in names", {
  expect_error(register_method("cps", "wor", sample_fn = identity), "built-in")
  expect_error(
    register_method("brewer", "wor", sample_fn = identity),
    "built-in"
  )
  expect_error(
    register_method("chromy", "wr", sample_fn = identity),
    "built-in"
  )
  expect_error(register_method("srs", "wor", sample_fn = identity), "built-in")
  expect_error(
    register_method("systematic", "wor", sample_fn = identity),
    "built-in"
  )
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

# Dispatch with a toy WOR method
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
    if (!is.null(sample_idx)) {
      pik <- pik[sample_idx]
    }
    J <- outer(pik, pik)
    diag(J) <- pik
    J
  }

  register_method(
    "toy_wor",
    "wor",
    sample_fn = toy_wor_sample,
    joint_fn = toy_joint
  )

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
    if (!is.null(sample_idx)) {
      pik <- pik[sample_idx]
    }
    J <- outer(pik, pik)
    diag(J) <- pik
    J
  }

  register_method(
    "toy_wor",
    "wor",
    sample_fn = toy_wor_sample,
    joint_fn = toy_joint
  )

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- unequal_prob_wor(pik, method = "toy_wor")

  delta <- sampling_cov(s)
  expect_equal(nrow(delta), 6L)
})

test_that("PRN warning for registered method without PRN support", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method(
    "toy_wor",
    "wor",
    sample_fn = toy_wor_sample,
    supports_prn = FALSE
  )

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

# Dispatch with a toy WR method
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

# Integration with sampling::UPtille
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

  register_method(
    "tille",
    "wor",
    sample_fn = tille_sample,
    joint_fn = tille_joint,
    fixed_size = TRUE
  )

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  set.seed(123)
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
    if (!is.null(sample_idx)) {
      pikl <- pikl[sample_idx, sample_idx, drop = FALSE]
    }
    pikl
  }

  register_method(
    "tille",
    "wor",
    sample_fn = tille_sample,
    joint_fn = tille_joint
  )

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

# Custom Poisson with PRN (random-size + sample coordination)
poisson_prn_sample <- function(pik, n = NULL, prn = NULL, ...) {
  N <- length(pik)
  if (is.null(prn)) {
    prn <- runif(N)
  }
  which(prn < pik)
}

poisson_prn_joint <- function(pik, sample_idx = NULL, ...) {
  # Independent draws: pi_ij = pi_i * pi_j
  if (!is.null(sample_idx)) {
    pik <- pik[sample_idx]
  }
  J <- outer(pik, pik)
  diag(J) <- pik
  J
}

test_that("registered Poisson with PRN produces correct sample", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    joint_fn = poisson_prn_joint,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

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
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    joint_fn = poisson_prn_joint,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

  pik <- c(0.2, 0.8, 0.5, 0.9, 0.1)
  prn <- c(0.15, 0.5, 0.3, 0.85, 0.05)

  s1 <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn)
  s2 <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn)
  expect_equal(s1$sample, s2$sample)
})

test_that("PRN coordination: different PRN can change sample", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    joint_fn = poisson_prn_joint,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

  pik <- c(0.5, 0.5, 0.5, 0.5)
  prn_a <- c(0.1, 0.6, 0.1, 0.6) # selects units 1, 3
  prn_b <- c(0.6, 0.1, 0.6, 0.1) # selects units 2, 4

  s_a <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn_a)
  s_b <- unequal_prob_wor(pik, method = "poisson_prn", prn = prn_b)
  expect_equal(s_a$sample, c(1L, 3L))
  expect_equal(s_b$sample, c(2L, 4L))
})

test_that("PRN + nrep > 1 errors for registered method", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

  pik <- c(0.3, 0.7)
  expect_error(
    unequal_prob_wor(pik, method = "poisson_prn", prn = c(0.5, 0.5), nrep = 2),
    "prn and nrep > 1 cannot be used together"
  )
})

test_that("random-size batch mode returns list of varying-length samples", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

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
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    joint_fn = poisson_prn_joint,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

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
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    joint_fn = poisson_prn_joint,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

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

# method_spec() integration with registration
test_that("method_spec tracks register/unregister lifecycle", {
  on.exit(unregister_method("test_spec"), add = TRUE)

  expect_null(method_spec("test_spec"))

  register_method(
    "test_spec",
    "wor",
    sample_fn = identity,
    fixed_size = FALSE,
    supports_prn = TRUE
  )
  spec <- method_spec("test_spec")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)

  unregister_method("test_spec")
  expect_null(method_spec("test_spec"))
})

test_that("method_spec returns correct metadata for toy WOR sampler", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  spec <- method_spec("toy_wor")
  expect_equal(spec$type, "wor")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("method_spec returns correct metadata for toy WR sampler", {
  on.exit(unregister_method("toy_wr"), add = TRUE)
  register_method("toy_wr", "wr", sample_fn = toy_wr_sample)

  spec <- method_spec("toy_wr")
  expect_equal(spec$type, "wr")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
})

test_that("method_spec returns correct metadata for random-size PRN method", {
  on.exit(unregister_method("poisson_prn"), add = TRUE)
  register_method(
    "poisson_prn",
    "wor",
    sample_fn = poisson_prn_sample,
    fixed_size = FALSE,
    supports_prn = TRUE
  )

  spec <- method_spec("poisson_prn")
  expect_equal(spec$type, "wor")
  expect_false(spec$fixed_size)
  expect_true(spec$supports_prn)
})

# Registration validation for balanced methods
test_that("register_method validates balanced-specific flags", {
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_strata = NA),
    "TRUE or FALSE"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_strata = TRUE),
    "only applies to type"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_spread = NA),
    "TRUE or FALSE"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_spread = TRUE),
    "only applies to type"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_aux = NA),
    "TRUE or FALSE"
  )
  expect_error(
    register_method("x", "wor", sample_fn = identity, supports_aux = FALSE),
    "only applies to type"
  )
  expect_error(
    register_method("x", "balanced", sample_fn = identity, supports_prn = TRUE),
    "cannot use prn"
  )
  expect_error(
    register_method("cube", "balanced", sample_fn = identity),
    "built-in"
  )
})

# Dispatch with a toy balanced method
toy_balanced_sample <- function(pik, n = NULL, aux = NULL, ...) {
  # Deterministic: always pick the n units with largest pik
  order(pik, decreasing = TRUE)[seq_len(n)]
}

test_that("registered balanced method dispatches through balanced_wor", {
  on.exit(unregister_method("toy_bal"), add = TRUE)
  register_method("toy_bal", "balanced", sample_fn = toy_balanced_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- balanced_wor(pik, method = "toy_bal")

  expect_s3_class(s, "sondage_sample")
  expect_s3_class(s, "wor")
  expect_s3_class(s, "unequal_prob")
  expect_s3_class(s, "balanced")
  expect_equal(s$method, "toy_bal")
  expect_equal(s$N, 6L)
  expect_equal(s$n, 2L)
  expect_true(s$fixed_size)
  expect_equal(sort(s$sample), c(5L, 6L))
  expect_equal(s$pik, pik)
})

test_that("registered balanced method receives validated aux", {
  on.exit(unregister_method("toy_bal"), add = TRUE)

  seen_aux <- NULL
  register_method(
    "toy_bal",
    "balanced",
    sample_fn = function(pik, n = NULL, aux = NULL, ...) {
      seen_aux <<- aux
      seq_len(n)
    }
  )

  pik <- c(0.5, 0.5, 0.5, 0.5)
  balanced_wor(pik, aux = c(1, 2, 3, 4), method = "toy_bal")

  expect_true(is.matrix(seen_aux))
  expect_equal(dim(seen_aux), c(4L, 1L))
  expect_identical(storage.mode(seen_aux), "double")

  # aux validation fires before sample_fn is called
  expect_error(
    balanced_wor(pik, aux = matrix(1, 3, 1), method = "toy_bal"),
    "does not match"
  )
  expect_error(
    balanced_wor(pik, aux = matrix(NA_real_, 4, 1), method = "toy_bal"),
    "NA, NaN, or Inf"
  )
})

test_that("registered balanced method works in batch mode", {
  on.exit(unregister_method("toy_bal"), add = TRUE)
  register_method("toy_bal", "balanced", sample_fn = toy_balanced_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- balanced_wor(pik, method = "toy_bal", nrep = 5)

  expect_equal(dim(s$sample), c(2L, 5L))
  for (i in seq_len(5)) {
    expect_equal(sort(s$sample[, i]), c(5L, 6L))
  }
})

test_that("aux-only balanced method rejects strata and spread", {
  on.exit(unregister_method("toy_bal"), add = TRUE)
  register_method("toy_bal", "balanced", sample_fn = toy_balanced_sample)

  pik <- rep(0.5, 8)
  expect_error(
    balanced_wor(pik, strata = rep(1:2, each = 4), method = "toy_bal"),
    "does not support stratified"
  )
  expect_error(
    balanced_wor(pik, spread = matrix(runif(16), 8, 2), method = "toy_bal"),
    "does not support spatial spreading"
  )
})

test_that("built-in cube rejects spread", {
  pik <- rep(0.5, 8)
  expect_error(
    balanced_wor(pik, spread = matrix(runif(16), 8, 2)),
    "does not support spatial spreading"
  )
})

test_that("spread-only method with supports_aux = FALSE rejects aux", {
  on.exit(unregister_method("toy_spat"), add = TRUE)

  seen_aux <- "unset"
  register_method(
    "toy_spat",
    "balanced",
    sample_fn = function(pik, n = NULL, aux = NULL, spread = NULL, ...) {
      seen_aux <<- aux
      seq_len(n)
    },
    supports_aux = FALSE,
    supports_spread = TRUE
  )

  pik <- rep(0.5, 8)
  coords <- matrix(as.double(1:16), 8, 2)
  expect_error(
    balanced_wor(
      pik,
      aux = matrix(1, 8, 1),
      spread = coords,
      method = "toy_spat"
    ),
    "does not use auxiliary balancing variables"
  )

  # aux = NULL still flows through normally
  s <- balanced_wor(pik, spread = coords, method = "toy_spat")
  expect_null(seen_aux)
  expect_equal(s$n, 4L)
})

test_that("spread-capable balanced method receives validated spread", {
  on.exit(unregister_method("toy_spat"), add = TRUE)

  seen_spread <- NULL
  register_method(
    "toy_spat",
    "balanced",
    sample_fn = function(pik, n = NULL, aux = NULL, spread = NULL, ...) {
      seen_spread <<- spread
      seq_len(n)
    },
    supports_spread = TRUE
  )

  pik <- rep(0.5, 8)
  coords <- data.frame(x = as.double(1:8), y = as.double(8:1))
  s <- balanced_wor(pik, spread = as.matrix(coords), method = "toy_spat")

  expect_s3_class(s, "balanced")
  expect_true(is.matrix(seen_spread))
  expect_equal(dim(seen_spread), c(8L, 2L))
  expect_identical(storage.mode(seen_spread), "double")

  # spread validation fires before sample_fn is called
  expect_error(
    balanced_wor(pik, spread = matrix(1, 3, 2), method = "toy_spat"),
    "does not match"
  )
  expect_error(
    balanced_wor(pik, spread = matrix(NA_real_, 8, 2), method = "toy_spat"),
    "NA, NaN, or Inf"
  )
  expect_error(
    balanced_wor(pik, spread = matrix(0, 8, 0), method = "toy_spat"),
    "at least one column"
  )
})

test_that("strata and spread capabilities compose", {
  on.exit(unregister_method("toy_spat_str"), add = TRUE)

  seen <- list()
  register_method(
    "toy_spat_str",
    "balanced",
    sample_fn = function(
      pik,
      n = NULL,
      aux = NULL,
      strata = NULL,
      spread = NULL,
      ...
    ) {
      seen <<- list(strata = strata, spread = spread)
      unlist(
        lapply(split(seq_along(pik), strata), function(i) {
          i[seq_len(round(sum(pik[i])))]
        }),
        use.names = FALSE
      )
    },
    supports_strata = TRUE,
    supports_spread = TRUE
  )

  pik <- rep(0.5, 8)
  coords <- matrix(as.double(1:16), 8, 2)
  s <- balanced_wor(
    pik,
    strata = rep(1:2, each = 4),
    spread = coords,
    method = "toy_spat_str"
  )

  expect_equal(seen$strata, rep(1:2, each = 4))
  expect_equal(seen$spread, coords)
  expect_true(s$fixed_size)
  expect_equal(s$n, 4L)
})

test_that("strata-capable balanced method receives dense strata labels", {
  on.exit(unregister_method("toy_bal_str"), add = TRUE)

  seen_strata <- NULL
  stratified_sampler <- function(
    pik,
    n = NULL,
    aux = NULL,
    strata = NULL,
    ...
  ) {
    seen_strata <<- strata
    # Take the first round(sum(pik)) units of each stratum
    unlist(
      lapply(split(seq_along(pik), strata), function(i) {
        i[seq_len(round(sum(pik[i])))]
      }),
      use.names = FALSE
    )
  }
  register_method(
    "toy_bal_str",
    "balanced",
    sample_fn = stratified_sampler,
    supports_strata = TRUE
  )

  pik <- rep(0.5, 8)
  s <- balanced_wor(
    pik,
    strata = c(10, 10, 10, 10, 42, 42, 42, 42),
    method = "toy_bal_str"
  )

  expect_equal(seen_strata, rep(1:2, each = 4))
  expect_true(s$fixed_size)
  expect_equal(s$n, 4L)
  expect_equal(s$sample, c(1L, 2L, 5L, 6L))
})

test_that("non-integer per-stratum sums demote fixed_size with warning", {
  on.exit(unregister_method("toy_bal_str"), add = TRUE)

  stratified_sampler <- function(
    pik,
    n = NULL,
    aux = NULL,
    strata = NULL,
    ...
  ) {
    unlist(
      lapply(split(seq_along(pik), strata), function(i) {
        i[seq_len(floor(sum(pik[i])))]
      }),
      use.names = FALSE
    )
  }
  register_method(
    "toy_bal_str",
    "balanced",
    sample_fn = stratified_sampler,
    supports_strata = TRUE
  )

  pik <- rep(0.5, 6)
  expect_warning(
    s <- balanced_wor(pik, strata = rep(1:2, each = 3), method = "toy_bal_str"),
    "not close to an integer"
  )
  expect_false(s$fixed_size)
  expect_equal(s$n, 3)
})

test_that("registered balanced method with joint_fn works", {
  on.exit(unregister_method("toy_bal"), add = TRUE)

  toy_joint <- function(pik, sample_idx = NULL, ...) {
    if (!is.null(sample_idx)) {
      pik <- pik[sample_idx]
    }
    J <- outer(pik, pik)
    diag(J) <- pik
    J
  }
  register_method(
    "toy_bal",
    "balanced",
    sample_fn = toy_balanced_sample,
    joint_fn = toy_joint
  )

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- balanced_wor(pik, method = "toy_bal")

  pikl <- joint_inclusion_prob(s)
  expect_equal(dim(pikl), c(6L, 6L))
  expect_equal(diag(pikl), pik)

  pikl_s <- joint_inclusion_prob(s, sampled_only = TRUE)
  expect_equal(nrow(pikl_s), 2L)
})

test_that("registered balanced method without joint_fn errors", {
  on.exit(unregister_method("toy_bal"), add = TRUE)
  register_method("toy_bal", "balanced", sample_fn = toy_balanced_sample)

  pik <- c(0.2, 0.3, 0.15, 0.35, 0.5, 0.5)
  s <- balanced_wor(pik, method = "toy_bal")
  expect_error(joint_inclusion_prob(s), "not implemented for method 'toy_bal'")
})

# Type gating across dispatchers
test_that("registered methods are rejected by the wrong dispatcher", {
  on.exit(
    {
      unregister_method("gate_wor")
      unregister_method("gate_wr")
      unregister_method("gate_bal")
    },
    add = TRUE
  )

  register_method("gate_wor", "wor", sample_fn = toy_wor_sample)
  register_method(
    "gate_wr",
    "wr",
    sample_fn = function(pik, n = NULL, prn = NULL, ...) {
      sample.int(length(pik), n, replace = TRUE, prob = pik)
    }
  )
  register_method("gate_bal", "balanced", sample_fn = toy_balanced_sample)

  pik <- c(0.5, 0.5, 0.5, 0.5)
  expect_error(
    unequal_prob_wor(pik, method = "gate_wr"),
    "registered as type 'wr'; use unequal_prob_wr"
  )
  expect_error(
    unequal_prob_wor(pik, method = "gate_bal"),
    "registered as type 'balanced'; use balanced_wor"
  )
  expect_error(
    unequal_prob_wr(c(1, 1, 2), method = "gate_wor"),
    "registered as type 'wor'; use unequal_prob_wor"
  )
  expect_error(
    balanced_wor(pik, method = "gate_wor"),
    "registered as type 'wor'"
  )
  expect_error(balanced_wor(pik, method = "gate_wr"), "registered as type 'wr'")
})

test_that("method_spec returns correct metadata for registered balanced method", {
  on.exit(unregister_method("toy_bal"), add = TRUE)
  register_method(
    "toy_bal",
    "balanced",
    sample_fn = toy_balanced_sample,
    supports_strata = TRUE
  )

  spec <- method_spec("toy_bal")
  expect_equal(spec$type, "balanced")
  expect_true(spec$fixed_size)
  expect_false(spec$supports_prn)
  expect_true(spec$supports_aux)
  expect_true(spec$supports_strata)
  expect_false(spec$supports_spread)
})

test_that("method_spec reports supports_spread for spread-capable method", {
  on.exit(unregister_method("toy_spat"), add = TRUE)
  register_method(
    "toy_spat",
    "balanced",
    sample_fn = toy_balanced_sample,
    supports_aux = FALSE,
    supports_spread = TRUE
  )

  spec <- method_spec("toy_spat")
  expect_equal(spec$type, "balanced")
  expect_false(spec$supports_aux)
  expect_false(spec$supports_strata)
  expect_true(spec$supports_spread)
})

test_that("method_spec reports supports_aux = FALSE for wor/wr methods", {
  on.exit(unregister_method("toy_wor"), add = TRUE)
  register_method("toy_wor", "wor", sample_fn = toy_wor_sample)

  expect_false(method_spec("toy_wor")$supports_aux)
  expect_false(method_spec("brewer")$supports_aux)
  expect_true(method_spec("cube")$supports_aux)
})

# Unknown method still errors
test_that("unknown method (not built-in, not registered) still errors", {
  expect_error(unequal_prob_wor(c(0.5, 0.5), method = "nonexistent"))
  expect_error(unequal_prob_wr(c(1, 1), method = "nonexistent"))
})

# variance_family declaration
test_that("register_method validates variance_family values", {
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      variance_family = "banana"
    ),
    "'variance_family' must be NULL or one of"
  )
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      variance_family = 1
    ),
    "'variance_family' must be NULL or one of"
  )
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      variance_family = c("srs", "wr")
    ),
    "'variance_family' must be NULL or one of"
  )
  expect_false(is_registered_method("vf_bad"))
})

test_that("variance_family consistency matrix is enforced", {
  # srs / pps_brewer require fixed_size = TRUE
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      fixed_size = FALSE,
      variance_family = "srs"
    ),
    "requires fixed_size = TRUE"
  )
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      fixed_size = FALSE,
      variance_family = "pps_brewer"
    ),
    "requires fixed_size = TRUE"
  )
  # poisson requires type = "wor" and fixed_size = FALSE
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      fixed_size = TRUE,
      variance_family = "poisson"
    ),
    "fixed_size = FALSE"
  )
  expect_error(
    register_method(
      "vf_bad",
      "wr",
      sample_fn = toy_wr_sample,
      fixed_size = FALSE,
      variance_family = "poisson"
    ),
    'type = "wor"'
  )
  expect_error(
    register_method(
      "vf_bad",
      "balanced",
      sample_fn = toy_wor_sample,
      fixed_size = FALSE,
      variance_family = "poisson"
    ),
    'type = "wor"'
  )
  # wr family requires type = "wr"
  expect_error(
    register_method(
      "vf_bad",
      "wor",
      sample_fn = toy_wor_sample,
      variance_family = "wr"
    ),
    'requires type = "wr"'
  )
  # type = "wr" allows only wr / unsupported
  expect_error(
    register_method(
      "vf_bad",
      "wr",
      sample_fn = toy_wr_sample,
      variance_family = "srs"
    ),
    '"wr" or "unsupported"'
  )
  # balanced allows only pps_brewer / unsupported
  expect_error(
    register_method(
      "vf_bad",
      "balanced",
      sample_fn = toy_wor_sample,
      variance_family = "srs"
    ),
    '"pps_brewer" or "unsupported"'
  )
  expect_false(is_registered_method("vf_bad"))
})

test_that("valid variance_family combinations register and are reported", {
  combos <- list(
    list(name = "vf_srs", type = "wor", fixed = TRUE, family = "srs"),
    list(name = "vf_brw", type = "wor", fixed = TRUE, family = "pps_brewer"),
    list(name = "vf_poi", type = "wor", fixed = FALSE, family = "poisson"),
    list(name = "vf_wun", type = "wor", fixed = FALSE, family = "unsupported"),
    list(name = "vf_wr", type = "wr", fixed = TRUE, family = "wr"),
    list(name = "vf_run", type = "wr", fixed = TRUE, family = "unsupported"),
    list(
      name = "vf_bal",
      type = "balanced",
      fixed = TRUE,
      family = "pps_brewer"
    ),
    list(
      name = "vf_bun",
      type = "balanced",
      fixed = TRUE,
      family = "unsupported"
    )
  )
  for (cmb in combos) {
    on.exit(unregister_method(cmb$name), add = TRUE)
    register_method(
      cmb$name,
      cmb$type,
      sample_fn = toy_wor_sample,
      fixed_size = cmb$fixed,
      variance_family = cmb$family
    )
    expect_identical(
      method_spec(cmb$name)$variance_family,
      cmb$family,
      info = cmb$name
    )
  }
})

test_that("variance_family defaults to NULL (undeclared), field present", {
  on.exit(unregister_method("vf_default"), add = TRUE)
  register_method("vf_default", "wor", sample_fn = toy_wor_sample)

  spec <- method_spec("vf_default")
  expect_true("variance_family" %in% names(spec))
  expect_null(spec$variance_family)
})

test_that("declared variance_family does not affect sampling dispatch", {
  on.exit(unregister_method("vf_draw"), add = TRUE)
  register_method(
    "vf_draw",
    "wor",
    sample_fn = toy_wor_sample,
    variance_family = "pps_brewer"
  )

  pik <- inclusion_prob(c(1, 2, 3, 4, 5), n = 2)
  s <- unequal_prob_wor(pik, method = "vf_draw")
  expect_s3_class(s, "sondage_sample")
  expect_length(s$sample, 2L)
})
