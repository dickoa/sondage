# Tests for Brewer's method (up_brewer)

test_that("up_brewer returns correct number of indices", {
    pik <- c(0.2, 0.4, 0.6, 0.8)  # n = 2
    idx <- up_brewer(pik)
    expect_length(idx, 2)
})

test_that("up_brewer indices are in valid range", {
    pik <- c(0.2, 0.4, 0.6, 0.8)
    idx <- up_brewer(pik)
    expect_true(all(idx >= 1 & idx <= 4))
})

test_that("up_brewer has no duplicates", {
    pik <- c(0.2, 0.4, 0.6, 0.8)
    idx <- up_brewer(pik)
    expect_equal(length(unique(idx)), length(idx))
})

test_that("up_brewer achieves correct inclusion probabilities", {
    pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
    N <- length(pik)
    n_sim <- 5000

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- up_brewer(pik)
        counts[idx] <- counts[idx] + 1
    }
    pi_hat <- counts / n_sim

    expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("up_brewer handles certainty selections", {
    # Unit with pik very close to 1 should always be selected
    pik <- c(0.999999, 0.5, 0.5)
    set.seed(123)
    for (i in 1:50) {
        idx <- up_brewer(pik)
        expect_true(1 %in% idx)
    }
})

test_that("up_brewer excludes near-zero units", {
    # Unit with pik very close to 0 should never be selected
    pik <- c(1e-8, 0.5, 0.5)
    set.seed(123)
    for (i in 1:50) {
        idx <- up_brewer(pik)
        expect_false(1 %in% idx)
    }
})

test_that("up_brewer is reproducible with set.seed", {
    pik <- c(0.2, 0.3, 0.5)

    set.seed(999)
    idx1 <- up_brewer(pik)
    set.seed(999)
    idx2 <- up_brewer(pik)

    expect_identical(idx1, idx2)
})

test_that("up_brewer works with larger populations", {
    set.seed(456)
    N <- 100
    n <- 20
    x <- runif(N)
    pik <- n * x / sum(x)

    idx <- up_brewer(pik)
    expect_length(idx, n)
    expect_true(all(idx >= 1 & idx <= N))
})

test_that("up_brewer rejects invalid input", {
    expect_error(up_brewer(c(0.5, NA, 0.5)), "missing values")
    expect_error(up_brewer(c("a", "b")), "numeric")
    expect_error(up_brewer(numeric(0)), "empty")
})
