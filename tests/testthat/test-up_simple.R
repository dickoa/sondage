# Tests for up_poisson, up_multinomial, up_systematic

# ============================================================
# up_poisson tests
# ============================================================

test_that("up_poisson returns indices in valid range", {
    pik <- c(0.2, 0.5, 0.8)
    idx <- up_poisson(pik)
    expect_true(all(idx >= 1 & idx <= 3))
})

test_that("up_poisson has no duplicates", {
    pik <- c(0.2, 0.5, 0.8)
    idx <- up_poisson(pik)
    expect_equal(length(unique(idx)), length(idx))
})

test_that("up_poisson achieves correct inclusion probabilities", {
    pik <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    N <- length(pik)
    n_sim <- 5000

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- up_poisson(pik)
        counts[idx] <- counts[idx] + 1
    }
    pi_hat <- counts / n_sim

    expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("up_poisson handles boundary cases", {
    # All 0 - select nothing
    idx <- up_poisson(c(0, 0, 0))
    expect_length(idx, 0)

    # All 1 - select everything
    idx <- up_poisson(c(1, 1, 1))
    expect_equal(sort(idx), 1:3)
})

test_that("up_poisson is reproducible with set.seed", {
    pik <- c(0.2, 0.5, 0.8)

    set.seed(123)
    idx1 <- up_poisson(pik)
    set.seed(123)
    idx2 <- up_poisson(pik)

    expect_identical(idx1, idx2)
})

test_that("up_poisson rejects invalid input", {
    expect_error(up_poisson(c(0.5, NA)), "missing values")
    expect_error(up_poisson(c(-0.1, 0.5)), "must be in")
    expect_error(up_poisson(c(0.5, 1.1)), "must be in")
})

# ============================================================
# up_multinomial tests
# ============================================================

test_that("up_multinomial returns correct number of indices", {
    pik <- c(0.2, 0.4, 0.6, 0.8)  # sum = 2
    idx <- up_multinomial(pik)
    expect_length(idx, 2)
})

test_that("up_multinomial indices are in valid range", {
    pik <- c(1, 2, 3, 4)  # sum = 10
    idx <- up_multinomial(pik)
    expect_true(all(idx >= 1 & idx <= 4))
})

test_that("up_multinomial can have duplicates", {
    set.seed(42)
    pik <- c(1, 2, 3, 4)  # sum = 10 draws
    idx <- up_multinomial(pik)
    # With 10 draws from 4 units, we expect duplicates
    expect_true(length(unique(idx)) < 10)
})

test_that("up_multinomial achieves correct proportions", {
    pik <- c(1, 2, 3, 4)  # Should get 10%, 20%, 30%, 40%
    n_sim <- 5000
    n_draws <- round(sum(pik))  # 10
    N <- length(pik)

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- up_multinomial(pik)
        tab <- tabulate(idx, nbins = N)
        counts <- counts + tab
    }
    proportions <- counts / (n_sim * n_draws)
    expected <- pik / sum(pik)

    expect_equal(proportions, expected, tolerance = 0.02)
})

test_that("up_multinomial is reproducible with set.seed", {
    pik <- c(1, 2, 3, 4)

    set.seed(123)
    idx1 <- up_multinomial(pik)
    set.seed(123)
    idx2 <- up_multinomial(pik)

    expect_identical(idx1, idx2)
})

test_that("up_multinomial rejects invalid input", {
    expect_error(up_multinomial(c(1, NA, 2)), "missing values")
    expect_error(up_multinomial(c(1, -1, 2)), "non-negative")
    expect_error(up_multinomial(c(0, 0, 0)), "sum of pik must be positive")
})

# ============================================================
# up_systematic tests
# ============================================================

test_that("up_systematic returns correct number of indices", {
    pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)  # n = 2
    idx <- up_systematic(pik)
    expect_length(idx, 2)
})

test_that("up_systematic indices are in valid range", {
    pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
    idx <- up_systematic(pik)
    expect_true(all(idx >= 1 & idx <= 5))
})

test_that("up_systematic has no duplicates", {
    pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)
    idx <- up_systematic(pik)
    expect_equal(length(unique(idx)), length(idx))
})

test_that("up_systematic achieves correct inclusion probabilities", {
    pik <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5)  # n = 2
    N <- length(pik)
    n_sim <- 10000

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- up_systematic(pik)
        counts[idx] <- counts[idx] + 1
    }
    pi_hat <- counts / n_sim

    expect_equal(pi_hat, pik, tolerance = 0.02)
})

test_that("up_systematic handles certainty selections", {
    # pik near 1 should always be selected
    pik <- c(0.999999, 0.5, 0.5)
    set.seed(123)
    for (i in 1:50) {
        idx <- up_systematic(pik)
        expect_true(1 %in% idx)
    }
})

test_that("up_systematic is reproducible with set.seed", {
    pik <- c(0.2, 0.3, 0.5, 0.4, 0.6)

    set.seed(123)
    idx1 <- up_systematic(pik)
    set.seed(123)
    idx2 <- up_systematic(pik)

    expect_identical(idx1, idx2)
})

test_that("up_systematic works with larger populations", {
    set.seed(42)
    N <- 500
    n <- 100
    pik <- runif(N)
    pik <- pik / sum(pik) * n

    idx <- up_systematic(pik)

    expect_length(idx, n)
    expect_true(all(idx >= 1 & idx <= N))
})

test_that("up_systematic rejects invalid input", {
    expect_error(up_systematic(c(0.5, NA)), "missing values")
    expect_error(up_systematic(c("a", "b")), "numeric")
    expect_error(up_systematic(numeric(0)), "empty")
})
