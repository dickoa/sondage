# Tests for equal probability sampling: srs(), sys)(, bernoulli()

# ============================================================
# srs() tests
# ============================================================

test_that("srs returns correct number of indices", {
    idx <- srs(5, 100)
    expect_length(idx, 5)
})

test_that("srs indices are in valid range", {
    idx <- srs(10, 50)
    expect_true(all(idx >= 1 & idx <= 50))
})

test_that("srs without replacement has no duplicates", {
    idx <- srs(20, 100)
    expect_equal(length(unique(idx)), 20)
})

test_that("srs with replacement can have duplicates", {
    set.seed(42)
    # With high probability, 50 draws from 10 will have duplicates
    idx <- srs(50, 10, replace = TRUE)
    expect_true(length(unique(idx)) < 50)
})

test_that("srs achieves uniform selection", {
    n_sim <- 5000
    N <- 10
    n <- 3

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- srs(n, N)
        counts[idx] <- counts[idx] + 1
    }

    # Each unit should be selected about n_sim * n/N times
    expected <- n_sim * n / N
    expect_true(all(abs(counts - expected) < 0.1 * expected))
})

test_that("srs is reproducible with set.seed", {
    set.seed(123)
    idx1 <- srs(5, 100)
    set.seed(123)
    idx2 <- srs(5, 100)
    expect_identical(idx1, idx2)
})

test_that("srs rejects invalid input", {
    expect_error(srs(10, 5), "cannot exceed")  # n > N without replacement
    expect_error(srs(-1, 10), "non-negative")
    expect_error(srs(5, 0), "positive")
})

test_that("srs handles edge cases", {
    # n = 0
    idx <- srs(0, 10)
    expect_length(idx, 0)

    # n = N (select all)
    idx <- srs(5, 5)
    expect_equal(sort(idx), 1:5)
})

# ============================================================
# systematic() tests
# ============================================================

test_that("sys returns correct number of indices", {
    idx <- systematic(5, 100)
    expect_length(idx, 5)
})

test_that("sys indices are in valid range", {
    idx <- systematic(10, 50)
    expect_true(all(idx >= 1 & idx <= 50))
})

test_that("sys has no duplicates", {
    idx <- systematic(20, 100)
    expect_equal(length(unique(idx)), 20)
})

test_that("sys produces evenly spaced indices", {
    idx <- systematic(5, 100)  # interval = 20
    diffs <- diff(idx)
    # All differences should be close to 20
    expect_true(all(diffs >= 19 & diffs <= 21))
})

test_that("sys achieves uniform selection", {
    n_sim <- 5000
    N <- 12
    n <- 3

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- systematic(n, N)
        counts[idx] <- counts[idx] + 1
    }

    expected <- n_sim * n / N
    expect_true(all(abs(counts - expected) < 0.15 * expected))
})

test_that("sys is reproducible with set.seed", {
    set.seed(123)
    idx1 <- systematic(5, 100)
    set.seed(123)
    idx2 <- systematic(5, 100)
    expect_identical(idx1, idx2)
})

test_that("sys rejects invalid input", {
    expect_error(systematic(10, 5), "cannot exceed")
    expect_error(systematic(-1, 10), "non-negative")
})

# ============================================================
# bernoulli() tests
# ============================================================

test_that("bern returns indices in valid range", {
    idx <- bernoulli(0.5, 100)
    expect_true(all(idx >= 1 & idx <= 100))
})

test_that("bern has no duplicates", {
    idx <- bernoulli(0.5, 100)
    expect_equal(length(unique(idx)), length(idx))
})

test_that("bern achieves correct expected size", {
    n_sim <- 1000
    N <- 100
    p <- 0.3

    set.seed(42)
    sizes <- replicate(n_sim, length(bernoulli(p, N)))

    expect_equal(mean(sizes), N * p, tolerance = 2)
})

test_that("bern achieves correct inclusion probability", {
    n_sim <- 5000
    N <- 20
    p <- 0.4

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- bernoulli(p, N)
        counts[idx] <- counts[idx] + 1
    }

    pi_hat <- counts / n_sim
    expect_equal(pi_hat, rep(p, N), tolerance = 0.03)
})

test_that("bern is reproducible with set.seed", {
    set.seed(123)
    idx1 <- bernoulli(0.5, 100)
    set.seed(123)
    idx2 <- bernoulli(0.5, 100)
    expect_identical(idx1, idx2)
})

test_that("bern handles boundary cases", {
    # p = 0 (select nothing)
    idx <- bernoulli(0, 100)
    expect_length(idx, 0)

    # p = 1 (select all)
    idx <- bernoulli(1, 100)
    expect_equal(sort(idx), 1:100)
})

test_that("bern rejects invalid input", {
    expect_error(bernoulli(-0.1, 100), "probability")
    expect_error(bernoulli(1.1, 100), "probability")
})
