# Tests for maximum entropy sampling (up_maxent)

test_that("up_maxent returns correct number of indices", {
    pik <- c(0.2, 0.4, 0.6, 0.8)  # n = 2
    idx <- up_maxent(pik)
    expect_length(idx, 2)
})

test_that("up_maxent indices are in valid range", {
    pik <- c(0.2, 0.4, 0.6, 0.8)
    idx <- up_maxent(pik)
    expect_true(all(idx >= 1 & idx <= 4))
})

test_that("up_maxent has no duplicates", {
    pik <- c(0.2, 0.4, 0.6, 0.8)
    idx <- up_maxent(pik)
    expect_equal(length(unique(idx)), length(idx))
})

test_that("up_maxent achieves correct inclusion probabilities", {
    pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
    N <- length(pik)
    n_sim <- 5000

    set.seed(42)
    counts <- integer(N)
    for (i in 1:n_sim) {
        idx <- up_maxent(pik)
        counts[idx] <- counts[idx] + 1
    }
    pi_hat <- counts / n_sim

    expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("up_maxent handles certainty selections", {
    pik <- c(0.999999, 0.5, 0.5)
    set.seed(123)
    for (i in 1:50) {
        idx <- up_maxent(pik)
        expect_true(1 %in% idx)
    }
})

test_that("up_maxent is reproducible with set.seed", {
    pik <- c(0.2, 0.3, 0.5)

    set.seed(999)
    idx1 <- up_maxent(pik)
    set.seed(999)
    idx2 <- up_maxent(pik)

    expect_identical(idx1, idx2)
})

test_that("up_maxent batch mode returns correct dimensions", {
    pik <- c(0.2, 0.4, 0.6, 0.8)  # n = 2
    samples <- up_maxent(pik, nrep = 100)

    expect_true(is.matrix(samples))
    expect_equal(nrow(samples), 2)  # n rows
    expect_equal(ncol(samples), 100)  # nrep columns
})

test_that("up_maxent batch achieves correct inclusion probabilities", {
    pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
    N <- length(pik)
    
    set.seed(42)
    samples <- up_maxent(pik, nrep = 5000)
    
    # Convert indices to indicators and compute means
    counts <- integer(N)
    for (j in 1:ncol(samples)) {
        counts[samples[, j]] <- counts[samples[, j]] + 1
    }
    pi_hat <- counts / ncol(samples)

    expect_equal(pi_hat, pik, tolerance = 0.03)
})

test_that("up_maxent batch is faster than loop", {
    pik <- c(0.07, 0.17, 0.41, 0.61, 0.83, 0.91)
    nrep <- 500
    
    # Batch mode
    t_batch <- system.time({
        set.seed(42)
        samples_batch <- up_maxent(pik, nrep = nrep)
    })[3]
    
    # Loop mode
    t_loop <- system.time({
        set.seed(42)
        samples_loop <- replicate(nrep, up_maxent(pik))
    })[3]
    
    # Batch should be at least 2x faster (usually much more)
    expect_true(t_batch < t_loop / 2 || t_batch < 0.1)  # Allow fast machines
})

test_that("up_maxent works with larger populations", {
    set.seed(123)
    N <- 500
    n <- 100
    pik <- runif(N)
    pik <- pik / sum(pik) * n

    idx <- up_maxent(pik)
    expect_length(idx, n)
    expect_true(all(idx >= 1 & idx <= N))
})

test_that("up_maxent rejects invalid input", {
    expect_error(up_maxent(c(0.5, NA, 0.5)), "missing values")
    expect_error(up_maxent(c("a", "b")), "numeric")
    expect_error(up_maxent(numeric(0)), "empty")
    expect_error(up_maxent(c(0.5, 0.5), nrep = 0), "positive")
})
