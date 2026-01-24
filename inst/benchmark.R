#!/usr/bin/env Rscript
#
# sondage Package Benchmark
# ══════════════════════════════════════════════════════════════════════════════
#
# Compares performance of sondage functions vs sampling package equivalents.
#
# Usage:
#   Rscript benchmark.R
#   # or source("benchmark.R") in R
#

cat("\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
cat("  sondage Package Benchmark\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

# Setup
library(sondage)
has_sampling <- requireNamespace("sampling", quietly = TRUE)
if (has_sampling) {
    library(sampling)
} else {
    cat("Note: 'sampling' package not installed. Showing sondage times only.\n\n")
}

# Helper for timing
time_func <- function(expr, reps = 50) {
    times <- replicate(reps, system.time(expr)["elapsed"])
    c(mean = mean(times) * 1000, sd = sd(times) * 1000)
}

# Generate test inclusion probabilities
make_pik <- function(N, n, seed = 42) {
    set.seed(seed)
    pik <- runif(N, 0.01, 0.99)
    pik <- pik / sum(pik) * n
    pmax(pmin(pik, 0.999), 0.001)
}

# ──────────────────────────────────────────────────────────────────────────────
# Benchmark 1: Brewer's Method
# ──────────────────────────────────────────────────────────────────────────────

cat("BENCHMARK 1: Brewer's Method\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

for (N in c(100, 500, 1000)) {
    n <- N / 5
    pik <- make_pik(N, n)
    
    t_sondage <- time_func(up_brewer(pik))
    
    if (has_sampling) {
        t_sampling <- time_func(UPbrewer(pik))
        speedup <- t_sampling["mean"] / t_sondage["mean"]
        cat(sprintf("N=%4d: sondage %6.2f ms, sampling %6.2f ms (%.1fx faster)\n",
                    N, t_sondage["mean"], t_sampling["mean"], speedup))
    } else {
        cat(sprintf("N=%4d: sondage %6.2f ms\n", N, t_sondage["mean"]))
    }
}

# ──────────────────────────────────────────────────────────────────────────────
# Benchmark 2: Maximum Entropy (CPS)
# ──────────────────────────────────────────────────────────────────────────────

cat("\nBENCHMARK 2: Maximum Entropy Sampling\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

for (N in c(100, 500, 1000)) {
    n <- N / 5
    pik <- make_pik(N, n)
    
    t_sondage <- time_func(up_maxent(pik))
    
    if (has_sampling) {
        t_sampling <- time_func(UPmaxentropy(pik))
        speedup <- t_sampling["mean"] / t_sondage["mean"]
        cat(sprintf("N=%4d: sondage %6.2f ms, sampling %6.2f ms (%.1fx faster)\n",
                    N, t_sondage["mean"], t_sampling["mean"], speedup))
    } else {
        cat(sprintf("N=%4d: sondage %6.2f ms\n", N, t_sondage["mean"]))
    }
}

# ──────────────────────────────────────────────────────────────────────────────
# Benchmark 3: Batch Mode Performance
# ──────────────────────────────────────────────────────────────────────────────

cat("\nBENCHMARK 3: Repeated Sampling (batch mode vs loop)\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

N <- 500
n <- 100
nrep <- 100
pik <- make_pik(N, n)

# Batch mode
t_batch <- system.time({
    set.seed(42)
    samples_batch <- up_maxent(pik, nrep = nrep)
})[3] * 1000

# Loop mode
t_loop <- system.time({
    set.seed(42)
    samples_loop <- replicate(nrep, up_maxent(pik))
})[3] * 1000

cat(sprintf("N=%d, n=%d, nrep=%d:\n", N, n, nrep))
cat(sprintf("  Batch mode: %6.1f ms\n", t_batch))
cat(sprintf("  Loop mode:  %6.1f ms\n", t_loop))
cat(sprintf("  Speedup:    %.1fx\n", t_loop / t_batch))

# ──────────────────────────────────────────────────────────────────────────────
# Benchmark 4: Equal Probability Methods
# ──────────────────────────────────────────────────────────────────────────────

cat("\nBENCHMARK 4: Equal Probability Methods (N=100,000)\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

N <- 100000
n <- 1000

t_srs <- time_func(srs(n, N), reps = 100)
t_sys <- time_func(sys(n, N), reps = 100)
t_bern <- time_func(bern(n/N, N), reps = 100)

cat(sprintf("srs():  %6.2f ms\n", t_srs["mean"]))
cat(sprintf("sys():  %6.2f ms\n", t_sys["mean"]))
cat(sprintf("bern(): %6.2f ms\n", t_bern["mean"]))

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  Benchmark complete\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")
