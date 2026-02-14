# sondage 0.1.0

Initial CRAN release.

## Sampling dispatchers

Four dispatchers returning S3 design objects with class
`c("<prob_class>", "<replacement>", "sondage_sample")`:

* `equal_prob_wor(N, n, method=)` - SRS, systematic, Bernoulli
* `equal_prob_wr(N, n, method=)` - SRS with replacement
* `unequal_prob_wor(pik, method=)` - CPS (max entropy), Brewer, systematic PPS, Poisson
* `unequal_prob_wr(hits, method=)` - Chromy (minimum replacement), multinomial PPS

All dispatchers support `nrep` for batch sampling and `u` for sample coordination.

## Generics

* `inclusion_prob()` - Compute or extract first-order inclusion probabilities
* `expected_hits()` - Compute or extract expected hits (WR analogue)
* `joint_inclusion_prob()` - Joint inclusion probabilities (WOR)
* `joint_expected_hits()` - Pairwise expectations E(n_i n_j) (WR)
* `sampling_cov()` - Sampling covariance matrix (WOR and WR)
* `sampling_cov(scaled = TRUE)` - SYG check quantities

## Features

* All functions return design objects for direct use with `s$sample`
* Batch sampling via `nrep` argument for simulations
* C implementations for performance-critical algorithms
* Unified variance estimation framework following Chromy (2009)
