# sondage 0.1.0

Initial CRAN release.

## Sampling methods

* `srs()` - Simple random sampling with or without replacement
* `systematic()` - Systematic sampling with random start
* `bernoulli()` - Bernoulli sampling (random sample size)
* `up_brewer()` - Brewer's method for unequal probability sampling
* `up_maxent()` - Maximum entropy / Conditional Poisson sampling
* `up_systematic()` - Systematic sampling with unequal probabilities
* `up_poisson()` - Poisson sampling (random sample size)
* `up_multinomial()` - Unequal probability sampling with replacement

## Utilities

* `inclusion_prob()` - Compute inclusion probabilities from measure of size

## Features

* All functions return indices for direct use with `df[idx, ]`
* Batch sampling via `up_maxent(pik, nrep = 1000)` for simulations
* C implementations for performance-critical algorithms
