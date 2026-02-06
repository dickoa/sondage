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
* `up_chromy()` - Chromy's sequential PPS with minimum replacement

## Joint inclusion probabilities
* `up_maxent_jip()` - Exact CPS joint probabilities (Aires' formula)
* `up_brewer_jip()` - Brewer & Donadio approximation
* `up_systematic_jip()` - Exact systematic joint probabilities
* `up_poisson_jip()` - Independent selections

## Utilities
* `inclusion_prob()` - Compute inclusion probabilities from measure of size

## Features
* All functions return indices for direct use with `df[idx, ]`
* Batch sampling via `nrep` argument for simulations (maxent/cps)
* C implementations for performance-critical algorithms
