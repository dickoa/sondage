# sondage 0.2.9999

## Sampling dispatchers

Four dispatchers returning S3 design objects with class
`c("<prob_class>", "<replacement>", "sondage_sample")`:

* `equal_prob_wor(N, n, method=)` - SRS, systematic, Bernoulli
* `equal_prob_wr(N, n, method=)` - SRS with replacement
* `unequal_prob_wor(pik, method=)` - CPS (max entropy), Brewer, systematic PPS, Poisson, SPS, Pareto
* `unequal_prob_wr(hits, method=)` - Chromy (minimum replacement), multinomial PPS

All dispatchers support `nrep` for batch sampling and `prn` for sample coordination
(Bernoulli, Poisson, SPS, Pareto).

## Generics

* `inclusion_prob()` - Compute or extract first-order inclusion probabilities
* `expected_hits()` - Compute or extract expected hits (WR analogue)
* `joint_inclusion_prob()` - Joint inclusion probabilities (WOR)
* `joint_expected_hits()` - Pairwise expectations E(n_i n_j) (WR)
* `sampling_cov()` - Sampling covariance matrix (WOR and WR)
* `sampling_cov(weighted = TRUE)` - SYG check quantities

## Features

* All functions return design objects for direct use with `s$sample`
* Batch sampling via `nrep` argument for simulations
* C implementations for performance-critical algorithms
* Unified variance estimation framework following Chromy (2009)
* High-entropy joint inclusion probabilities (Brewer & Donadio, 2003)
  handle certainty units correctly and clamp to valid bounds;
  a warning is issued when the marginal defect exceeds 5% of n
* Equal-probability systematic sampling returns exact joint inclusion
  probabilities (with structural zeros), not the SRS approximation
* N-size guard on dense joint probability matrices (N > 10,000)
* `joint_expected_hits()` diagonal is E(n_i^2), consistent with
  the WOR convention where diagonal is pi_i
* `sampling_cov(weighted = TRUE)` returns `NA` (not `NaN`) for
  undefined entries (zero joint probabilities or zero selection
  probabilities)
* `print()` handles non-integer expected sample sizes
  (Poisson, Bernoulli)
