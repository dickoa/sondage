# sondage 0.8

Initial CRAN release.

## Sampling

Five dispatchers, 13 built-in methods:

* `equal_prob_wor(N, n, method=)`:  `"srs"`, `"systematic"`, `"bernoulli"`.
* `equal_prob_wr(N, n, method=)`:  `"srs"`.
* `unequal_prob_wor(pik, method=)`:  `"cps"` (conditional Poisson /
  maximum entropy), `"brewer"`, `"systematic"`, `"poisson"`, `"sps"`
  (sequential Poisson), `"pareto"`.
* `unequal_prob_wr(hits, method=)`:  `"chromy"` (minimum replacement),
  `"multinomial"`.
* `balanced_wor(pik, aux, strata, method=)`: `"cube"` with optional
  stratification.

All sampling functions return S3 design objects with class
`c(prob_class, wor_or_wr, "sondage_sample")` (cube designs additionally
carry `"balanced"`).

## Design queries

* `inclusion_prob()`: first-order inclusion probabilities (from size
  measures, or extracted from a WOR design).
* `expected_hits()`: expected number of selections (WR analogue).
* `joint_inclusion_prob()`:  exact for `cps`, `systematic`, `poisson`,
  `srs`, `bernoulli`; high-entropy approximation for `brewer`, `sps`,
  `pareto`, `cube`.
* `joint_expected_hits()`: exact analytic for `multinomial` / `srs`,
  simulation-based for `chromy`.
* `sampling_cov()`: sampling covariance; `weighted = TRUE` returns
  Sen-Yates-Grundy check quantities.

All generics accept `sampled_only = TRUE` to return only the
sampled-units submatrix (useful for large populations).

## Extensibility

* `register_method()` / `unregister_method()` / `registered_methods()` /
  `is_registered_method()` / `method_spec()` register custom
  unequal-probability methods that flow through the existing
  dispatchers and generics.
* `he_jip()` (Brewer & Donadio 2003 high-entropy approximation) and
  `hajek_jip()` (Hajek 1964) are exported and can be passed directly
  as `joint_fn` to `register_method()`.

## Other features

* Batch sampling via `nrep` for Monte Carlo simulations. Fixed-size
  designs return a matrix; random-size designs return a list.
* Permanent random numbers (`prn`) for sample coordination (Bernoulli,
  Poisson, SPS, Pareto).
* C implementations for all sampling algorithms.
* Vignette "Extending sondage with Custom Methods".
