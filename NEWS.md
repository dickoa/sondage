# sondage 0.8.6

Initial CRAN release.

## Sampling

Five dispatchers, 15 built-in methods:

* `equal_prob_wor(N, n, method=)`:  `"srs"`, `"systematic"`, `"bernoulli"`.
* `equal_prob_wr(N, n, method=)`:  `"srs"`.
* `unequal_prob_wor(pik, method=)`:  `"cps"` (conditional Poisson /
  maximum entropy), `"brewer"`, `"systematic"`, `"poisson"`, `"sps"`
  (sequential Poisson), `"pareto"`.
* `unequal_prob_wr(hits, method=)`:  `"chromy"` (minimum replacement),
  `"multinomial"`.
* `balanced_wor(pik, aux, strata, spread, bounds, method=)`: `"cube"`
  with optional stratification, and optional linear inequality
  constraints on the realized sample (`bounds = list(B, lower, upper)`;
  Tripet & Tillé 2026). Inequality bounds enable controlled selection à
  la Goodman & Kish: category counts, possibly overlapping (e.g. the
  margins of a two-way control table), are kept within the integers
  adjacent to their expectations while `E(s) = pik` holds exactly.
  They also support controlled matrix rounding and minimum group sizes.
  `"lpm2"` (local pivotal method 2; Grafström, Lundström & Schelin
  2012) draws spatially balanced, well-spread samples on the
  coordinates in `spread`. `"scps"` implements Grafström's (2012)
  maximal-weight spatially correlated Poisson sampling. Its C core
  uses weighted quickselect rather than sorting all remaining units at
  every step, for expected O(N^2 d) time and O(N) workspace.

All sampling functions return S3 design objects with class
`c(prob_class, wor_or_wr, "sondage_sample")` (balanced designs
additionally carry `"balanced"`).

## Design queries

* `inclusion_prob()`: first-order inclusion probabilities (from size
  measures, or extracted from a WOR design).
* `expected_hits()`: expected number of selections (WR analogue).
* `joint_inclusion_prob()`:  exact for `cps`, `systematic`, `poisson`,
  `srs`, `bernoulli`; high-entropy approximation for `brewer`, `sps`,
  `pareto`, `cube`. Not available for `lpm2` or `scps`: well-spread designs are
  deliberately low-entropy, so no tractable approximation applies.
  Their `method_spec()` metadata reports `variance_family = "unsupported"`
  rather than suggesting a high-entropy PPS variance treatment.
* `joint_expected_hits()`: exact analytic for `multinomial` / `srs`,
  simulation-based for `chromy`.
* `sampling_cov()`: sampling covariance; `weighted = TRUE` returns
  Sen-Yates-Grundy check quantities.

The matrix-valued generics accept `sampled_only = TRUE` to return only the
sampled-units submatrix (useful for large populations).

## Extensibility

* `register_method()` / `unregister_method()` / `registered_methods()` /
  `is_registered_method()` / `method_spec()` register custom
  unequal-probability and balanced methods that flow through the
  existing dispatchers and generics. Balanced methods (`type =
  "balanced"`) dispatch through `balanced_wor()` and opt into
  stratification with `supports_strata = TRUE` or spatial spreading
  with `supports_spread = TRUE` (well-spread designs such as the
  local pivotal method, SCPS, or the local cube receive the
  coordinate matrix passed to `balanced_wor(spread = )`), the same
  way WOR/WR methods opt into permanent random numbers with
  `supports_prn`. Spread-only methods can declare
  `supports_aux = FALSE` so that passing `aux` errors instead of
  being silently ignored.
* Registered methods can declare a `variance_family` (`"srs"`,
  `"pps_brewer"`, `"poisson"`, `"wr"`, `"unsupported"`) describing
  the design-based variance treatment downstream packages should
  apply; `method_spec()` reports it for built-in and registered
  methods.
* `he_jip()` (Brewer & Donadio 2003 high-entropy approximation) and
  `hajek_jip()` (Hajek 1964) are exported and can be passed directly
  as `joint_fn` to `register_method()`.

## Other features

* Batch sampling via `nrep` for Monte Carlo simulations. Fixed-size
  designs return a matrix; random-size designs return a list.
* Permanent random numbers (`prn`) for sample coordination (Bernoulli,
  Poisson, SPS, Pareto).
* C implementations for all built-in sampling algorithms.
* Vignette "Extending sondage with Custom Methods".
