# sondage 0.9.0

Initial CRAN release.

## Sampling

Five dispatchers, 16 built-in methods:

* `equal_prob_wor(N, n, method=)`:  `"srs"`, `"systematic"`, `"bernoulli"`.
* `equal_prob_wr(N, n, method=)`:  `"srs"`.
* `unequal_prob_wor(pik, method=)`:  `"cps"` (conditional Poisson /
  maximum entropy), `"sampford"` (exact fixed-size PPS with exact joint
  inclusion probabilities), `"brewer"`, `"systematic"`, `"poisson"`, `"sps"`
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
`c(prob_class, {wor|wr}, "sondage_sample")` (balanced designs
additionally carry `"balanced"`).

## Design queries

* `inclusion_prob()`: first-order inclusion probabilities (from size
  measures, or extracted from a WOR design).
* `expected_hits()`: expected number of selections (WR analogue).
* `joint_inclusion_prob()`:  exact for `cps`, `sampford`, `systematic`, `poisson`,
  `srs`, `bernoulli`; high-entropy approximation for `brewer`, `sps`,
  `pareto`, `cube`. Not available for `lpm2` or `scps`: well-spread designs are
  deliberately low-entropy, so no tractable approximation applies.
  Their `method_spec()` metadata reports `variance_family = "unsupported"`
  rather than suggesting a high-entropy PPS variance treatment.
* `joint_expected_hits()`: exact analytic for `multinomial` / `srs`,
  simulation-based for `chromy`.
* `sampling_cov()`: sampling covariance; `weighted = TRUE` returns
  Sen-Yates-Grundy check quantities.

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
* `register_method()` now rejects an already registered custom name instead
  of silently replacing its implementation. Deliberate replacements require
  an explicit call to `unregister_method()` first.
* Registered methods can declare a `variance_family` (`"srs"`,
  `"pps_brewer"`, `"poisson"`, `"wr"`, `"unsupported"`) describing
  the design-based variance treatment downstream packages should
  apply; `method_spec()` reports it for built-in and registered
  methods.
* Registered methods declare where they sit in the first-order
  probability taxonomy with `probabilities`: `"exact"` (realized
  inclusion probabilities, or expected hits for `"wr"`, equal the
  `pik` or `hits` handed to the method), `"approximate"` (honored to a
  documented approximation, as Pareto and sequential Poisson order
  sampling do), or `"unknown"` (the default: `pik` is a selection
  weight only, so design weights `1/pik` would be biased). The
  default is deliberately strict; downstream packages may refuse to
  draw with an `"unknown"` method, while sampling through sondage
  itself is never affected. `method_spec()` reports the tier for
  every method: built-ins are `"exact"` except `"sps"` and
  `"pareto"`, which report `"approximate"`.
* Custom WR callback contracts are documented with `hits`, consistently
  with the values actually passed to `sample_fn` and `joint_fn`. Validation
  errors now distinguish joint expected hits from joint inclusion
  probabilities.
* Capability arguments in `register_method()` now default to `NULL`, meaning
  unspecified. Explicit capabilities are type-specific: WOR/WR methods may
  declare `supports_prn`, while balanced methods may declare `supports_aux`,
  `supports_strata`, and `supports_spread`. Supplying an irrelevant capability
  now errors instead of being silently normalized.
* `method_spec()` also returns `sample_fn` and `joint_fn`, the
  implementation functions of a registered method (`NULL` for
  built-ins). Downstream packages use them to fingerprint the
  implementation a saved design was executed with.
* `method_spec()` now identifies the public `dispatcher` for every method.
  Shared built-in names such as `"srs"` and `"systematic"` require an explicit
  dispatcher instead of silently selecting one variant.
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
