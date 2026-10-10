# sondage 0.11.0

* `inclusion_prob()` accepts a fractional `n`, the expected size of a
  random-size design, from 0 up to the number of units with positive size.
  Fixed-size samplers still require probabilities that sum to a whole
  number. A value near a whole number is no longer rounded to it, so
  `inclusion_prob(x, 6.99995)` now sums to 6.99995, not 7.

* `inclusion_prob()` returns the correct probabilities when the positive
  sizes span 1e200 or more. Once the large units were certain, the small
  ones could get `NaN`, all be made certain so that the probabilities summed
  to more than `n`, or be treated as zero-size: `c(1e-300, 1e-300, 1e10)`
  with `n = 2` gave `1 1 1` instead of `0.5 0.5 1`. Inputs with a smaller
  span give the same results as before, bit for bit.

# sondage 0.10.0

* `balanced_wor()` gains `prn` for sample coordination with permanent
  random numbers, supported by `method = "scps"` (Grafström & Matei, 2018).
  Units are then visited in row order, so the sample depends only on `pik`,
  `spread`, `prn` and the row order. Drawing a second survey with `1 - prn`
  coordinates the two negatively. Spreading on a measure of response burden
  instead of coordinates gives the adapted SCP sampling of Matei, Smith,
  Smeets & Klingwort (2023).

* `balanced_wor(method = "scps")` is faster when many units are equally
  distant from the unit being decided, as with a 0/1 burden indicator in
  `spread`. A draw from 17,000 units takes about a tenth of the time. The
  weights those units receive are unchanged. Samples drawn with `prn` are
  unchanged. Without `prn`, a seed can give a different sample of the same
  design when many units are tied, because the order in which decided units
  leave the pool changes.

* `register_method()` accepts `supports_prn = TRUE` for `type = "balanced"`.
  Such methods receive `prn` when the caller supplies it.

* `balanced_wor(method = "scps")` no longer fails with "SCPS maximal weights
  are numerically infeasible" when `sum(pik)` misses an integer by a residue
  the input check accepts. At 1e-11 every draw failed, and residues of a few
  ulps, as `inclusion_prob()` can return, made some draws fail. Draws that
  succeeded before are unchanged.

* `unequal_prob_wor(method = "cps")` now draws a single unit, or all but one
  unit, from the exact design. Calibration did not converge when two units
  had inclusion probabilities near 0.5, so it warned after 500 iterations and
  the realized probabilities missed their targets by up to 7e-3. One draw has
  a closed form, odds proportional to `pik`, which is now used. Samples of
  every other size are unchanged.

# sondage 0.9.1

* Allocated the `long double` scratch buffers used by the exact Sampford
  joint inclusion probabilities with `R_allocLD()` instead of `R_alloc()`.
  `R_alloc()` guarantees only the alignment required by `double`, so on
  platforms where `long double` needs 16-byte alignment the accesses were
  undefined behavior, as reported by CRAN's gcc-UBSan check. Computed
  values are unchanged.

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
