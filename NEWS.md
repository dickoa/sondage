# sondage 0.15.0

Changes since 0.9.1, the previous CRAN release.

## New features

* `unequal_prob_wr()` gains `method = "systematic_mr"`, systematic PPS with
  minimum replacement. A unit larger than the sampling interval is hit once
  per interval it spans, where `unequal_prob_wor(method = "systematic")`
  takes it once, with certainty. Each unit receives the whole part of its
  expected hits or one more, and when every expected hit is below 1 the two
  methods draw the same units from the same seed. `joint_expected_hits()` is
  exact for it.

* `balanced_wor()` gains `prn` for sample coordination with permanent
  random numbers, supported by `method = "scps"` (Grafström & Matei, 2018).
  Units are then visited in row order, so the sample depends only on `pik`,
  `spread`, `prn` and the row order. Drawing a second survey with `1 - prn`
  coordinates the two negatively. Spreading on a measure of response burden
  instead of coordinates gives the adapted SCP sampling of Matei, Smith,
  Smeets & Klingwort (2023).

* `register_method()` accepts `supports_prn = TRUE` for `type = "balanced"`.
  Such methods receive `prn` when the caller supplies it.

* `inclusion_prob()` accepts a fractional `n`, the expected size of a
  random-size design, from 0 up to the number of units with positive size.
  Fixed-size samplers still require probabilities that sum to a whole
  number. A value near a whole number is no longer rounded to it, so
  `inclusion_prob(x, 6.99995)` now sums to 6.99995, not 7.

## Bug fixes

* `joint_inclusion_prob()` and `sampling_cov()` return correct joint
  inclusion probabilities for `method = "sampford"`. With moderate
  populations and samples of about a third to a half of the eligible units,
  cancellation gave wrong values inside the valid bounds, without a warning.
  Queries in that range are also much faster, 1.5 s instead of 6 minutes for
  `N = 1000, n = 500`.

* `joint_inclusion_prob()` and `sampling_cov()` return exact joint inclusion
  probabilities for `method = "cps"` when inclusion probabilities are nearly
  equal. Aires' formula was applied to the target probabilities, which
  calibration meets only to 1e-9, and dividing that residue by the small gap
  between two units' odds gave errors up to 4e-5. Rows then missed the
  fixed-size identity and the covariance had a negative eigenvalue.
  Probabilities within about 2e-7 of each other also got the SRS matrix.
  Close pairs are now computed without subtraction, which also makes them
  fast: `N = 800` with nearly equal probabilities takes 0.06 s instead of
  12.5 s.

* `balanced_wor(method = "cube", condition_aux = TRUE)` keeps the balancing
  equations of `aux`. Conditioning centered the raw columns, which changed
  the constraints when `pik` varies: in a design that balances exactly, 18%
  of draws lost balance. It also dropped a constant column, which is a real
  constraint under unequal `pik`, and any column whose spread fell below
  `qr_tol`, so `1e-8 * aux` lost its constraint. Columns are now centered
  and scaled in the `aux / pik` space the cube balances, within each stratum
  when `strata` is given, and a column is dropped only when its spread is
  small relative to its own size. A column constant within strata, such as
  a region dummy, is now dropped as redundant instead of costing the exact
  balance of the other columns.

* `balanced_wor()` with `strata` treats a stratum sum a small distance from a
  whole number as random-size, with the documented warning. The check
  accepted any sum within 1e-4 of a whole number, so a stratum summing to
  5e-5 was taken as fixed-size and a batch of draws could fail. It now uses
  the floating-point tolerance of the `sum(pik)` check.

* A sample size above `.Machine$integer.max`, from `sum(hits)` or from `N` or
  `n` in the equal-probability functions, is refused with an error. It used
  to become `NA` with a coercion warning and reach the compiled code.

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

* `inclusion_prob()` returns the correct probabilities when the positive
  sizes span 1e200 or more. Once the large units were certain, the small
  ones could get `NaN`, all be made certain so that the probabilities summed
  to more than `n`, or be treated as zero-size: `c(1e-300, 1e-300, 1e10)`
  with `n = 2` gave `1 1 1` instead of `0.5 0.5 1`. Inputs with a smaller
  span give the same results as before, bit for bit.

* `joint_inclusion_prob()` with `method = "cps"` no longer reads uninitialized
  memory when every unit below certainty is so close to 1 that the
  probabilities still sum to their count within the input tolerance. Its
  pair values could exceed the marginal probabilities.

* The high-entropy and Hajek approximations, in `he_jip()`, `hajek_jip()` and
  `joint_inclusion_prob()`, give a pair of near-certainty units the smaller
  of their two probabilities. The value used to depend on the unit order and
  could exceed that bound.

* Arguments passed through `...` now reach a registered `sample_fn` or
  `joint_fn` unchanged. A name such as `q` or `m` could be read as one of
  the package's internal arguments, which gave `joint_fn` a wrong value
  without an error or stopped with a confusing message.

* A balanced method registered as fixed-size now receives a whole `n` when
  non-integer stratum sums make the draw random-size. A floating-point
  residue in `sum(pik)` could give `n = 2.9999999999999991`, which
  `seq_len(n)` truncates to 2.

* `unequal_prob_wor()`, `unequal_prob_wr()` and `balanced_wor()` now stop
  with an error when called with an argument named `n`. These functions take
  the sample size from `sum(pik)` or `sum(hits)`, and R partially matched
  `n = 2` to `nrep`, so the call silently returned two replicates of the
  full sample size.

## Other changes

* `unequal_prob_wr()` requires `sum(hits)` to be a whole number of at least
  1 to floating-point accuracy, as `unequal_prob_wor()` already does for
  `sum(pik)`. It used to accept sums up to 1e-4 away from a whole number
  and round them, so `c(1e-6, 1e-6)` gave an empty sample with `n = 0`.
  Methods registered with `type = "wr"` get the same check.

* `he_jip()` and `hajek_jip()` name the rows and columns of their result
  after `names(pik)`, and label a submatrix of unnamed `pik` with the
  indices in `sample_idx`, as `joint_inclusion_prob()` does. They now
  refuse a result above 10,000 rows, the limit of `joint_inclusion_prob()`,
  and point to `sample_idx`. An argument named `sampled_only`, which they
  used to ignore, is now an error.

* `unequal_prob_wor()` returns `$sample` in ascending unit order for
  `method = "brewer"`, and for `method = "cps"` when some `pik` equal 1, as
  the other built-in without-replacement methods already did. The selected
  units and the random number stream are unchanged.

## Speed

* CPS calibration is about 2.5 times faster and needs O(n) workspace
  instead of O(Nn). Sampford draws that use the non-rejective fallback get
  the same speed-up. This change leaves samples and CPS joint inclusion
  probabilities unchanged, bit for bit.

* `balanced_wor(method = "lpm2")` finds nearest neighbors with a k-d tree
  when `spread` has at most eight columns. A draw from 100,000 units with
  two spreading variables takes 0.06 s instead of 4 s, and one from 50,000
  units with six takes 0.3 s instead of 6.5 s. Samples are unchanged
  unless `spread` has tied distances, as with grid coordinates, duplicated
  points or a 0/1 variable. A seed can then give a different sample of the
  same design.

* `balanced_wor(method = "cube")` computes the direction of each step by
  Gaussian elimination and back substitution instead of a full Gauss-Jordan
  reduction. Draws are about 1.15 times faster with 15 to 30 balancing
  variables and 1.4 times faster with 100. The direction is computed with
  different rounding, so a seed can in rare cases give a different sample.
  None of 600 seeded test draws changed.

* `balanced_wor(method = "scps")` is faster when many units are equally
  distant from the unit being decided, as with a 0/1 burden indicator in
  `spread`. A draw from 17,000 units takes about a tenth of the time. The
  weights those units receive are unchanged, but a seed can give a different
  sample of the same design when many units are tied, because the order in
  which decided units leave the pool changes.

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
  constraints on the realized sample (`bounds = list(B, lower, upper)`,
  Tripet & Tillé 2026). Inequality bounds enable controlled selection in
  the sense of Goodman & Kish: category counts, possibly overlapping (e.g. the
  margins of a two-way control table), are kept within the integers
  adjacent to their expectations while `E(s) = pik` holds exactly.
  They also support controlled matrix rounding and minimum group sizes.
  `"lpm2"` (local pivotal method 2 of Grafström, Lundström & Schelin,
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
* `sampling_cov()`: sampling covariance. `weighted = TRUE` returns
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
  apply, and `method_spec()` reports it for built-in and registered
  methods.
* Registered methods declare where they sit in the first-order
  probability taxonomy with `probabilities`: `"exact"` (realized
  inclusion probabilities, or expected hits for `"wr"`, equal the
  `pik` or `hits` handed to the method), `"approximate"` (honored to a
  documented approximation, as Pareto and sequential Poisson order
  sampling do), or `"unknown"` (the default: `pik` is a selection
  weight only, so design weights `1/pik` would be biased). The
  default is deliberately strict. Downstream packages may refuse to
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
  designs return a matrix, and random-size designs return a list.
* Permanent random numbers (`prn`) for sample coordination (Bernoulli,
  Poisson, SPS, Pareto).
* C implementations for all built-in sampling algorithms.
* Vignette "Extending sondage with Custom Methods".
