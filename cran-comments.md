## Purpose of this submission

This release follows 0.9.1 on CRAN. It fixes defects present in 0.9.1, adds
three features, and speeds up several samplers. The defects are listed first.

* `joint_inclusion_prob()` and `sampling_cov()` returned wrong joint inclusion
  probabilities for `unequal_prob_wor(method = "sampford")` with moderate
  populations and samples of about a third to a half of the eligible units.
  Cancellation in a polynomial deletion gave values inside the valid bounds,
  so no warning was issued. The computation no longer cancels, and it agrees
  with a 500-digit evaluation of the design definition.
* `joint_inclusion_prob()` and `sampling_cov()` returned joint inclusion
  probabilities for `method = "cps"` that were off by up to 4e-5 when
  inclusion probabilities were nearly equal, so rows missed the fixed-size
  identity and the covariance had a negative eigenvalue. Aires' formula
  divided the calibration residue by the small gap between two units' odds.
  It now uses the probabilities the calibrated odds imply, and close pairs
  are computed without subtraction. The result agrees with exact
  enumeration to 2e-11.
* `balanced_wor(method = "cube", condition_aux = TRUE)` changed the
  balancing equations when `pik` varied, so a design that balances exactly
  lost balance in 18% of draws. Conditioning now works in the `aux / pik`
  space the cube balances.
* `balanced_wor()` with `strata` took a stratum sum within 1e-4 of a whole
  number as fixed-size, so a batch of draws could fail. It now uses the
  floating-point tolerance of the `sum(pik)` check.
* A sample size above `.Machine$integer.max` became `NA` with a coercion
  warning and reached the compiled code. It is now refused.
* `balanced_wor(method = "scps")` failed with "SCPS maximal weights are
  numerically infeasible" whenever `sum(pik)` missed an integer by a residue
  the input check accepts. `inclusion_prob()` output routinely carries such
  residues, so draws from realistic strata failed, in some cases on every
  draw. The tolerance now scales with the residue.
* `unequal_prob_wor(method = "cps")` did not converge when drawing one unit
  (or all but one) with two inclusion probabilities near 0.5, warned, and
  drew off target by up to 7e-3. One draw now uses its exact closed form.
* `inclusion_prob()` could return `NaN`, or probabilities summing to more
  than `n`, when the positive sizes spanned 1e200 or more, because the small
  units became subnormal once the large ones were certain. Their values are
  now recomputed from the raw sizes in that case, and the result is checked
  to sum to `n` before it is returned. Inputs with a smaller span give the
  same results as before, bit for bit.
* `unequal_prob_wr()` accepted sums up to 1e-4 away from a whole number and
  rounded them, so `c(1e-6, 1e-6)` gave an empty sample. It now requires
  `sum(hits)` to be a whole number of at least 1 to floating-point accuracy,
  as `unequal_prob_wor()` already did for `sum(pik)`.
* Arguments passed through `...` could be misread on their way to a method
  registered with `register_method()`, which gave a wrong value without an
  error or stopped with a confusing message.
* `joint_inclusion_prob()` with `method = "cps"` read uninitialized memory
  when every unit below certainty was so close to 1 that the probabilities
  still summed to their count within the input tolerance.
* The high-entropy and Hajek approximations gave a pair of near-certainty
  units a value that depended on the unit order and could exceed the smaller
  marginal probability.

Three features are new.

* `unequal_prob_wr()` gains `method = "systematic_mr"`, systematic PPS with
  minimum replacement, with exact joint expected hits.
* `balanced_wor()` gains `prn` for coordinating samples with permanent
  random numbers, supported by `method = "scps"`, and registered balanced
  methods may declare PRN support.
* `inclusion_prob()` accepts a fractional `n`, the expected size of a
  random-size design such as Poisson sampling, from 0 up to the number of
  units with positive size. Fixed-size samplers still refuse probabilities
  that do not sum to a whole number. A value near a whole number is no
  longer rounded to it.

The remaining changes concern speed and output order.

* CPS calibration is about 2.5 times faster and uses O(n) workspace instead
  of O(Nn). Sampford draws that use the non-rejective fallback get the same
  speed-up. The speed-up leaves samples as they were.
* `balanced_wor(method = "lpm2")` uses a k-d tree when `spread` has at most
  eight columns, `method = "cube"` uses Gaussian elimination and back
  substitution, and `method = "scps"` is faster when many units share a
  spreading value. A seed can give a different sample of the same design for
  `lpm2` and `scps` when many units are tied and, in rare cases, for `cube`.
* `unequal_prob_wor()` returns `$sample` in ascending unit order for
  `method = "brewer"`, and for `method = "cps"` when some `pik` equal 1, as
  its other built-in methods do. The selected units are unchanged.
* `he_jip()` and `hajek_jip()` label their result as
  `joint_inclusion_prob()` does, and refuse a result above 10,000 rows and
  an argument named `sampled_only`.

The compiled code changed in `src/cps.c`, `src/cps_core.h`, `src/cube.c`,
`src/inclusion_prob.c`, `src/init.c`, `src/joint_probs.c`,
`src/kdtree_core.h`, `src/lpm.c`, `src/scps.c`, `src/spatial_core.h`,
`src/up_brewer.c`, `src/up_order.c` and `src/up_sampford.c`, and
`src/systematic_core.h` is new. The full test suite and a stress run of the new code
paths are clean under gcc `-fsanitize=undefined` and `-fsanitize=address`
builds of the package.

## R CMD check results

0 errors | 0 warnings | 1 note

* The `BugReports` URL intentionally points to GitLab's unified work-items
  page. `R CMD check --as-cran` heuristically suggests appending `/issues`,
  but that is not the tracker URL used by this project or by GitLab projects
  in general.

Local checks on Arch Linux also note non-portable compilation flags. They
come from that R installation's own build configuration. The package sets
no compiler flags.

## Test environments

* Local: Arch Linux, R 4.6.2 Patched
* Local: gcc UBSan and ASan builds of the package's compiled code
* GitLab CI: Linux (rocker/r-ver:4.6.0, R-release)
* GitLab CI: Linux (rocker/r-devel, R-devel)
* GitHub Actions: macOS-latest (R-release)
* GitHub Actions: windows-latest (R-release)
* GitHub Actions: ubuntu-latest (R-release, R-devel, R-oldrel-1)

## Copyright

src/cube.c is a C port of the cube method from BalancedSampling
2.0.6, released under GPL (>= 2) (the package moved to AGPL-3 only
at 2.1.1, after the version ported). The original author, Wilmer
Prentius, is credited as ctb and cph in Authors@R. All other C code
is original, implemented from the published algorithms cited in the
sources.

## Downstream dependencies

None.
