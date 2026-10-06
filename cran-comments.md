## Purpose of this submission

This release adds one feature and fixes two defects present in 0.9.1.

* `balanced_wor(method = "scps")` failed with "SCPS maximal weights are
  numerically infeasible" whenever `sum(pik)` missed an integer by a residue
  the input check accepts. `inclusion_prob()` output routinely carries such
  residues, so draws from realistic strata failed, in some cases on every
  draw. The tolerance now scales with the residue.
* `unequal_prob_wor(method = "cps")` did not converge when drawing one unit
  (or all but one) with two inclusion probabilities near 0.5, warned, and
  drew off target by up to 7e-3. One draw now uses its exact closed form.
* New: `balanced_wor()` gains `prn` for coordinating samples with permanent
  random numbers, supported by `method = "scps"`, and registered balanced
  methods may declare PRN support. `scps` is also faster when many units
  share a spreading value.

The compiled code changed in `src/scps.c`, `src/cps_core.h` and
`src/init.c`. The full test suite and a stress run of the new code paths are
clean under gcc `-fsanitize=undefined` and `-fsanitize=address` builds of
the package.

## R CMD check results

0 errors | 0 warnings | 1 note

* The `BugReports` URL intentionally points to GitLab's unified work-items
  page. `R CMD check --as-cran` heuristically suggests appending `/issues`,
  but that is not the tracker URL used by this project or by GitLab projects
  in general.

## Test environments

* Local: Arch Linux, R 4.6.1 Patched
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
