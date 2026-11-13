## Purpose of this submission

This release follows 0.9.1 on CRAN. It fixes defects present in 0.9.1, the
most serious being wrong joint inclusion probabilities for Sampford designs
at moderate sizes and for CPS designs with nearly equal probabilities. It
also adds three features and speeds up several samplers. `NEWS.md` lists
every change.

The full test suite and a stress run of the new code paths are clean under
gcc `-fsanitize=undefined` and `-fsanitize=address` builds of the package.

## R CMD check results

0 errors | 0 warnings | 1 note

* The `BugReports` URL intentionally points to GitLab's work-items page.
  `R CMD check --as-cran` suggests appending `/issues`, but that is not the
  tracker URL GitLab uses.

Local checks on Arch Linux also note non-portable compilation flags. They
come from that R installation's build configuration. The package sets no
compiler flags.

## Test environments

* Local: Arch Linux, R 4.6.2 Patched
* Local: gcc UBSan and ASan builds of the package's compiled code
* GitLab CI: Linux (rocker/r-ver:4.6.0, R-release)
* GitLab CI: Linux (rocker/r-devel, R-devel)
* GitHub Actions: macOS-latest (R-release)
* GitHub Actions: windows-latest (R-release)
* GitHub Actions: ubuntu-latest (R-release, R-devel, R-oldrel-1)

## Copyright

`src/cube.c` is a C port of the cube method from BalancedSampling 2.0.6,
released under GPL (>= 2). The package moved to AGPL-3 only at 2.1.1, after
the ported version. Its author, Wilmer Prentius, is credited as ctb and cph
in Authors@R. All other C code is original.

## Downstream dependencies

None.
