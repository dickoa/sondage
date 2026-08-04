## Purpose of this submission

This release fixes the undefined behavior reported by the CRAN gcc-UBSan
check on 0.9.0. `src/joint_probs.c` allocated its `long double` scratch
buffers with `R_alloc()`, which guarantees only the alignment required by
`double`; on x86-64 `long double` requires 16-byte alignment, so the loads
and stores in the Sampford joint-probability code were misaligned. They now
use `R_allocLD()`. Computed values are unchanged.

The short interval since 0.9.0 is because this addresses that report.

I reproduced all 33 reported locations locally with a `-fsanitize=undefined`
build, confirmed the fix silences every one of them, and confirmed the full
test suite is also clean under `-fsanitize=address`.

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
