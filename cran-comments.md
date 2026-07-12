## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission. The `BugReports` URL intentionally points to
  GitLab's unified work-items page. `R CMD check --as-cran` heuristically
  suggests appending `/issues`, but that is not the tracker URL used by this
  project.
* The local Arch Linux toolchain injects non-portable compilation flags
  (`-march=x86-64`, hardening, and format-security flags). The package has
  no `Makevars` file and does not set compiler flags.

## Test environments

* Local: Arch Linux, R 4.6.1 Patched

## Copyright

src/cube.c is a C port of the cube method from BalancedSampling
2.0.6, released under GPL (>= 2) (the package moved to AGPL-3 only
at 2.1.1, after the version ported). The original author, Wilmer
Prentius, is credited as ctb and cph in Authors@R. All other C code
is original, implemented from the published algorithms cited in the
sources.

## Downstream dependencies

None (new package).
