# sondage 0.7.0

Initial CRAN release.

## Sampling dispatchers

Five dispatchers returning S3 design objects with class
`c("<prob_class>", "<replacement>", "sondage_sample")`:

* `equal_prob_wor(N, n, method=)` - Simple random sampling, systematic, Bernoulli (random size).
* `equal_prob_wr(N, n, method=)` - Simple random sampling with replacement.
* `unequal_prob_wor(pik, method=)` - Conditional Poisson / maximum entropy (CPS), Brewer, systematic PPS, Poisson (random size), Sequential Poisson (SPS), Pareto.
* `unequal_prob_wr(hits, method=)` - Chromy's minimum replacement, multinomial PPS.
* `balanced_wor(pik, aux, strata, method=)` - Cube method for balanced sampling (Deville & Tille, 2004), with optional stratification (Chauvet & Tille, 2006; Chauvet, 2009).

## Design queries

Five generics for variance estimation quantities:

* `inclusion_prob()` - Compute from size measures or extract from a WOR design.
* `expected_hits()` - Compute from size measures or extract from a WR design.
* `joint_inclusion_prob()` - N x N matrix of joint inclusion probabilities (WOR). Exact for CPS, systematic, Poisson, SRS, Bernoulli; high-entropy approximation (Brewer & Donadio, 2003) for Brewer, SPS, Pareto, cube. `sampled_only = TRUE` returns only the n x n submatrix for sampled units.
* `joint_expected_hits()` - N x N matrix of pairwise expectations E(n_i n_j) (WR). Exact analytic for multinomial and SRS and simulation-based for Chromy (`nsim` parameter). `sampled_only = TRUE` returns only the submatrix for selected units.
* `sampling_cov()` - Sampling covariance matrix. `weighted = TRUE` returns Sen-Yates-Grundy check quantities. `sampled_only = TRUE` computes only the submatrix for sampled units.

## Sampling methods (13 total)

| Method      | Dispatcher       | Fixed size | PRN | Joint probs         |
|-------------|------------------|------------|-----|---------------------|
| srs         | equal_prob_wor   | yes        | no  | Exact               |
| systematic  | equal_prob_wor   | yes        | no  | Exact (some = 0)    |
| bernoulli   | equal_prob_wor   | no         | yes | Exact (independent) |
| srs         | equal_prob_wr    | yes        | no  | Exact               |
| cps         | unequal_prob_wor | yes        | no  | Exact               |
| brewer      | unequal_prob_wor | yes        | no  | HE approximation    |
| systematic  | unequal_prob_wor | yes        | no  | Exact (some = 0)    |
| poisson     | unequal_prob_wor | no         | yes | Exact (independent) |
| sps         | unequal_prob_wor | yes        | yes | HE approximation    |
| pareto      | unequal_prob_wor | yes        | yes | HE approximation    |
| chromy      | unequal_prob_wr  | yes        | no  | Simulated           |
| multinomial | unequal_prob_wr  | yes        | no  | Exact               |
| cube        | balanced_wor     | yes        | no  | HE approximation    |

## Features

* All functions return design objects usable with `s$sample`, `s$pik`, `s$hits`, etc.
* Batch sampling via `nrep` argument for Monte Carlo simulations. Fixed-size designs return a matrix; random-size designs return a list.
* Permanent random numbers (`prn`) for sample coordination (Bernoulli, Poisson, SPS, Pareto).
* C implementations for CPS calibration and sampling, Brewer's draw-by-draw, Chromy's sequential PPS, systematic PPS joint probabilities, high-entropy joint probabilities, cube flight/landing phases, and inclusion probability capping.
* CPS batch optimisation: the O(N^2) design matrix is computed once, then each replicate is a cheap sequential draw.
* Cube batch optimisation: C-level batch entry points reuse the workspace allocation and the A = X/pi matrix across replicates, with precomputed X/pi helpers for stratified designs.
* Cube auxiliary conditioning: `condition_aux = TRUE` pre-conditions `aux` by weighted centering/scaling and QR-pivot rank pruning, improving numerical stability with ill-conditioned or collinear auxiliary variables.
* Stratified cube: within-stratum size constraints are preserved exactly when per-stratum `sum(pik)` is close to an integer; otherwise a warning is issued and the design is marked random-size.
* High-entropy joint inclusion probabilities handle certainty units (pi_k = 1) correctly and clamp values to valid bounds. A warning is issued when the marginal defect exceeds 5% of n.
* N-size guard: `joint_inclusion_prob()` and `joint_expected_hits()` refuse N > 10,000 to prevent accidental allocation of multi-GB dense matrices. `sampled_only = TRUE` bypasses this limit for methods that can compute the n x n submatrix directly (HE, Poisson, Bernoulli, SRS, multinomial); for CPS, systematic, and Chromy the guard still applies because the algorithms require an N x N intermediate.
* `sampling_cov(weighted = TRUE)` returns `NA` (not `NaN`) for undefined entries (zero joint probabilities or zero selection probabilities), with an informative warning.
* All long-running C loops call `R_CheckUserInterrupt()` so Ctrl-C works.
* CPS calibration warns on non-convergence.
* Brewer's denominator is clamped (not zeroed) for numerical safety with extreme inclusion probabilities.
