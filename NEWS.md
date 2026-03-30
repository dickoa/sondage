# sondage 0.8

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

## Sampling methods (13 built-in)

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

## Custom method registration

* `register_method()` lets users plug custom unequal-probability sampling
  algorithms into the existing dispatchers and generics. Registered methods
  work with `unequal_prob_wor()`, `unequal_prob_wr()`, `joint_inclusion_prob()`,
  `sampling_cov()`, batch mode (`nrep`), and all downstream tooling.
* `method_spec()` queries method metadata (type, fixed size, PRN support)
  for any method, built-in or registered. Covers all five dispatchers
  including balanced sampling.
* Helper functions `registered_methods()`, `is_registered_method()`, and
  `unregister_method()` manage the registry.
* `he_jip()` and `hajek_jip()` are exported joint inclusion probability
  approximations that can be passed directly as `joint_fn` to
  `register_method()`. `he_jip()` uses the Brewer & Donadio (2003)
  high-entropy formula (C implementation); `hajek_jip()` uses the
  Hajek (1964) conditional Poisson approximation.
* New vignette "Extending sondage with Custom Methods" with worked examples
  (Sampford's method with `he_jip()`, wrapping `sampling::UPtille`).

## Features

* All functions return design objects usable with `s$sample`, `s$pik`, `s$hits`, etc.
* Batch sampling via `nrep` argument for Monte Carlo simulations. Fixed-size designs return a matrix; random-size designs return a list.
* Permanent random numbers (`prn`) for sample coordination (Bernoulli, Poisson, SPS, Pareto).
* C implementations for all sampling algorithms.
* CPS batch optimisation: the O(N^2) design matrix is computed once, then each replicate is a cheap sequential draw.
* Cube batch optimisation: C-level batch entry points reuse the workspace allocation across replicates.
* Cube auxiliary conditioning: `condition_aux = TRUE` pre-conditions `aux` by weighted centering/scaling and QR-pivot rank pruning.
* Stratified cube: within-stratum size constraints are preserved exactly when per-stratum `sum(pik)` is close to an integer.
* High-entropy joint inclusion probabilities handle certainty units correctly and clamp values to valid bounds.
* N-size guard: `joint_inclusion_prob()` and `joint_expected_hits()` refuse N > 10,000 to prevent accidental allocation of multi-GB dense matrices. `sampled_only = TRUE` bypasses this limit.
* All long-running C loops call `R_CheckUserInterrupt()` so Ctrl-C works.
