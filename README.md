
# sondage

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sondage)](https://CRAN.R-project.org/package=sondage)
[![R-CMD-check](https://gitlab.com/dickoa/sondage/badges/main/pipeline.svg)](https://gitlab.com/dickoa/sondage/-/pipelines)
[![Codecov test
coverage](https://codecov.io/gl/dickoa/sondage/branch/main/graph/badge.svg)](https://app.codecov.io/gl/dickoa/sondage?branch=main)
<!-- badges: end -->

Fast survey sampling algorithms for R. Sampling functions return
**design objects** with generics for extracting inclusion probabilities,
joint inclusion probabilities, and variance estimation quantities.

## Installation

``` r
# From GitLab
pak::pkg_install("gitlab::dickoa/sondage")
```

## Usage

``` r
library(sondage)

# Use built-in US state data
data(state)
states <- as.data.frame(state.x77)

# Compute inclusion probabilities from population size
pik <- inclusion_prob(states$Population, n = 10)

# Draw a sample (Conditional Poisson Sampling)
s <- unequal_prob_wor(pik, method = "cps")
states[s$sample, ]
#>              Population Income Illiteracy Life Exp Murder HS Grad Frost   Area
#> California        21198   5114        1.1    71.71   10.3    62.6    20 156361
#> Georgia            4931   4091        2.0    68.54   13.9    40.6    60  58073
#> Michigan           9111   4751        0.9    70.63   11.1    52.8   125  56817
#> Mississippi        2341   3098        2.4    68.09   12.5    41.0    50  47296
#> Missouri           4767   4254        0.8    70.69    9.3    48.8   108  68995
#> Nebraska           1544   4508        0.6    72.60    2.9    59.3   139  76483
#> New York          18076   4903        1.4    70.55   10.9    52.7    82  47831
#> Pennsylvania      11860   4449        1.0    70.43    6.1    50.2   126  44966
#> Washington         3559   4864        0.6    71.72    4.3    63.5    32  66570
#> Wisconsin          4589   4468        0.7    72.48    3.0    54.5   149  54464
```

``` r
# Joint inclusion probabilities for variance estimation
pikl  <- joint_inclusion_prob(s)
delta <- sampling_cov(s)                # pi_ij - pi_i * pi_j
chk   <- sampling_cov(s, weighted = TRUE) # 1 - pi_i * pi_j / pi_ij
```

``` r
# Equal probability sampling
s <- equal_prob_wor(nrow(states), 10)
states[s$sample, ]
#>                Population Income Illiteracy Life Exp Murder HS Grad Frost
#> Minnesota            3921   4675        0.6    72.96    2.3    57.6   160
#> Colorado             2541   4884        0.7    72.06    6.8    63.9   166
#> South Carolina       2816   3635        2.3    67.96   11.6    37.8    65
#> Utah                 1203   4022        0.6    72.90    4.5    67.3   137
#> Missouri             4767   4254        0.8    70.69    9.3    48.8   108
#> Wisconsin            4589   4468        0.7    72.48    3.0    54.5   149
#> Rhode Island          931   4558        1.3    71.90    2.4    46.4   127
#> Tennessee            4173   3821        1.7    70.11   11.0    41.8    70
#> Vermont               472   3907        0.6    71.64    5.5    57.1   168
#> Mississippi          2341   3098        2.4    68.09   12.5    41.0    50
#>                  Area
#> Minnesota       79289
#> Colorado       103766
#> South Carolina  30225
#> Utah            82096
#> Missouri        68995
#> Wisconsin       54464
#> Rhode Island     1049
#> Tennessee       41328
#> Vermont          9267
#> Mississippi     47296
```

``` r
# PPS with minimum replacement (Chromy)
hits <- expected_hits(states$Population, n = 10)
s <- unequal_prob_wr(hits, method = "chromy")
```

``` r
# Balanced sampling (cube method)
pik <- inclusion_prob(states$Population, n = 10)
x <- matrix(states$Income)
s_bal <- balanced_wor(pik, aux = x)
s_bal
#> Unequal prob WOR [cube] (n=10, N=50): 5 10 13 14 18 24 25 32 38 44
```

``` r
# Batch sampling for simulations (design object with matrix $sample)
sim <- unequal_prob_wor(pik, method = "cps", nrep = 1000)
dim(sim$sample)   # 10 x 1000
#> [1]   10 1000
inclusion_prob(sim) # generics still work
#>  [1] 0.17026107 0.01719095 0.10418188 0.09937783 0.99839394 0.11967728
#>  [7] 0.14600534 0.02727003 0.38983426 0.23224269 0.04088150 0.03829108
#> [13] 0.52736187 0.25023432 0.13474880 0.10738457 0.15952261 0.17925688
#> [19] 0.04983021 0.19414000 0.27383066 0.42911441 0.18467321 0.11025758
#> [25] 0.22451854 0.03513548 0.07272008 0.02778811 0.03824398 0.34537328
#> [31] 0.05388068 0.85135243 0.25626292 0.03000174 0.50560237 0.12787242
#> [37] 0.10757297 0.55858818 0.04384870 0.13262937 0.03207408 0.19654203
#> [43] 0.57634431 0.05665949 0.02223049 0.23459761 0.16762355 0.08473020
#> [49] 0.21613500 0.01770903
```

## Sampling functions

**Equal probability without replacement** (`equal_prob_wor`):

- `equal_prob_wor(N, n, method = "srs")` - Simple random sampling
- `equal_prob_wor(N, n, method = "systematic")` - Systematic sampling
- `equal_prob_wor(N, n, method = "bernoulli")` - Bernoulli sampling
  (random size)

**Equal probability with replacement** (`equal_prob_wr`):

- `equal_prob_wr(N, n, method = "srs")` - Simple random sampling with
  replacement

**Unequal probability without replacement** (`unequal_prob_wor`):

- `unequal_prob_wor(pik, method = "cps")` - Conditional Poisson /
  maximum entropy
- `unequal_prob_wor(pik, method = "brewer")` - Brewer’s method
- `unequal_prob_wor(pik, method = "systematic")` - Systematic PPS
- `unequal_prob_wor(pik, method = "poisson")` - Poisson sampling (random
  size)
- `unequal_prob_wor(pik, method = "sps")` - Sequential Poisson sampling
  (order sampling)
- `unequal_prob_wor(pik, method = "pareto")` - Pareto sampling (order
  sampling)

**Unequal probability with replacement** (`unequal_prob_wr`):

- `unequal_prob_wr(hits, method = "chromy")` - PPS with minimum
  replacement
- `unequal_prob_wr(hits, method = "multinomial")` - Multinomial PPS

**Balanced sampling without replacement** (`balanced_wor`):

- `balanced_wor(pik, aux, method = "cube")` - Cube method (Deville &
  Tille, 2004)
- `balanced_wor(pik, aux, strata, method = "cube")` - Stratified cube
  (Chauvet, 2009)

## Design queries

- `inclusion_prob(x, n)` - Compute inclusion probabilities from size
  measures
- `inclusion_prob(s)` - Extract inclusion probabilities from a WOR
  design
- `expected_hits(x, n)` - Compute expected hits from size measures
- `expected_hits(s)` - Extract expected hits from a WR design
- `joint_inclusion_prob(s)` - Joint inclusion probabilities (WOR)
- `joint_inclusion_prob(s, sampled_only = TRUE)` - n x n submatrix for
  sampled units only (scales to large N)
- `joint_expected_hits(s)` - Pairwise expectations E(n_i n_j) (WR)
- `joint_expected_hits(s, sampled_only = TRUE)` - Submatrix for selected
  units only
- `sampling_cov(s)` - Sampling covariance matrix
- `sampling_cov(s, weighted = TRUE)` - Check quantities for SYG variance
  estimator
- `sampling_cov(s, sampled_only = TRUE)` - Covariance for sampled units
  only

## Method comparison

| Method        | Fixed n | Exact pi_i | Exact pi_ij       | PRN support |
|---------------|---------|------------|-------------------|-------------|
| `cps`         | yes     | yes        | yes               | no          |
| `brewer`      | yes     | yes        | approx (HE)       | no          |
| `systematic`  | yes     | yes        | yes (some = 0)    | no          |
| `poisson`     | no      | yes        | yes (independent) | yes         |
| `sps`         | yes     | approx\*   | approx (HE)       | yes         |
| `pareto`      | yes     | approx\*   | approx (HE)       | yes         |
| `multinomial` | yes     | yes        | yes (analytic)    | no          |
| `chromy`      | yes     | yes        | simulated         | no          |
| `cube`        | yes     | yes        | approx (HE)       | no          |

\*approx = target probabilities, exact asymptotically. HE = high-entropy
approximation.

## References

Brewer, K.R.W. and Donadio, M.E. (2003). The High Entropy Variance of
the Horvitz-Thompson Estimator. *Survey Methodology*, 29(2), 189-196.

Chauvet, G. (2009). Stratified balanced sampling. *Survey Methodology*,
35, 115-119.

Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson
estimator. *Memorial JSM*.

Deville, J.C. and Tille, Y. (2004). Efficient balanced sampling: the
cube method. *Biometrika*, 91(4), 893-912.

Tille, Y. (2006). *Sampling Algorithms*. Springer.
