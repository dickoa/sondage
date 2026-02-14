# sondage

Fast survey sampling algorithms for R. Sampling functions return **design objects** with generics for extracting inclusion probabilities, joint inclusion probabilities, and variance estimation quantities.

## Installation

```r
# From GitLab
remotes::install_gitlab("dickoa/sondage")
```

## Usage

```r
library(sondage)

# Use built-in US state data
data(state)
states <- as.data.frame(state.x77)

# Compute inclusion probabilities from population size
pik <- inclusion_prob(states$Population, n = 10)

# Draw a sample (Conditional Poisson Sampling)
s <- unequal_prob_wor(pik, method = "cps")
states[s$sample, ]

# Joint inclusion probabilities for variance estimation
pikl  <- joint_inclusion_prob(s)
delta <- sampling_cov(s)                # pi_ij - pi_i * pi_j
chk   <- sampling_cov(s, scaled = TRUE) # 1 - pi_i * pi_j / pi_ij

# Equal probability sampling
s <- equal_prob_wor(nrow(states), 10)
states[s$sample, ]

# PPS with minimum replacement (Chromy)
hits <- expected_hits(states$Population, n = 10)
s <- unequal_prob_wr(hits, method = "chromy")
joint_expected_hits(s, nsim = 10000)

# Batch sampling for simulations (returns n x nrep matrix)
sim <- unequal_prob_wor(pik, method = "cps", nrep = 10000)
dim(sim)   # 10 x 10000
```

## Sampling functions

**Equal probability without replacement** (`equal_prob_wor`):

- `equal_prob_wor(N, n, method = "srs")` - Simple random sampling
- `equal_prob_wor(N, n, method = "systematic")` - Systematic sampling
- `equal_prob_wor(N, n, method = "bernoulli")` - Bernoulli sampling (random size)

**Equal probability with replacement** (`equal_prob_wr`):

- `equal_prob_wr(N, n, method = "srs")` - Simple random sampling with replacement

**Unequal probability without replacement** (`unequal_prob_wor`):

- `unequal_prob_wor(pik, method = "cps")` - Conditional Poisson / maximum entropy
- `unequal_prob_wor(pik, method = "brewer")` - Brewer's method
- `unequal_prob_wor(pik, method = "systematic")` - Systematic PPS
- `unequal_prob_wor(pik, method = "poisson")` - Poisson sampling (random size)
- `unequal_prob_wor(pik, method = "sps")` - Sequential Poisson sampling (order sampling)
- `unequal_prob_wor(pik, method = "pareto")` - Pareto sampling (order sampling)

**Unequal probability with replacement** (`unequal_prob_wr`):

- `unequal_prob_wr(hits, method = "chromy")` - PPS with minimum replacement
- `unequal_prob_wr(hits, method = "multinomial")` - Multinomial PPS

## Design queries

- `inclusion_prob(x, n)` - Compute inclusion probabilities from size measures
- `inclusion_prob(s)` - Extract inclusion probabilities from a WOR design
- `expected_hits(x, n)` - Compute expected hits from size measures
- `expected_hits(s)` - Extract expected hits from a WR design
- `joint_inclusion_prob(s)` - Joint inclusion probabilities (WOR)
- `joint_expected_hits(s)` - Pairwise expectations E(n_i n_j) (WR)
- `sampling_cov(s)` - Sampling covariance matrix
- `sampling_cov(s, scaled = TRUE)` - Check quantities for SYG variance estimator

## Method comparison

| Method         | Fixed n | Exact pi | All pi_kl > 0         | PRN support |
|----------------|---------|----------|-----------------------|-------------|
| `cps`          | yes     | yes      | yes                   | no          |
| `brewer`       | yes     | yes      | yes                   | no          |
| `systematic`   | yes     | yes      | no                    | no          |
| `poisson`      | no      | yes      | yes                   | yes         |
| `sps`          | yes     | yes      | yes                   | yes         |
| `pareto`       | yes     | yes      | yes                   | yes         |
| `multinomial`  | yes     | --       | yes (with replacement)| no          |
| `chromy`       | yes     | yes      | yes                   | no          |

## References

Brewer, K.R.W. and Donadio, M.E. (2003). The High Entropy Variance of the Horvitz-Thompson Estimator. *Survey Methodology*, 29(2), 189-196.

Chromy, J.R. (2009). Some generalizations of the Horvitz-Thompson estimator. *Memorial JSM*.

Tille, Y. (2006). *Sampling Algorithms*. Springer.
