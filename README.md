# sondage

Fast survey sampling algorithms for R. All functions return **indices** for easy subsetting.

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

# Simple random sample: 10 states
idx <- srs(10, nrow(states))
states[idx, ]

# Unequal probability sampling: probability proportional to population
pik <- inclusion_prob(states$Population, n = 10)
idx <- up_maxent(pik)
states[idx, ]

# Batch sampling for simulations (returns n × nrep matrix)
sim_cps <- up_maxent(pik, nrep = 10000)
```

## Functions

**Equal probability:**

- `srs()` - Simple random sampling
- `systematic()` - Systematic sampling
- `bernoulli()` - Bernoulli sampling

**Unequal probability:**

- `up_maxent()` - Maximum entropy / Conditional Poisson
- `up_brewer()` - Brewer's method
- `up_systematic()` - Systematic PPS
- `up_poisson()` - Poisson sampling
- `up_multinomial()` - PPS with replacement

**Joint inclusion probabilities :**

- `up_maxent_joint()` - Exact CPS joint probabilities (Aires' formula)
- `up_brewer_joint()` - Brewer's approximation (equation 18)
- `up_systematic_joint()` - Exact systematic joint probabilities
- `up_poisson_joint()` - Independent selections (π_i × π_j)

**Utilities:**

- `inclusion_prob()` - Compute π from measure of size

## Method comparison

| Method           | Fixed n | Exact π | All π_kl > 0         |
|------------------|---------|---------|----------------------|
| `up_maxent`      | ✓       | ✓       | ✓                    |
| `up_brewer`      | ✓       | ✓       | ✓                    |
| `up_systematic`  | ✓       | ✓       | ✗                    |
| `up_poisson`     | ✗       | ✓       | ✓                    |
| `up_multinomial` | ✓       | —       | ✓ (with replacement) |

## References

Brewer, K.R.W. and Donadio, M.E. (2003). The High Entropy Variance of the Horvitz-Thompson Estimator. *Survey Methodology*, 29(2), 189-196.

Tillé, Y. (2006). *Sampling Algorithms*. Springer.
