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

# Simple random sample
idx <- srs(50, 1000)
sample_df <- df[idx, ]

# Unequal probability sampling
pik <- inclusion_prob(df$size, n = 50)
idx <- up_maxent(pik)
sample_df <- df[idx, ]

# Balanced sampling (cube method)
pik <- inclusion_prob(df$size, n = 50)
x <- cbind(pik, df$region, df$income)  # balancing variables
idx <- cube(pik, x)
sample_df <- df[idx, ]

# Batch sampling for simulations
samples <- up_maxent(pik, nrep = 10000)  # n × 10000 matrix
samples <- cube(pik, x, nrep = 10000)    # n × 10000 matrix
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

**Balanced sampling:**

- `cube()` - Cube method (Deville & Tillé)

**Utilities:**

- `inclusion_prob()` - Compute π from measure of size

## Method comparison

| Method          | Fixed n | Exact π | All π_kl > 0 | Balanced |
|-----------------|---------|---------|--------------|----------|
| `up_maxent`     | ✓       | ✓       | ✓            | ✗        |
| `up_brewer`     | ✓       | ✓       | ✓            | ✗        |
| `up_systematic` | ✓       | ✓       | ✗            | ✗        |
| `cube`          | ✓*      | ✓       | ✓            | ✓        |

\* Fixed n when inclusion probabilities are included as a balancing variable.

## References

Deville, J.C. and Tillé, Y. (2004). Efficient balanced sampling: the cube method. *Biometrika*, 91(4), 893-912.

Tillé, Y. (2006). *Sampling Algorithms*. Springer.
