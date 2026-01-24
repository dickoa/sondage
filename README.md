# sondage

Survey Sampling Algorithms for R

## Overview

Fast implementations of survey sampling algorithms for drawing samples from finite populations. All functions return **indices** for easy subsetting with `df[idx, ]` or `dplyr::slice()`.

### Equal Probability Methods
- `srs()` - Simple random sampling (with/without replacement)
- `systematic()` - Systematic sampling
- `bernouilli()` - Bernoulli sampling (random size)

### Unequal Probability Methods
- `up_maxent()` - Maximum entropy / Conditional Poisson Sampling
- `up_brewer()` - Brewer's method (Algorithm 6.10 from Tillé)
- `up_systematic()` - Systematic PPS
- `up_poisson()` - Poisson sampling (random size)
- `up_multinomial()` - PPS with replacement

### Utilities
- `inclusion_prob()` - Compute inclusion probabilities from size measure


## Why sondage?

The [sampling](https://CRAN.R-project.org/package=sampling) package is
excellent and comprehensive. sondage focuses on a narrower goal:

|                | sondage                     | sampling            |
|----------------|-----------------------------|---------------------|
| **Focus**      | Core drawing algorithms     | Full survey toolkit |
| **Output**     | Indices (`df[idx, ]`)       | 0/1 indicators      |
| **Batch mode** | `up_maxent(pik, nrep=1000)` | Loop required       |
| **Speed**      | Optimized C                 | Mixed               |

**Use sondage when:**
- You need fast repeated sampling (simulations, variance estimation)
- You want a simple, consistent API
- You prefer indices over indicators

**Use sampling when:**
- You need balanced sampling (cube method)
- You need calibration, stratification utilities
- You want a battle-tested CRAN package

## Installation

```r
# Install from source
remotes::install_gitlab("dickoa/sondage")
```

## Quick Start

```r
library(sondage)

# Simple random sample
idx <- srs(50, 1000)
sample_df <- df[idx, ]

# Unequal probability sampling
pik <- inclusion_prob(df$revenue, n = 50)
idx <- up_maxent(pik)
sample_df <- df[idx, ]

# With replacement (indices may repeat)
idx <- srs(100, 1000, replace = TRUE)
bootstrap_df <- df[idx, ]

# Simulation with batch mode (much faster than loop)
samples <- up_maxent(pik, nrep = 10000)  # n × 10000 matrix
```

## Method Comparison

| Property              | up_maxent | up_brewer | up_systematic |
|-----------------------|-----------|-----------|---------------|
| Fixed sample size     | ✓         | ✓         | ✓             |
| Exact inclusion probs | ✓         | ✓         | ✓             |
| All joint probs > 0   | ✓         | ✓         | ✗             |
| Order invariant       | ✓         | ✓         | ✗             |
| Maximum entropy       | ✓         | ✗         | ✗             |
| Speed                 | Fast      | Fast      | Very fast     |

## Why Maximum Entropy?

`up_maxent()` implements Conditional Poisson Sampling, the unique design that maximizes entropy subject to fixed inclusion probabilities:

1. **Optimal variance** - Minimizes sampling variance
2. **Tractable theory** - Joint inclusion probabilities have closed form
3. **Fast batch mode** - Design computed once, reused for all replicates

## References

- Tillé, Y. (2006). *Sampling Algorithms*. Springer.
- Chen, S.X., Dempster, A.P., and Liu, J.S. (1994). Weighted finite population
  sampling to maximize entropy. *Biometrika*, 81, 457-469.
- Brewer, K.R.W. (1963). A model of systematic sampling with unequal
  probabilities. *Australian Journal of Statistics*, 5, 5-13.

## License

MIT
