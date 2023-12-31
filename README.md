
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sondage

<!-- badges: start -->
<!-- badges: end -->

The goal of sondage is to …

## Installation

You can install the development version of sondage like so:

``` r
# install.packages("sondage")
pak::pkg_install("dickoa/sondage")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sondage)

set.seed(123)
p <- c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9) # pik
N <- length(p)
marginal <- vector(mode = "numeric", length = N)
sequential <- vector(mode = "numeric", length = N)
nsim <- 1e5
for(i in 1:nsim){
  s1 <- sample_int(N, prob = p,
                   replace = FALSE,
                   method = "sequential")
  s2 <- sample_int(N, prob = p,
                   replace = FALSE,
                   method = "marginal")
  sequential[s1] <- sequential[s1] + 1
  marginal[s2] <- marginal[s2] + 1
}

cbind(pik = p,
      sequential = sequential/nsim, marginal = marginal/nsim)
```
