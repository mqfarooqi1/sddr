# Ergodic (steady-state) distribution of a Markov chain

Computes the long-run stationary distribution of a row-stochastic
transition matrix, i.e. the probability vector \\\pi\\ satisfying \\\pi
P = \pi\\ with \\\sum_i \pi_i = 1\\. It is obtained as the left
eigenvector of `P` associated with the eigenvalue closest to 1.

## Usage

``` r
steady_state(x, ...)

# Default S3 method
steady_state(x, ...)

# S3 method for class 'sddr_markov'
steady_state(x, ...)
```

## Arguments

- x:

  A row-stochastic matrix, or an object of class `sddr_markov`.

- ...:

  Currently unused.

## Value

A named numeric vector giving the stationary distribution.

## Examples

``` r
P <- matrix(c(0.8, 0.2,
              0.3, 0.7), nrow = 2, byrow = TRUE)
steady_state(P)
#>   1   2 
#> 0.6 0.4 
```
