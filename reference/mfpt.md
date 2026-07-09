# Mean first passage times

The mean first passage time from class \\i\\ to class \\j\\ is the
expected number of periods for a unit currently in \\i\\ to first reach
\\j\\. Together with the transition matrix these times summarise how
quickly a distribution mixes. The diagonal holds the *mean recurrence
time* \\1/\pi_j\\ (the expected return time to a class), where \\\pi\\
is the ergodic distribution.

Computed with the Kemeny-Snell fundamental matrix \\Z = (I - P +
A)^{-1}\\, where \\A\\ has every row equal to \\\pi\\; then \\M\_{ij} =
(\delta\_{ij} - Z\_{ij} + Z\_{jj}) / \pi_j\\.

## Usage

``` r
mfpt(x, ...)

# Default S3 method
mfpt(x, ...)

# S3 method for class 'sddr_markov'
mfpt(x, ...)
```

## Arguments

- x:

  A row-stochastic matrix, or an object of class `sddr_markov`.

- ...:

  Currently unused.

## Value

A numeric matrix of mean first passage times.

## References

Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*. Springer.

## See also

[`steady_state()`](https://mqfarooqi1.github.io/sddr/reference/steady_state.md),
[`sojourn_time()`](https://mqfarooqi1.github.io/sddr/reference/sojourn_time.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2,
              0.3, 0.7), nrow = 2, byrow = TRUE)
mfpt(P)
#>          [,1] [,2]
#> [1,] 1.666667  5.0
#> [2,] 3.333333  2.5
```
