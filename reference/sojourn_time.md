# Sojourn times

The expected number of consecutive periods a unit remains in a class
before leaving it, \\1 / (1 - p\_{ii})\\ for each class \\i\\. Large
values indicate "sticky" classes with strong persistence.

## Usage

``` r
sojourn_time(x, ...)

# Default S3 method
sojourn_time(x, ...)

# S3 method for class 'sddr_markov'
sojourn_time(x, ...)
```

## Arguments

- x:

  A row-stochastic matrix, or an object of class `sddr_markov`.

- ...:

  Currently unused.

## Value

A named numeric vector of sojourn times (one per class). Absorbing
classes (\\p\_{ii} = 1\\) return `Inf`.

## See also

[`mfpt()`](https://mqfarooqi1.github.io/sddr/reference/mfpt.md),
[`steady_state()`](https://mqfarooqi1.github.io/sddr/reference/steady_state.md)

## Examples

``` r
P <- matrix(c(0.8, 0.2,
              0.3, 0.7), nrow = 2, byrow = TRUE)
sojourn_time(P)
#>        1        2 
#> 5.000000 3.333333 
```
