# LISA Markov chain

A Markov chain over the four quadrants of the Moran scatterplot (Rey,
2001). Each unit, in each period, is placed in a quadrant according to
the sign of its value (relative to the cross-sectional mean) and the
sign of its spatial lag:

- `HH`:

  high value, high neighbourhood (quadrant 1).

- `LH`:

  low value, high neighbourhood (quadrant 2).

- `LL`:

  low value, low neighbourhood (quadrant 3).

- `HL`:

  high value, low neighbourhood (quadrant 4).

The resulting 4 by 4 transition matrix describes how units move between
local spatial-association types over time — capturing the co-evolution
of a unit and its neighbourhood.

## Usage

``` r
lisa_markov(data, id, time, value, weights, row_standardize = TRUE)
```

## Arguments

- data:

  A data frame in long format with one row per unit and period. The
  panel must be balanced (every unit present in every period).

- id, time, value:

  Character scalars naming the unit, period, and numeric value columns.

- weights:

  A square numeric spatial weights matrix whose row and column names
  match the unit ids.

- row_standardize:

  Logical; row-standardise `weights` before computing the spatial lag
  (default `TRUE`).

## Value

An object of class `sddr_markov`: a 4 by 4 transition matrix over the
quadrants `HH`, `LH`, `LL`, `HL`, plus the per-unit-per-period
`quadrants` matrix.

## References

Rey, S. J. (2001). Spatial empirics for economic growth and convergence.
*Geographical Analysis*, 33(3), 195-214.
[doi:10.1111/j.1538-4632.2001.tb00444.x](https://doi.org/10.1111/j.1538-4632.2001.tb00444.x)

## See also

[`spatial_markov()`](https://mqfarooqi1.github.io/sddr/reference/spatial_markov.md),
[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md)

## Examples

``` r
set.seed(1)
n <- 9; periods <- 8
W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
for (i in 1:n) {                       # ring contiguity
  W[i, i %% n + 1] <- 1
  W[i, (i - 2) %% n + 1] <- 1
}
df <- data.frame(id = rep(1:n, times = periods),
                 time = rep(1:periods, each = n),
                 value = rnorm(n * periods))
lisa_markov(df, "id", "time", "value", weights = W)
#> <sddr> LISA Markov chain
#>   units: 9 | transitions: 63 | classes: 4 | breaks: fixed
#> 
#> Transition probability matrix (rows = from, cols = to):
#>       HH    LH    LL    HL
#> HH 0.316 0.368 0.263 0.053
#> LH 0.471 0.118 0.235 0.176
#> LL 0.308 0.154 0.308 0.231
#> HL 0.143 0.357 0.143 0.357
#> 
#> Ergodic (steady-state) distribution:
#>    HH    LH    LL    HL 
#> 0.321 0.251 0.245 0.183 
```
