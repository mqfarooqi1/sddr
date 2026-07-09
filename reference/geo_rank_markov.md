# Geographic-rank Markov chain

A Markov chain that tracks, for each rank position, *which unit*
occupies it over time (Rey, 2016). States are the geographic units, so
the `n` by `n` matrix records how units exchange rank positions — a
spatially explicit view of positional mobility.

## Usage

``` r
geo_rank_markov(data, id, time, value)
```

## Arguments

- data:

  A data frame in long format with one row per unit and period. The
  panel must be balanced (every unit present in every period).

- id, time, value:

  Character scalars naming the unit, period, and numeric value columns.

## Value

An object of class `sddr_markov` (an `n` by `n` transition matrix over
units, its counts, and ergodic distribution).

## References

Rey, S. J. (2016). Space-time patterns of rank concordance: Local
indicators of mobility association. *Annals of the AAG*, 106(4),
788-803.

## See also

[`full_rank_markov()`](https://mqfarooqi1.github.io/sddr/reference/full_rank_markov.md),
[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  id = rep(1:6, times = 8),
  time = rep(1:8, each = 6),
  value = rnorm(48)
)
geo_rank_markov(df, "id", "time", "value")
#> <sddr> geo-rank Markov chain
#>   units: 6 | transitions: 42 | classes: 6 | breaks: fixed
#> 
#> Transition probability matrix (rows = from, cols = to):
#>       1     2     3     4     5     6
#> 1 0.286 0.000 0.143 0.286 0.143 0.143
#> 2 0.000 0.143 0.286 0.143 0.143 0.286
#> 3 0.143 0.143 0.143 0.286 0.143 0.143
#> 4 0.143 0.571 0.000 0.000 0.286 0.000
#> 5 0.143 0.143 0.143 0.143 0.286 0.143
#> 6 0.286 0.000 0.286 0.143 0.000 0.286
#> 
#> Ergodic (steady-state) distribution:
#>     1     2     3     4     5     6 
#> 0.167 0.167 0.167 0.167 0.167 0.167 
```
