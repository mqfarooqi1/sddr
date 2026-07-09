# Full-rank Markov chain

A Markov chain over *ranks* rather than discretised classes: each unit's
rank within the cross-section (1 = highest value, `n` = lowest) becomes
its state, and the `n` by `n` matrix records how ranks transition from
one period to the next (Rey, 2014). Unlike
[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md) no
binning is needed — the full rank ordering is retained.

## Usage

``` r
full_rank_markov(data, id, time, value)
```

## Arguments

- data:

  A data frame in long format with one row per unit and period. The
  panel must be balanced (every unit present in every period).

- id, time, value:

  Character scalars naming the unit, period, and numeric value columns.

## Value

An object of class `sddr_markov` (an `n` by `n` rank transition matrix,
its counts, and ergodic distribution).

## References

Rey, S. J. (2014). Fast algorithms for a space-time concordance measure.
*Computational Statistics*, 29(3-4), 799-811.

## See also

[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md),
[`geo_rank_markov()`](https://mqfarooqi1.github.io/sddr/reference/geo_rank_markov.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  id = rep(1:6, times = 8),
  time = rep(1:8, each = 6),
  value = rnorm(48)
)
full_rank_markov(df, "id", "time", "value")
#> <sddr> full-rank Markov chain
#>   units: 6 | transitions: 42 | classes: 6 | breaks: fixed
#> 
#> Transition probability matrix (rows = from, cols = to):
#>       1     2     3     4     5     6
#> 1 0.286 0.000 0.143 0.143 0.143 0.286
#> 2 0.429 0.000 0.000 0.000 0.143 0.429
#> 3 0.143 0.286 0.143 0.143 0.143 0.143
#> 4 0.000 0.143 0.429 0.286 0.143 0.000
#> 5 0.143 0.286 0.000 0.000 0.429 0.143
#> 6 0.000 0.286 0.286 0.429 0.000 0.000
#> 
#> Ergodic (steady-state) distribution:
#>     1     2     3     4     5     6 
#> 0.167 0.167 0.167 0.167 0.167 0.167 
```
