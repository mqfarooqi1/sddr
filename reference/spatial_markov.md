# Spatial Markov chain for distribution dynamics

Estimates a spatial Markov chain (Rey, 2001): a set of class-transition
probability matrices, each conditioned on the *spatial context* of a
unit at the start of the period. Where
[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md) asks
"given a unit's own class, where does it move next?", `spatial_markov()`
asks "given a unit's own class **and the class of its spatial
neighbourhood**, where does it move next?". Comparing the conditional
matrices with the pooled (aspatial) matrix reveals whether, and how,
spatial context shapes distribution dynamics.

The neighbourhood is summarised by the *spatial lag*: a
(row-standardised) weighted average of neighbours' values. Supply either
a spatial weights matrix via `weights` (the lag is computed for you, per
period), or a precomputed spatial-lag column via `lag` (keeping the
estimator independent of any particular spatial-data stack).

## Usage

``` r
spatial_markov(
  data,
  id,
  time,
  value,
  weights = NULL,
  lag = NULL,
  k = 5,
  m = k,
  breaks = NULL,
  lag_breaks = NULL,
  fixed = TRUE,
  row_standardize = TRUE
)
```

## Arguments

- data:

  A data frame in long format with one row per unit and period.

- id:

  Name (character scalar) of the column identifying the unit.

- time:

  Name (character scalar) of the column giving the period. Periods are
  compared in sorted order within each unit.

- value:

  Name (character scalar) of the numeric column to analyse.

- weights:

  Optional square numeric spatial weights matrix whose row and column
  names match the unit ids. Used to compute the spatial lag of `value`
  within each period. Requires a balanced panel. Exactly one of
  `weights` or `lag` must be supplied.

- lag:

  Optional name (character scalar) of a column giving a precomputed
  spatial lag of `value`.

- k:

  Integer number of classes to discretise `value` into. Default `5`
  (quintiles). Ignored if `breaks` is supplied.

- m:

  Integer number of spatial-lag classes to condition on. Defaults to
  `k`.

- breaks:

  Optional numeric vector of *interior* class boundaries. If supplied
  these fixed cut points are used instead of data-driven quantiles,
  giving `length(breaks) + 1` classes.

- lag_breaks:

  Optional numeric vector of *interior* class boundaries for the spatial
  lag (the lag analogue of `breaks`). If supplied, these fixed cut
  points are used instead of data-driven quantiles. Supplying both
  `breaks` and `lag_breaks` makes the classification fully explicit,
  which is how exact agreement with `giddy` (`cutoffs` / `lag_cutoffs`)
  is obtained.

- fixed:

  Logical. If `TRUE` (default) a single set of quantile breaks is
  estimated from the pooled distribution and applied to every period
  (absolute mobility). If `FALSE`, breaks are re-estimated within each
  period (relative mobility). Ignored when `breaks` is supplied.

- row_standardize:

  Logical; row-standardise `weights` before computing the lag (default
  `TRUE`). Ignored when `lag` is supplied.

## Value

An object of class `sddr_spatial_markov`: a list with elements

- `matrices`:

  Length-`m` list of `k` by `k` conditional transition probability
  matrices, one per spatial-lag class.

- `transitions`:

  Length-`m` list of conditional count matrices.

- `pooled`:

  The `k` by `k` pooled (aspatial) transition matrix; this equals the
  matrix from
  [`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md)
  with the same class definition.

- `pooled_transitions`:

  Pooled count matrix.

- `steady_states`:

  Length-`m` list of conditional ergodic distributions (`NA` where a
  conditional matrix has unvisited classes).

- `classes`, `lag_classes`:

  `k` and `m`.

- `n`, `n_transitions`:

  Units contributing, and total transitions.

- `data`:

  Input augmented with `lagval`, value state and lag state.

## References

Rey, S. J. (2001). Spatial empirics for economic growth and convergence.
*Geographical Analysis*, 33(3), 195-214.
[doi:10.1111/j.1538-4632.2001.tb00444.x](https://doi.org/10.1111/j.1538-4632.2001.tb00444.x)

## See also

[`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md)

## Examples

``` r
set.seed(1)
n <- 12; periods <- 8
# Ring contiguity: each unit neighbours the next and previous on a circle.
W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
for (i in 1:n) {
  W[i, i %% n + 1] <- 1
  W[i, (i - 2) %% n + 1] <- 1
}
df <- do.call(rbind, lapply(seq_len(periods), function(p)
  data.frame(id = 1:n, time = p, value = rnorm(n))))
sm <- spatial_markov(df, "id", "time", "value", weights = W, k = 3)
sm
#> <sddr> spatial Markov chain
#>   units: 12 | transitions: 84 | value classes: 3 | lag classes: 3 | breaks: fixed
#> 
#> Pooled transition matrix (aspatial):
#>       1     2     3
#> 1 0.267 0.367 0.367
#> 2 0.500 0.214 0.286
#> 3 0.231 0.423 0.346
#> 
#> Conditional on spatial-lag class 1 (31 transitions):
#>       1     2     3
#> 1 0.308 0.154 0.538
#> 2 0.700 0.200 0.100
#> 3 0.375 0.250 0.375
#> 
#> Conditional on spatial-lag class 2 (29 transitions):
#>       1     2     3
#> 1 0.333 0.556 0.111
#> 2 0.100 0.300 0.600
#> 3 0.100 0.600 0.300
#> 
#> Conditional on spatial-lag class 3 (24 transitions):
#>       1     2     3
#> 1 0.125 0.500 0.375
#> 2 0.750 0.125 0.125
#> 3 0.250 0.375 0.375
```
