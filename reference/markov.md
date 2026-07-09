# Classic Markov chain for distribution dynamics

Discretises a continuous longitudinal variable into ordered classes and
estimates a discrete-time, first-order Markov transition probability
matrix together with its ergodic (steady-state) distribution. This is
the classic distribution-dynamics estimator of Quah (1993): the
transition matrix summarises how units move between parts of the
distribution from one period to the next, and the steady-state
distribution describes where the process settles in the long run.

Input is *tidy*, long-format panel data with one row per unit-period.
The estimator is agnostic to what the units are (regions, firms,
households, ...) and to whether space is involved; spatial conditioning
is added by a companion
[`spatial_markov()`](https://mqfarooqi1.github.io/sddr/reference/spatial_markov.md)
function (available in a later release).

## Usage

``` r
markov(data, id, time, value, k = 5, breaks = NULL, fixed = TRUE)
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

- k:

  Integer number of classes to discretise `value` into. Default `5`
  (quintiles). Ignored if `breaks` is supplied.

- breaks:

  Optional numeric vector of *interior* class boundaries. If supplied
  these fixed cut points are used instead of data-driven quantiles,
  giving `length(breaks) + 1` classes.

- fixed:

  Logical. If `TRUE` (default) a single set of quantile breaks is
  estimated from the pooled distribution and applied to every period
  (absolute mobility). If `FALSE`, breaks are re-estimated within each
  period (relative mobility). Ignored when `breaks` is supplied.

## Value

An object of class `sddr_markov`: a list with elements

- `matrix`:

  The `k` by `k` row-stochastic transition probability matrix `P`.

- `transitions`:

  The `k` by `k` matrix of observed transition counts.

- `steady_state`:

  The ergodic (stationary) distribution of `P`.

- `classes`:

  Number of classes `k`.

- `breaks`:

  The class boundaries used (when `fixed`).

- `n`:

  Number of units contributing at least one transition.

- `n_transitions`:

  Total number of observed transitions.

- `data`:

  The input data augmented with the integer `state` column.

## References

Quah, D. (1993). Empirical cross-section dynamics in economic growth.
*European Economic Review*, 37(2-3), 426-434.

## See also

[`steady_state()`](https://mqfarooqi1.github.io/sddr/reference/steady_state.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  id = rep(1:50, each = 4),
  time = rep(2000:2003, times = 50),
  value = rnorm(200)
)
m <- markov(df, id = "id", time = "time", value = "value", k = 4)
m
#> <sddr> classic Markov chain
#>   units: 50 | transitions: 150 | classes: 4 | breaks: fixed
#> 
#> Transition probability matrix (rows = from, cols = to):
#>       1     2     3     4
#> 1 0.225 0.225 0.275 0.275
#> 2 0.229 0.257 0.229 0.286
#> 3 0.306 0.333 0.111 0.250
#> 4 0.205 0.154 0.333 0.308
#> 
#> Ergodic (steady-state) distribution:
#>     1     2     3     4 
#> 0.240 0.239 0.241 0.281 
m$steady_state
#>         1         2         3         4 
#> 0.2396739 0.2387903 0.2408205 0.2807152 
```
