# Theta rank-mobility statistic

The Theta statistic (Rey, 2004) measures how much of the overall rank
movement in a system is *within-regime* exchange. For each
period-to-period transition it is the sum over regimes of the absolute
*net* within-regime rank change, divided by the total absolute rank
change. Values near zero indicate that rank movement is dominated by
within-regime shuffling that nets out; larger values indicate coherent
regime-level movement.

## Usage

``` r
theta(y, regime)
```

## Arguments

- y:

  A numeric matrix of values, units in rows and time periods in columns
  (at least two columns).

- regime:

  A vector of length `nrow(y)` assigning each unit to a regime.

## Value

An object of class `sddr_theta`: a list with `theta` (one value per
period-to-period transition), the per-period `total` absolute rank
change, and the column rank matrix `ranks`.

## References

Rey, S. J. (2004). Spatial dependence in the evolution of regional
income distributions. In *Spatial Econometrics and Spatial Statistics*,
194-213.

## See also

[`tau()`](https://mqfarooqi1.github.io/sddr/reference/tau.md)

## Examples

``` r
set.seed(1)
y <- matrix(rnorm(40), nrow = 10)    # 10 units, 4 periods
regime <- rep(c("a", "b"), each = 5)
theta(y, regime)
#> <sddr> Theta rank-mobility statistic
#>   2 regimes | 3 period transitions
#>   theta by transition:
#> [1] 0.1429 0.6364 0.3750
```
