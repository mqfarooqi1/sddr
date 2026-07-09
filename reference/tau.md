# Kendall's tau rank correlation

Kendall's \\\tau\\ measures the rank concordance between two variables —
for distribution dynamics, typically a set of units' values at two
points in time. It is the (concordant minus discordant) pair count,
normalised, and is a natural summary of *positional* (exchange)
mobility: \\\tau = 1\\ means ranks are perfectly preserved, \\\tau =
-1\\ a complete reversal, and \\\tau \approx 0\\ extensive rank
shuffling. Ties are handled with the tau-b correction.

## Usage

``` r
tau(x, y)
```

## Arguments

- x, y:

  Numeric vectors of equal length.

## Value

An object of class `sddr_tau`: a list with the statistic `tau`, the
`concordant` and `discordant` pair counts, an asymptotic (normal)
two-sided `p_value`, and `n`.

## References

Kendall, M. G. (1938). A new measure of rank correlation. *Biometrika*,
30(1/2), 81-93. [doi:10.2307/2332226](https://doi.org/10.2307/2332226)

## See also

[`theta()`](https://mqfarooqi1.github.io/sddr/reference/theta.md)

## Examples

``` r
set.seed(1)
x <- rnorm(30)
y <- x + rnorm(30)          # correlated ranks
tau(x, y)
#> <sddr> Kendall's tau
#>   tau = 0.5126  (concordant 329, discordant 106, n = 30)
#>   asymptotic two-sided p-value = 6.934e-05
```
