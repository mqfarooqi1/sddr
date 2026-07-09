# Local Kendall's tau (local indicator of mobility association)

Decomposes Kendall's tau into a per-observation contribution. Each
unit's local tau is the average sign of its pairwise rank concordance
with all other units, \\\sum\_{j \ne i} S\_{ij} / (n - 1)\\, where
\\S\_{ij} = \mathrm{sign}((x_i - x_j)(y_i - y_j))\\. Units with local
tau near 1 keep their relative position; near -1 they move against the
grain. The overall (tau-a) statistic is the mean of the off-diagonal
`S`.

## Usage

``` r
tau_local(x, y)
```

## Arguments

- x, y:

  Numeric vectors of equal length.

## Value

An object of class `sddr_tau_local`: a list with the per-observation
`tau_local` vector, the global `tau` (tau-a), the sign matrix `S`, and
`n`.

## See also

[`tau()`](https://mqfarooqi1.github.io/sddr/reference/tau.md)

## Examples

``` r
set.seed(1)
x <- rnorm(20); y <- x + rnorm(20)
tau_local(x, y)$tau_local
#>  [1]  0.05263158  0.15789474  0.89473684 -0.36842105  0.15789474  0.89473684
#>  [7]  0.36842105 -0.15789474  0.05263158  0.26315789  0.89473684  0.26315789
#> [13]  0.36842105  1.00000000 -0.26315789  0.57894737  0.57894737  0.36842105
#> [19]  0.68421053  0.57894737
```
