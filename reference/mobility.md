# Markov mobility indices

Scalar summaries of how much mixing a transition matrix implies — larger
values mean more mobility (a faster-churning distribution), values near
zero mean a near-immobile system. Several classical indices are
provided:

- `prais`:

  Prais index \\(k - \mathrm{tr}(P)) / (k - 1)\\ — average probability
  of leaving a class.

- `determinant`:

  \\1 - \|\det(P)\|\\.

- `L2`:

  \\1 - \|\lambda_2\|\\, one minus the second-largest eigenvalue modulus
  (the mixing rate).

- `shorrock1`:

  \\(k - k \sum_i \pi_i p\_{ii}) / (k - 1)\\, weighted by an initial
  distribution `ini`.

- `shorrock2`:

  \\\sum_i \sum_j \pi_i p\_{ij} \|i - j\| / (k - 1)\\, penalising longer
  jumps.

## Usage

``` r
mobility(x, measure = "all", ini = NULL)
```

## Arguments

- x:

  A row-stochastic matrix, or an object of class `sddr_markov`.

- measure:

  One of `"prais"`, `"determinant"`, `"L2"`, `"shorrock1"`,
  `"shorrock2"`, or `"all"` (default) to return every index.

- ini:

  Optional initial distribution (length `k`) for the Shorrock indices;
  defaults to uniform.

## Value

A named numeric vector of the requested index or indices.

## References

Prais, S. J. (1955). Measuring social mobility. *JRSS A*, 118(1), 56-66.
Shorrocks, A. F. (1978). The measurement of mobility. *Econometrica*,
46(5), 1013-1024. [doi:10.2307/1911433](https://doi.org/10.2307/1911433)

## Examples

``` r
P <- matrix(c(0.7, 0.2, 0.1,
              0.2, 0.6, 0.2,
              0.1, 0.3, 0.6), nrow = 3, byrow = TRUE)
mobility(P)
#>       prais determinant          L2   shorrock1   shorrock2 
#>   0.5500000   0.8100000   0.4381966   0.5500000   0.2166667 
mobility(P, "prais")
#> prais 
#>  0.55 
```
