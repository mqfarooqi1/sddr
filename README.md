# sddr

<!-- badges: start -->
<!-- badges: end -->

**Spatial distribution dynamics for R** — how a cross-sectional distribution of
values (incomes, prices, rates, indices) evolves over time.

`sddr` is a tidy, long-format toolkit for distribution-dynamics analysis. It is
built to be a **superset** of the methods in PySAL's
[`giddy`](https://pysal.org/giddy/) and the R package
[`griddy`](https://github.com/dshkol/griddy), and to go beyond both with
continuous stochastic kernels, continuous-time transitions, and modern
inference in later releases.

Unlike matrix-position APIs, `sddr` works directly from long panel data with
explicit `id`, `time`, and `value` columns.

## Installation

``` r
# development version
# install.packages("pak")
pak::pak("mqfarooqi1/sddr")
```

## Quick start

``` r
library(sddr)

# Long-format panel: one row per unit per period.
df <- data.frame(
  id    = rep(1:200, each = 5),
  time  = rep(2000:2004, times = 200),
  value = rnorm(1000)
)

# Classic Markov chain over relative-distribution quintiles.
m <- markov(df, id = "id", time = "time", value = "value", k = 5)
m

# Long-run (ergodic) distribution.
steady_state(m)
```

## Roadmap

| Phase | Contents |
|------|-----------|
| **1 (now)** | Classic Markov, ergodic steady-state, tidy design |
| 2 | Spatial Markov, LISA / rank Markov, Tau family, mobility indices, sequences (full `giddy` parity) |
| 3 | Continuous **stochastic kernels**, continuous-time Markov, modern inference (bootstrap/Bayes, order & stationarity & spatial-dependence tests) |
| 4 | Compiled backend, simulation engine, animated/interactive visualisation |

Every shared method is validated for numerical parity against the reference
`giddy` implementation.

## License

MIT © Muhammad Farooqi
