# sddr

<!-- badges: start -->
<!-- badges: end -->

**The complete toolkit for distribution dynamics in R** — analysing how a
distribution of values across units (regions, firms, markets, assets) *evolves
over time*, and where it is heading in the long run.

`sddr` brings distribution-dynamics analysis into one tidy, long-format
framework and pushes it **beyond what any existing package offers** — with
continuous stochastic kernels, continuous-time transitions, and modern
inference that current R (and Python) tooling simply does not provide.

## What makes `sddr` different — the new ideas

These are the capabilities that set `sddr` apart. (Status: ✅ available now ·
🔜 in active development.)

- 🔜 **Continuous stochastic kernels** — density-based distribution dynamics
  (Quah), estimating the full conditional density with *no* arbitrary class
  binning. The modern approach, almost absent from R.
- 🔜 **Continuous-time & irregular-interval Markov** — for panels not observed
  on a regular time grid.
- 🔜 **Modern inference** — bootstrap and Bayesian intervals for transition
  probabilities, plus formal tests for Markov order, time-stationarity, and
  *whether space matters at all*.
- 🔜 **Change-point / non-stationarity detection** in transition dynamics.
- 🔜 **Simulation engine, animated & interactive visualisation, `broom`
  tidiers**, and an optional compiled backend for large panels.

## Established methods, done right

A complete, tidy, `sf`-friendly implementation of the classical toolkit:

- ✅ **Classic Markov** chains and **ergodic** (steady-state) analysis
- ✅ **Spatial Markov** — transitions conditioned on neighbourhood context
- 🔜 LISA / full-rank / geo-rank Markov · rank & exchange mobility (Tau family)
  · mobility indices (Prais, Shorrocks) · directional LISA · sequence analysis
  · mean first-passage & sojourn times

All from long-format `id` / `time` / `value` data — no transition-matrix
bookkeeping.

## Installation

``` r
# install.packages("pak")
pak::pak("mqfarooqi1/sddr")
```

## Quick start

``` r
library(sddr)

df <- data.frame(
  id    = rep(1:200, each = 5),
  time  = rep(2000:2004, times = 200),
  value = rnorm(1000)
)

# Classic Markov chain over distribution quintiles.
m <- markov(df, id = "id", time = "time", value = "value", k = 5)
m
steady_state(m)              # long-run distribution
```

## Roadmap

| Release | Focus |
|--------|--------|
| **v0.1** | The complete classical toolkit in one tidy API (Markov, spatial Markov, rank mobility, ergodic, sequences) |
| v0.2+ | The new ideas: continuous **stochastic kernels**, continuous-time Markov, modern inference, change-point detection |
| later | Compiled backend, simulation engine, animated/interactive visualisation |

## Validation & prior work

`sddr`'s methods build directly on the distribution-dynamics literature —
Quah (1993) for distributional convergence and Rey (2001) for spatial Markov
dynamics.

As a correctness guarantee, where methods overlap with **PySAL's `giddy`** (the
established reference implementation) `sddr` is checked for numerical parity:
the classic and spatial Markov estimators currently reproduce `giddy` to
machine precision. `sddr` extends well past `giddy` and the R package `griddy`
with the continuous-kernel, continuous-time, and inference methods above.

## License

MIT © Muhammad Farooqi
