# Changelog

## sddr 0.0.0.9000

- Initial scaffold.
- [`markov()`](https://mqfarooqi1.github.io/sddr/reference/markov.md):
  classic discrete-time, first-order Markov transition estimation from
  tidy long-format `id`/`time`/`value` panel data, with quantile-based
  (fixed or per-period) class discretisation.
- [`full_rank_markov()`](https://mqfarooqi1.github.io/sddr/reference/full_rank_markov.md)
  and
  [`geo_rank_markov()`](https://mqfarooqi1.github.io/sddr/reference/geo_rank_markov.md):
  rank-based Markov chains (ranks-as-states, and
  units-exchanging-rank-positions) — no binning required.
- [`lisa_markov()`](https://mqfarooqi1.github.io/sddr/reference/lisa_markov.md):
  Markov chain over the four Moran-scatterplot quadrants (HH/LH/LL/HL),
  capturing the joint dynamics of a unit and its neighbourhood.
- [`spatial_markov()`](https://mqfarooqi1.github.io/sddr/reference/spatial_markov.md):
  spatial Markov chain (Rey 2001) with transition matrices conditioned
  on the spatial-lag class of the neighbourhood. Accepts either a
  spatial weights matrix (lag computed internally) or a precomputed lag
  column, and explicit `breaks` / `lag_breaks` cut points. With matching
  cut points it reproduces PySAL `giddy`’s spatial Markov matrices to
  machine precision.
- [`steady_state()`](https://mqfarooqi1.github.io/sddr/reference/steady_state.md):
  ergodic / stationary distribution of a transition matrix or an
  `sddr_markov` object.
- [`mfpt()`](https://mqfarooqi1.github.io/sddr/reference/mfpt.md): mean
  first passage times (Kemeny-Snell fundamental matrix), with mean
  recurrence times on the diagonal.
- [`sojourn_time()`](https://mqfarooqi1.github.io/sddr/reference/sojourn_time.md):
  expected persistence time in each class.
- [`tau()`](https://mqfarooqi1.github.io/sddr/reference/tau.md):
  Kendall’s tau rank correlation (positional / exchange mobility), with
  concordant/discordant counts and an asymptotic p-value.
- [`theta()`](https://mqfarooqi1.github.io/sddr/reference/theta.md):
  Theta rank-mobility statistic (Rey 2004) decomposing rank change by
  regime.
- [`tau_local()`](https://mqfarooqi1.github.io/sddr/reference/tau_local.md):
  per-observation (local) Kendall’s tau decomposition.
- [`mobility()`](https://mqfarooqi1.github.io/sddr/reference/mobility.md):
  Markov mobility indices (Prais, determinant, second-eigenvalue L2, and
  the Shorrocks B1/B2 indices).
- [`print()`](https://rdrr.io/r/base/print.html) methods for
  `sddr_markov` and `sddr_spatial_markov`.
- Hex logo, pkgdown site configuration, and GitHub Actions R-CMD-check
  CI.
