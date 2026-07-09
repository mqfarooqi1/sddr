# sddr 0.0.0.9000

* Initial scaffold.
* `markov()`: classic discrete-time, first-order Markov transition estimation
  from tidy long-format `id`/`time`/`value` panel data, with quantile-based
  (fixed or per-period) class discretisation.
* `spatial_markov()`: spatial Markov chain (Rey 2001) with transition matrices
  conditioned on the spatial-lag class of the neighbourhood. Accepts either a
  spatial weights matrix (lag computed internally) or a precomputed lag column,
  and explicit `breaks` / `lag_breaks` cut points. With matching cut points it
  reproduces PySAL `giddy`'s spatial Markov matrices to machine precision.
* `steady_state()`: ergodic / stationary distribution of a transition matrix
  or an `sddr_markov` object.
* `mfpt()`: mean first passage times (Kemeny-Snell fundamental matrix), with
  mean recurrence times on the diagonal.
* `sojourn_time()`: expected persistence time in each class.
* `tau()`: Kendall's tau rank correlation (positional / exchange mobility),
  with concordant/discordant counts and an asymptotic p-value.
* `theta()`: Theta rank-mobility statistic (Rey 2004) decomposing rank change
  by regime.
* `print()` methods for `sddr_markov` and `sddr_spatial_markov`.
* Hex logo, pkgdown site configuration, and GitHub Actions R-CMD-check CI.
