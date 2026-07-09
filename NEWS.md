# sddr 0.0.0.9000

* Initial scaffold.
* `markov()`: classic discrete-time, first-order Markov transition estimation
  from tidy long-format `id`/`time`/`value` panel data, with quantile-based
  (fixed or per-period) class discretisation.
* `steady_state()`: ergodic / stationary distribution of a transition matrix
  or an `sddr_markov` object.
* `print()` method for `sddr_markov`.
