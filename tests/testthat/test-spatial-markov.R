test_that("pooled spatial Markov equals the classic Markov matrix", {
  # Summing the lag-conditioned matrices must recover the aspatial matrix,
  # for any lag partition. This validates the transition-counting logic.
  set.seed(7)
  n <- 60; periods <- 5
  df <- do.call(rbind, lapply(seq_len(periods), function(p) {
    data.frame(id = 1:n, time = p, value = rnorm(n), lag = rnorm(n))
  }))
  m  <- markov(df, "id", "time", "value", k = 4)
  sm <- spatial_markov(df, "id", "time", "value", lag = "lag", k = 4, m = 3)

  expect_s3_class(sm, "sddr_spatial_markov")
  expect_equal(sm$pooled, m$matrix)
  expect_equal(sm$pooled_transitions, m$transitions)
})

test_that("spatial lag is the row-standardised weighted mean of neighbours", {
  # Line graph 1 - 2 - 3 (rook contiguity).
  W <- matrix(c(0, 1, 0,
                1, 0, 1,
                0, 1, 0),
              nrow = 3, byrow = TRUE, dimnames = list(1:3, 1:3))
  df <- data.frame(id = 1:3, time = 1, value = c(10, 20, 60))
  sm <- spatial_markov(df, "id", "time", "value", weights = W, k = 2, m = 2)

  # unit 1: mean(nbr = {2})    = 20
  # unit 2: mean(nbr = {1, 3}) = mean(10, 60) = 35
  # unit 3: mean(nbr = {2})    = 20
  lag_by_id <- sm$data$lagval[order(sm$data$id)]
  expect_equal(lag_by_id, c(20, 35, 20))
})

test_that("conditional matrices are row-stochastic (or empty)", {
  set.seed(3)
  n <- 20; periods <- 10
  W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
  for (i in 1:n) {
    W[i, i %% n + 1] <- 1
    W[i, (i - 2) %% n + 1] <- 1
  }
  df <- do.call(rbind, lapply(seq_len(periods), function(p) {
    data.frame(id = 1:n, time = p, value = rnorm(n))
  }))
  sm <- spatial_markov(df, "id", "time", "value", weights = W, k = 3, m = 3)

  for (P in sm$matrices) {
    rs <- rowSums(P)
    expect_true(all(abs(rs - 1) < 1e-8 | rs == 0))
  }
})

test_that("explicit breaks and lag_breaks drive the classification", {
  # 3 value classes via breaks; lag partitioned by an explicit lag break.
  df <- data.frame(
    id   = rep(1:2, each = 3),
    time = rep(1:3, times = 2),
    value = c(1, 5, 9, 9, 5, 1),
    lg    = c(2, 2, 8, 8, 2, 2)  # precomputed spatial lag
  )
  sm <- spatial_markov(df, "id", "time", "value", lag = "lg",
                       breaks = c(3, 7), lag_breaks = 5)
  expect_equal(sm$classes, 3L)       # value: (-Inf,3],(3,7],(7,Inf]
  expect_equal(sm$lag_classes, 2L)   # lag:   (-Inf,5],(5,Inf]
  # lag states: c(2,2,8,8,2,2) with break 5 -> c(1,1,2,2,1,1)
  expect_equal(sm$data$lstate[order(sm$data$id, sm$data$time)],
               c(1, 1, 2, 2, 1, 1))
})

test_that("exactly one of weights or lag is required", {
  df <- data.frame(id = 1:4, time = 1:4, value = c(1, 2, 3, 4))
  expect_error(spatial_markov(df, "id", "time", "value"), "exactly one")
  expect_error(
    spatial_markov(df, "id", "time", "value",
                   weights = diag(4), lag = "value"),
    "exactly one"
  )
})

test_that("unbalanced panels are rejected when using a weights matrix", {
  W <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(1:2, 1:2))
  df <- data.frame(id = c(1, 2, 1), time = c(1, 1, 2), value = c(5, 6, 7))
  expect_error(
    spatial_markov(df, "id", "time", "value", weights = W, k = 2),
    "balanced panel"
  )
})
