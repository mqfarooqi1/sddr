test_that("classic markov counts and normalises transitions correctly", {
  # Two classes split at 10. Construct known state sequences:
  #   a: 5, 5, 15, 15  -> states 1,1,2,2 -> transitions 1->1, 1->2, 2->2
  #   b: 15, 5, 5, 15  -> states 2,1,1,2 -> transitions 2->1, 1->1, 1->2
  df <- data.frame(
    id    = rep(c("a", "b"), each = 4),
    time  = rep(1:4, times = 2),
    value = c(5, 5, 15, 15,
              15, 5, 5, 15)
  )
  m <- markov(df, id = "id", time = "time", value = "value", breaks = 10)

  expect_s3_class(m, "sddr_markov")
  expect_equal(m$classes, 2L)

  # Observed transition counts.
  expect_equal(m$transitions[1, 1], 2)  # 1->1  (a, b)
  expect_equal(m$transitions[1, 2], 2)  # 1->2  (a, b)
  expect_equal(m$transitions[2, 1], 1)  # 2->1  (b)
  expect_equal(m$transitions[2, 2], 1)  # 2->2  (a)

  # Row-stochastic probability matrix.
  expect_equal(unname(m$matrix[1, ]), c(0.5, 0.5))
  expect_equal(unname(m$matrix[2, ]), c(0.5, 0.5))
  expect_equal(unname(rowSums(m$matrix)), c(1, 1))
})

test_that("steady state solves pi P = pi and sums to one", {
  P <- matrix(c(0.8, 0.2,
                0.3, 0.7), nrow = 2, byrow = TRUE)
  ss <- steady_state(P)

  expect_equal(sum(ss), 1)
  expect_equal(as.numeric(ss %*% P), as.numeric(ss))
  # Analytic stationary distribution of this chain is (0.6, 0.4).
  expect_equal(unname(ss), c(0.6, 0.4), tolerance = 1e-8)
})

test_that("quantile classification produces the requested number of classes", {
  set.seed(42)
  df <- data.frame(
    id    = rep(1:100, each = 3),
    time  = rep(1:3, times = 100),
    value = rnorm(300)
  )
  m <- markov(df, id = "id", time = "time", value = "value", k = 5)
  expect_equal(dim(m$matrix), c(5L, 5L))
  expect_true(all(abs(rowSums(m$matrix) - 1) < 1e-8))
})

test_that("missing values are rejected", {
  df <- data.frame(id = 1:4, time = 1:4, value = c(1, NA, 3, 4))
  expect_error(
    markov(df, id = "id", time = "time", value = "value", k = 2),
    "missing values"
  )
})
