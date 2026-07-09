test_that("sojourn_time is 1 / (1 - p_ii)", {
  P <- matrix(c(0.8, 0.2,
                0.3, 0.7), nrow = 2, byrow = TRUE)
  expect_equal(unname(sojourn_time(P)), c(5, 10 / 3))
})

test_that("mfpt matches the closed form for a two-state chain", {
  P <- matrix(c(0.8, 0.2,
                0.3, 0.7), nrow = 2, byrow = TRUE)
  M <- mfpt(P)
  # Diagonal = mean recurrence time = 1 / pi, with pi = (0.6, 0.4).
  expect_equal(diag(M), c(1 / 0.6, 1 / 0.4))
  # Off-diagonal of a two-state chain = 1 / p_ij.
  expect_equal(M[1, 2], 1 / 0.2)
  expect_equal(M[2, 1], 1 / 0.3)
})

test_that("mfpt and sojourn_time dispatch on sddr_markov objects", {
  set.seed(1)
  df <- data.frame(
    id   = rep(1:40, each = 4),
    time = rep(1:4, times = 40),
    value = rnorm(160)
  )
  m <- markov(df, "id", "time", "value", k = 3)
  expect_equal(mfpt(m), mfpt(m$matrix))
  expect_equal(sojourn_time(m), sojourn_time(m$matrix))
})
