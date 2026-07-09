ring_weights <- function(n) {
  W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
  for (i in 1:n) {
    W[i, i %% n + 1] <- 1
    W[i, (i - 2) %% n + 1] <- 1
  }
  W
}

test_that("spatial_markov accepts an spdep listw and matches the matrix path", {
  skip_if_not_installed("spdep")
  set.seed(5)
  n <- 9L; periods <- 12L
  W <- ring_weights(n)
  listw <- spdep::mat2listw(W, style = "B")
  df <- data.frame(id = rep(1:n, times = periods),
                   time = rep(1:periods, each = n), value = rnorm(n * periods))

  sm_mat <- spatial_markov(df, "id", "time", "value", weights = W, k = 3)
  sm_lw  <- spatial_markov(df, "id", "time", "value", weights = listw, k = 3)

  expect_equal(sm_mat$pooled, sm_lw$pooled)
  for (i in seq_along(sm_mat$matrices)) {
    expect_equal(sm_mat$matrices[[i]], sm_lw$matrices[[i]])
  }
})

test_that("lisa_markov accepts listw and nb objects", {
  skip_if_not_installed("spdep")
  set.seed(5)
  n <- 9L; periods <- 12L
  W <- ring_weights(n)
  listw <- spdep::mat2listw(W, style = "B")
  nb <- listw$neighbours
  df <- data.frame(id = rep(1:n, times = periods),
                   time = rep(1:periods, each = n), value = rnorm(n * periods))

  m_mat <- lisa_markov(df, "id", "time", "value", weights = W)
  m_lw  <- lisa_markov(df, "id", "time", "value", weights = listw)
  m_nb  <- lisa_markov(df, "id", "time", "value", weights = nb)

  expect_equal(m_mat$matrix, m_lw$matrix)
  expect_equal(m_mat$matrix, m_nb$matrix)
})

test_that("invalid weights specifications are rejected", {
  df <- data.frame(id = rep(1:3, 2), time = rep(1:2, each = 3), value = rnorm(6))
  expect_error(
    lisa_markov(df, "id", "time", "value", weights = "not-weights"),
    "matrix"
  )
})
