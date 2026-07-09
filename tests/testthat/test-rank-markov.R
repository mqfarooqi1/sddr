stable_panel <- function() {
  data.frame(
    id    = rep(1:3, times = 4),
    time  = rep(1:4, each = 3),
    value = rep(c(1, 2, 3), times = 4) + rep(1:4, each = 3) * 0.01
  )
}

test_that("full-rank Markov of a stable ordering is the identity", {
  m <- full_rank_markov(stable_panel(), "id", "time", "value")
  expect_s3_class(m, "sddr_markov")
  expect_equal(m$kind, "full-rank")
  expect_equal(unname(m$matrix), diag(3))
})

test_that("geo-rank Markov of a stable ordering is the identity", {
  m <- geo_rank_markov(stable_panel(), "id", "time", "value")
  expect_equal(m$kind, "geo-rank")
  expect_equal(unname(m$matrix), diag(3))
})

test_that("rank Markov matrices are row-stochastic", {
  set.seed(2)
  df <- data.frame(id = rep(1:6, times = 8), time = rep(1:8, each = 6),
                   value = rnorm(48))
  mats <- list(full_rank_markov(df, "id", "time", "value"),
               geo_rank_markov(df, "id", "time", "value"))
  for (m in mats) {
    rs <- rowSums(m$matrix)
    expect_true(all(abs(rs - 1) < 1e-8 | rs == 0))
  }
})

test_that("unbalanced panels are rejected", {
  df <- data.frame(id = c(1, 2, 1), time = c(1, 1, 2), value = c(3, 4, 5))
  expect_error(full_rank_markov(df, "id", "time", "value"), "balanced")
})
