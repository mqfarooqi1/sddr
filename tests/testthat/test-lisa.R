test_that("lisa_markov returns a 4x4 quadrant chain", {
  set.seed(1)
  n <- 9; periods <- 8
  W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
  for (i in 1:n) {
    W[i, i %% n + 1] <- 1
    W[i, (i - 2) %% n + 1] <- 1
  }
  df <- data.frame(id = rep(1:n, times = periods),
                   time = rep(1:periods, each = n), value = rnorm(n * periods))
  m <- lisa_markov(df, "id", "time", "value", weights = W)

  expect_s3_class(m, "sddr_markov")
  expect_equal(m$kind, "LISA")
  expect_equal(dim(m$matrix), c(4L, 4L))
  expect_equal(rownames(m$matrix), c("HH", "LH", "LL", "HL"))
  rs <- rowSums(m$matrix)
  expect_true(all(abs(rs - 1) < 1e-8 | rs == 0))
})

test_that("two anti-correlated units sit in the HL and LH quadrants", {
  W <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(1:2, 1:2))
  # unit 1 always the higher value; its single neighbour is always lower.
  df <- data.frame(id = rep(1:2, 3), time = rep(1:3, each = 2),
                   value = c(10, 2, 8, 1, 12, 3))
  m <- lisa_markov(df, "id", "time", "value", weights = W)

  expect_equal(m$quadrants[1, ], c(4L, 4L, 4L))   # HL
  expect_equal(m$quadrants[2, ], c(2L, 2L, 2L))   # LH
  expect_equal(m$matrix["HL", "HL"], 1)
  expect_equal(m$matrix["LH", "LH"], 1)
})
