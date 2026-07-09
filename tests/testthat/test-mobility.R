test_that("prais and determinant are 0 for the identity (no mobility)", {
  I2 <- diag(2)
  expect_equal(mobility(I2, "prais")[["prais"]], 0)
  expect_equal(mobility(I2, "determinant")[["determinant"]], 0)
})

test_that("mobility indices match closed forms for a symmetric 2x2", {
  P <- matrix(0.5, 2, 2)
  expect_equal(mobility(P, "prais")[["prais"]], 1)
  expect_equal(mobility(P, "determinant")[["determinant"]], 1)  # det(P) = 0
  # shorrock2 = sum_ij pi_i p_ij |i-j| / (k-1) = 0.5 for uniform ini.
  expect_equal(mobility(P, "shorrock2")[["shorrock2"]], 0.5)
})

test_that("mobility('all') returns all five indices", {
  P <- matrix(c(0.7, 0.2, 0.1,
                0.2, 0.6, 0.2,
                0.1, 0.3, 0.6), nrow = 3, byrow = TRUE)
  m <- mobility(P)
  expect_named(m, c("prais", "determinant", "L2", "shorrock1", "shorrock2"))
  expect_length(m, 5)
})

test_that("mobility dispatches on an sddr_markov object", {
  set.seed(1)
  df <- data.frame(id = rep(1:50, each = 4), time = rep(1:4, times = 50),
                   value = rnorm(200))
  m <- markov(df, "id", "time", "value", k = 3)
  expect_equal(mobility(m), mobility(m$matrix))
})
