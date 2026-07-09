test_that("tau is +1 for identical rankings and -1 for reversed", {
  t1 <- tau(1:5, 1:5)
  expect_s3_class(t1, "sddr_tau")
  expect_equal(t1$tau, 1)
  expect_equal(t1$concordant, 10)
  expect_equal(t1$discordant, 0)

  t2 <- tau(1:5, 5:1)
  expect_equal(t2$tau, -1)
  expect_equal(t2$discordant, 10)
})

test_that("tau agrees with stats::cor(method = 'kendall')", {
  set.seed(2)
  x <- rnorm(50); y <- x + rnorm(50)
  expect_equal(tau(x, y)$tau, cor(x, y, method = "kendall"))
})

test_that("theta is 0 for within-regime exchange, 1 for coherent regime moves", {
  regime <- c("a", "a", "b", "b")
  y_within <- cbind(c(1, 2, 3, 4), c(2, 1, 4, 3))   # swaps inside each regime
  expect_equal(theta(y_within, regime)$theta, 0)

  y_coher <- cbind(c(1, 2, 3, 4), c(3, 4, 1, 2))    # regimes move together
  expect_equal(theta(y_coher, regime)$theta, 1)
})

test_that("theta values lie in [0, 1]", {
  set.seed(9)
  y <- matrix(rnorm(60), nrow = 15)
  reg <- rep(1:3, each = 5)
  th <- theta(y, reg)$theta
  expect_true(all(th >= 0 & th <= 1))
})
