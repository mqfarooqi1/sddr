# Numerical parity against the reference PySAL `giddy` implementation.
#
# These tests are skipped unless `reticulate` and a Python with `giddy` are
# available, so the package checks cleanly on machines (and on CRAN) without a
# Python toolchain. When giddy *is* importable, they assert exact agreement.

skip_on_cran()
skip_if_not_installed("reticulate")

giddy_available <- function() {
  tryCatch({
    reticulate::import("giddy", delay_load = FALSE)
    reticulate::import("numpy", delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

test_that("classic Markov transition matrix matches PySAL giddy exactly", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")

  set.seed(1)
  n <- 50L; periods <- 8L; k <- 4L
  classes <- matrix(sample(seq_len(k), n * periods, replace = TRUE),
                    nrow = n, ncol = periods)

  # giddy: rows = units, cols = periods, entries = integer states.
  gm <- giddy$markov$Markov(np$array(classes))
  gp <- unname(as.matrix(gm$p))

  # sddr: long format; map integer class c -> state c via fixed breaks.
  df <- data.frame(
    id    = rep(seq_len(n), times = periods),
    time  = rep(seq_len(periods), each = n),
    value = as.vector(classes)
  )
  sm <- markov(df, "id", "time", "value", breaks = seq(1.5, k - 0.5, by = 1))

  expect_equal(unname(sm$matrix), gp, tolerance = 1e-9)
  expect_equal(unname(sm$transitions), unname(as.matrix(gm$transitions)),
               tolerance = 1e-9)
})
