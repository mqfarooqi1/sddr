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

test_that("spatial Markov matches PySAL giddy with identical cutoffs", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")
  have_libpysal <- tryCatch({ reticulate::import("libpysal"); TRUE },
                            error = function(e) FALSE)
  skip_if_not(have_libpysal, "libpysal not available")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")
  libpysal <- reticulate::import("libpysal")

  set.seed(11)
  nr <- 6L; nc <- 6L; n <- nr * nc; periods <- 15L
  w <- libpysal$weights$lat2W(nr, nc)
  w$transform <- "r"
  y <- matrix(rnorm(n * periods), n, periods)

  wf <- w$full(); Wmat <- wf[[1]]; ids <- as.integer(wf[[2]])
  rsum <- rowSums(Wmat); rsum[rsum == 0] <- 1; Wmat <- Wmat / rsum
  dimnames(Wmat) <- list(ids, ids)

  lagm <- Wmat %*% y
  cc <- as.numeric(stats::quantile(as.vector(y),    c(.25, .5, .75)))
  lc <- as.numeric(stats::quantile(as.vector(lagm), c(.25, .5, .75)))

  Sp <- giddy$markov$Spatial_Markov(
    np$array(y), w, cutoffs = np$array(cc), lag_cutoffs = np$array(lc),
    fixed = TRUE, variable_name = "y"
  )
  gP <- Sp$P

  dfl <- data.frame(id = rep(ids, times = periods),
                    time = rep(seq_len(periods), each = n),
                    value = as.vector(y))
  smk <- spatial_markov(dfl, "id", "time", "value", weights = Wmat,
                        k = 4L, m = 4L, breaks = cc, lag_breaks = lc,
                        row_standardize = FALSE)

  for (i in seq_len(dim(gP)[1])) {
    expect_equal(gP[i, , ], unname(smk$matrices[[i]]), tolerance = 1e-9)
  }
})

test_that("ergodic quantities match PySAL giddy", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")

  set.seed(5)
  k <- 5L
  P <- matrix(runif(k * k), k, k); P <- P / rowSums(P)

  expect_equal(as.numeric(steady_state(P)),
               as.numeric(giddy$ergodic$steady_state(np$array(P))),
               tolerance = 1e-9)
  expect_equal(unname(mfpt(P)),
               unname(as.matrix(giddy$ergodic$mfpt(np$array(P)))),
               tolerance = 1e-7)
  expect_equal(as.numeric(sojourn_time(P)),
               as.numeric(giddy$markov$sojourn_time(np$array(P))),
               tolerance = 1e-9)
})

test_that("Tau and Theta match PySAL giddy", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")

  set.seed(3)
  n <- 40L
  x <- rnorm(n); y <- x + rnorm(n)
  gt <- giddy$rank$Tau(np$array(x), np$array(y))
  st <- tau(x, y)
  expect_equal(st$tau, as.numeric(gt$tau), tolerance = 1e-9)
  expect_equal(st$concordant, as.numeric(gt$concordant))
  expect_equal(st$discordant, as.numeric(gt$discordant))

  ym <- matrix(rnorm(40), nrow = 10)
  reg <- rep(c(0L, 1L), each = 5)
  gth <- giddy$rank$Theta(np$array(ym), np$array(reg), permutations = 0L)
  expect_equal(theta(ym, reg)$theta, as.numeric(gth$theta), tolerance = 1e-9)
})

test_that("mobility indices and local tau match PySAL giddy", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")

  set.seed(4)
  k <- 5L
  P <- matrix(runif(k * k), k, k); P <- P / rowSums(P)
  mm <- giddy$mobility$markov_mobility
  codes <- c(prais = "P", determinant = "D", L2 = "L2",
             shorrock1 = "B1", shorrock2 = "B2")
  for (nm in names(codes)) {
    expect_equal(mobility(P, nm)[[1]],
                 as.numeric(mm(np$array(P), measure = codes[[nm]])),
                 tolerance = 1e-7)
  }

  n <- 30L
  x <- rnorm(n); y <- x + rnorm(n)
  gt <- giddy$rank$Tau_Local(np$array(x), np$array(y))
  st <- tau_local(x, y)
  expect_equal(st$tau_local, as.numeric(gt$tau_local), tolerance = 1e-9)
  expect_equal(st$tau, as.numeric(gt$tau), tolerance = 1e-9)
})

test_that("full-rank and geo-rank Markov match PySAL giddy", {
  skip_if_not(giddy_available(), "PySAL giddy not importable via reticulate")

  giddy <- reticulate::import("giddy")
  np <- reticulate::import("numpy")

  set.seed(7)
  n <- 8L; periods <- 14L
  y <- matrix(rnorm(n * periods), n, periods)
  df <- data.frame(id = rep(1:n, times = periods),
                   time = rep(1:periods, each = n),
                   value = as.vector(y))

  frm <- giddy$markov$FullRank_Markov(np$array(y), summary = FALSE)
  grm <- giddy$markov$GeoRank_Markov(np$array(y), summary = FALSE)

  expect_equal(unname(full_rank_markov(df, "id", "time", "value")$matrix),
               unname(as.matrix(frm$p)), tolerance = 1e-9)
  expect_equal(unname(geo_rank_markov(df, "id", "time", "value")$matrix),
               unname(as.matrix(grm$p)), tolerance = 1e-9)
})
