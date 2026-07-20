# ============================================================
#  sddr — comprehensive demo
#  Distribution dynamics: how a cross-sectional distribution of
#  values evolves over time, and where it is heading long run.
# ============================================================

# install.packages("sddr")
library(sddr)

# ------------------------------------------------------------
# 1. Simulate demo data
#    Regional incomes on a 6x6 grid of 36 regions, observed
#    over 12 years, with spatial spillovers + mean reversion.
# ------------------------------------------------------------
set.seed(2026)

grid_n <- 6L
n      <- grid_n * grid_n          # 36 regions
years  <- 2010:2021                # 12 periods
np     <- length(years)
ids    <- seq_len(n)
coords <- expand.grid(row = 1:grid_n, col = 1:grid_n)

# Rook-contiguity spatial weights matrix (base R; ids as names).
W <- matrix(0, n, n, dimnames = list(ids, ids))
for (i in ids) for (j in ids) {
  if (i != j &&
      abs(coords$row[i] - coords$row[j]) +
      abs(coords$col[i] - coords$col[j]) == 1L) W[i, j] <- 1
}
Ws <- W / pmax(rowSums(W), 1)      # row-standardised (for the simulation)

# Spatially-correlated income process.
Y <- matrix(NA_real_, n, np, dimnames = list(ids, years))
Y[, 1] <- 100 + 15 * scale(coords$row + coords$col)[, 1]
for (t in 2:np) {
  lag <- as.numeric(Ws %*% Y[, t - 1])
  Y[, t] <- 0.55 * Y[, t - 1] + 0.25 * lag +
            0.20 * mean(Y[, t - 1]) + rnorm(n, 0, 7)
}

# Tidy long format: id / year / income  (the sddr input shape).
panel <- data.frame(
  id     = rep(ids,   times = np),
  year   = rep(years, each  = n),
  income = as.vector(Y)
)
head(panel)

# ------------------------------------------------------------
# 2. Classic Markov chain over income quintiles
# ------------------------------------------------------------
m <- markov(panel, id = "id", time = "year", value = "income", k = 5)
print(m)                              # transition matrix + steady state

# ------------------------------------------------------------
# 3. Ergodic / long-run analysis
# ------------------------------------------------------------
steady_state(m)                       # long-run class shares
mfpt(m)                               # mean first passage / recurrence times
sojourn_time(m)                       # expected periods spent per class
steady_state(m$matrix)                # also works on a bare matrix

# ------------------------------------------------------------
# 4. Markov mobility indices (how much churn?)
# ------------------------------------------------------------
mobility(m)                           # Prais, determinant, L2, Shorrocks B1/B2
mobility(m, "prais")                  # a single index

# ------------------------------------------------------------
# 5. Spatial Markov (transitions conditioned on the neighbourhood)
# ------------------------------------------------------------
sm <- spatial_markov(panel, "id", "year", "income", weights = W, k = 5)
print(sm)                             # pooled + one matrix per lag class

# Same call, but with sf-derived weights via spdep (optional):
if (requireNamespace("spdep", quietly = TRUE)) {
  listw <- spdep::mat2listw(W, style = "B")
  sm_lw <- spatial_markov(panel, "id", "year", "income", weights = listw, k = 5)
  all.equal(sm$pooled, sm_lw$pooled)  # TRUE — identical to the matrix path
}

# ------------------------------------------------------------
# 6. LISA Markov (Moran-quadrant transitions: HH / LH / LL / HL)
# ------------------------------------------------------------
lm <- lisa_markov(panel, "id", "year", "income", weights = W)
print(lm)                             # 4x4 quadrant transition matrix
lm$quadrants[1, ]                     # region 1's quadrant path over time

# ------------------------------------------------------------
# 7. Rank-based Markov chains (no binning needed)
# ------------------------------------------------------------
fr <- full_rank_markov(panel, "id", "year", "income")  # rank = state (n x n)
gr <- geo_rank_markov(panel,  "id", "year", "income")  # unit per rank slot
dim(fr$matrix)
head(steady_state(gr))

# ------------------------------------------------------------
# 8. Rank-mobility measures (positional / exchange mobility)
# ------------------------------------------------------------
first <- Y[, 1]; last <- Y[, np]
tau(first, last)                      # Kendall's tau, first vs last year
head(tau_local(first, last)$tau_local)  # per-region contribution

# theta(): within-regime rank mobility. Regimes = north vs south half.
regime <- ifelse(coords$row <= grid_n / 2, "north", "south")
theta(Y, regime)                      # theta per year-to-year transition
