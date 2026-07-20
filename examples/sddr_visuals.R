# ============================================================
#  sddr — visual analysis (maps + heatmaps) with ggplot2
#  Saves 6 figures covering every sddr method.
# ============================================================
library(sddr)
library(ggplot2)

OUT <- "sddr_figures"; dir.create(OUT, showWarnings = FALSE)  # figures saved here

# ---- demo data: regional incomes on a 6x6 grid, 12 years ----
set.seed(2026)
grid_n <- 6L; n <- grid_n^2; years <- 2010:2021; np <- length(years)
ids <- seq_len(n); coords <- expand.grid(row = 1:grid_n, col = 1:grid_n)
W <- matrix(0, n, n, dimnames = list(ids, ids))
for (i in ids) for (j in ids)
  if (i != j && abs(coords$row[i]-coords$row[j]) + abs(coords$col[i]-coords$col[j]) == 1L)
    W[i, j] <- 1
Ws <- W / pmax(rowSums(W), 1)
Y <- matrix(NA_real_, n, np, dimnames = list(ids, years))
Y[, 1] <- 100 + 15 * scale(coords$row + coords$col)[, 1]
for (t in 2:np) Y[, t] <- 0.55*Y[, t-1] + 0.25*as.numeric(Ws %*% Y[, t-1]) +
                          0.20*mean(Y[, t-1]) + rnorm(n, 0, 7)
panel <- data.frame(id = rep(ids, np), year = rep(years, each = n),
                    income = as.vector(Y))

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid = element_blank(),
                plot.title = element_text(face = "bold")))
# tidy a matrix -> long data for heatmaps
heat <- function(mat)
  transform(as.data.frame.table(mat, stringsAsFactors = FALSE),
            from = Var1, to = Var2, v = Freq)

# ---- 1. Income choropleth maps over time --------------------
mp <- merge(panel, cbind(id = ids, coords), by = "id")
g1 <- ggplot(mp, aes(col, row, fill = income)) +
  geom_tile(color = "white", linewidth = .3) +
  facet_wrap(~year, nrow = 3) + coord_equal() + scale_y_reverse() +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Regional income over space and time",
       subtitle = "Each panel = one year; the distribution shifts and clusters spatially",
       x = NULL, y = NULL, fill = "income")
ggsave(file.path(OUT, "fig1_income_maps.png"), g1, width = 10, height = 6.5, dpi = 120)

# ---- 2. Classic Markov transition matrix heatmap ------------
m <- markov(panel, "id", "year", "income", k = 5)
g2 <- ggplot(heat(m$matrix),
             aes(to, factor(from, levels = rev(sort(unique(from)))), fill = v)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", v)), size = 4) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
  labs(title = "Classic Markov transition matrix (income quintiles)",
       subtitle = "P(next-year class | this-year class); strong diagonal = persistence",
       x = "to class", y = "from class", fill = "prob")
ggsave(file.path(OUT, "fig2_transition.png"), g2, width = 6.6, height = 5.4, dpi = 120)

# ---- 3. Spatial Markov: conditional matrices ----------------
sm <- spatial_markov(panel, "id", "year", "income", weights = W, k = 5)
mats <- c(list(Pooled = sm$pooled), sm$matrices)
smdf <- do.call(rbind, Map(function(nm, mat) cbind(panel = nm, heat(mat)),
                           names(mats), mats))
smdf$panel <- factor(smdf$panel, levels = c("Pooled", names(sm$matrices)))
g3 <- ggplot(smdf, aes(to, factor(from, levels = rev(sort(unique(from)))), fill = v)) +
  geom_tile(color = "white") + facet_wrap(~panel, nrow = 2) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
  labs(title = "Spatial Markov: transitions conditioned on the neighbourhood",
       subtitle = "If space did not matter, every panel would equal 'Pooled'",
       x = "to class", y = "from class", fill = "prob")
ggsave(file.path(OUT, "fig3_spatial_markov.png"), g3, width = 11, height = 6, dpi = 120)

# ---- 4. LISA Moran-quadrant maps over time ------------------
lm <- lisa_markov(panel, "id", "year", "income", weights = W)
qlab <- c("HH", "LH", "LL", "HL")
qdf <- merge(data.frame(id = rep(ids, np), year = rep(years, each = n),
                        quad = qlab[as.vector(lm$quadrants)]),
             cbind(id = ids, coords), by = "id")
g4 <- ggplot(subset(qdf, year %in% years[seq(1, np, 2)]),
             aes(col, row, fill = factor(quad, levels = qlab))) +
  geom_tile(color = "white", linewidth = .3) +
  facet_wrap(~year, nrow = 2) + coord_equal() + scale_y_reverse() +
  scale_fill_manual(values = c(HH="#d7191c", LH="#abd9e9", LL="#2c7bb6", HL="#fdae61")) +
  labs(title = "LISA Moran-scatterplot quadrant of each region",
       subtitle = "HH/LL = spatial clusters, LH/HL = spatial outliers; tracked by lisa_markov()",
       x = NULL, y = NULL, fill = "quadrant")
ggsave(file.path(OUT, "fig4_lisa_maps.png"), g4, width = 9.5, height = 5.2, dpi = 120)

# ---- 5. Long-run dynamics: steady state, mobility, theta ----
regime <- ifelse(coords$row <= grid_n/2, "north", "south")
ss <- steady_state(m); mob <- mobility(m); th <- theta(Y, regime)$theta
dyn <- rbind(
  data.frame(panel = "Long-run class shares", x = names(ss),  y = as.numeric(ss)),
  data.frame(panel = "Mobility indices",       x = names(mob), y = as.numeric(mob)),
  data.frame(panel = "Theta by year",          x = names(th),  y = as.numeric(th)))
dyn$panel <- factor(dyn$panel, levels = unique(dyn$panel))
g5 <- ggplot(dyn, aes(x, y, fill = panel)) + geom_col(show.legend = FALSE) +
  facet_wrap(~panel, scales = "free", nrow = 1) +
  scale_fill_viridis_d(end = .8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Long-run behaviour & mobility summaries", x = NULL, y = NULL)
ggsave(file.path(OUT, "fig5_dynamics.png"), g5, width = 11, height = 4, dpi = 120)

# ---- 6. Rank mobility: first vs last year -------------------
first <- Y[, 1]; last <- Y[, np]
rk <- data.frame(r0 = rank(-first), r1 = rank(-last),
                 local_tau = tau_local(first, last)$tau_local)
g6 <- ggplot(rk, aes(r0, r1, color = local_tau)) +
  geom_abline(slope = 1, linetype = 2, color = "grey60") +
  geom_point(size = 3) +
  scale_color_gradient2(low = "#d7191c", mid = "grey85", high = "#2c7bb6", midpoint = 0) +
  scale_x_reverse() + scale_y_reverse() + coord_equal() +
  labs(title = sprintf("Positional (rank) mobility  —  Kendall's tau = %.2f",
                       tau(first, last)$tau),
       x = "rank in 2010 (1 = highest)", y = "rank in 2021", color = "local tau")
ggsave(file.path(OUT, "fig6_rank_mobility.png"), g6, width = 6.4, height = 5.6, dpi = 120)

cat("Saved 6 figures to", normalizePath(OUT), "\n")
