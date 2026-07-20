# ============================================================
#  sddr on REAL data & a REAL world map
#  Life expectancy across ~125 countries, 1952-2007 (gapminder),
#  on spData::world geometry, with k-nearest-neighbour weights.
#
#  Extra packages (beyond sddr): ggplot2, sf, spData, gapminder, spdep
#    install.packages(c("ggplot2","sf","spData","gapminder","spdep"))
# ============================================================
suppressMessages({
  library(sddr); library(ggplot2); library(sf)
  library(spData); library(gapminder); library(spdep)
})
OUT <- "sddr_figures"; dir.create(OUT, showWarnings = FALSE)
sf_use_s2(FALSE)

# ---- 1. Join real panel data to real country geometry -------
w  <- spData::world
gm <- subset(gapminder, as.character(country) %in% w$name_long)
gm$country <- as.character(gm$country)
countries <- sort(unique(gm$country))            # 125 countries
yrs       <- sort(unique(gm$year))               # 12 five-year periods
geom <- w[match(countries, w$name_long), c("name_long", "continent")]

# lifeExp as a country x year matrix (for tau / theta)
Y <- matrix(NA_real_, length(countries), length(yrs),
            dimnames = list(countries, yrs))
Y[cbind(match(gm$country, countries), match(gm$year, yrs))] <- gm$lifeExp

# tidy long panel (the sddr input)
panel <- data.frame(country = gm$country, year = gm$year, lifeExp = gm$lifeExp)

# ---- 2. Spatial weights: 5 nearest neighbours (great-circle) ----
xy <- suppressWarnings(st_coordinates(st_centroid(st_geometry(geom))))
nb <- knn2nb(knearneigh(xy, k = 5, longlat = TRUE), row.names = countries)
lw <- nb2listw(nb, style = "W")

# ---- 3. sddr analysis on the real panel ---------------------
m   <- markov(panel, "country", "year", "lifeExp", k = 5)
sm  <- spatial_markov(panel, "country", "year", "lifeExp", weights = lw, k = 5)
lm  <- lisa_markov(panel, "country", "year", "lifeExp", weights = lw)
mob <- mobility(m)
th  <- theta(Y, geom$continent)$theta          # within-continent mobility
tau(Y[, 1], Y[, ncol(Y)])                       # 1952 vs 2007 rank correlation

theme_set(theme_minimal(base_size = 12) +
          theme(plot.title = element_text(face = "bold")))
heat <- function(mat)
  transform(as.data.frame.table(mat, stringsAsFactors = FALSE),
            from = Var1, to = Var2, v = Freq)
robin <- "+proj=robin +datum=WGS84"
bg <- st_transform(st_geometry(w), robin)

# ---- FIG 1: life-expectancy world maps over time ------------
sel <- c(1952, 1972, 1992, 2007)
wl <- do.call(rbind, lapply(sel, function(yr) {
  g <- geom; g$year <- yr; g$lifeExp <- Y[, as.character(yr)]; g
}))
wl <- st_transform(wl, robin)
g1 <- ggplot() +
  geom_sf(data = bg, fill = "grey93", color = "grey80", linewidth = .1) +
  geom_sf(data = wl, aes(fill = lifeExp), color = "white", linewidth = .05) +
  facet_wrap(~year) + scale_fill_viridis_c(option = "viridis") +
  labs(title = "Life expectancy across the world, 1952-2007",
       subtitle = "Real gapminder data on real country geometry; the distribution converges upward",
       fill = "years") +
  theme(panel.grid = element_blank(), axis.text = element_blank())
ggsave(file.path(OUT, "world1_lifeexp.png"), g1, width = 11, height = 6, dpi = 120)

# ---- FIG 2: LISA Moran-quadrant world maps ------------------
qlab <- c("HH", "LH", "LL", "HL")
ql <- do.call(rbind, lapply(c(1952, 2007), function(yr) {
  g <- geom; g$year <- yr
  g$quad <- qlab[lm$quadrants[, which(yrs == yr)]]; g
}))
ql <- st_transform(ql, robin)
g2 <- ggplot() +
  geom_sf(data = bg, fill = "grey93", color = "grey80", linewidth = .1) +
  geom_sf(data = ql, aes(fill = factor(quad, levels = qlab)),
          color = "white", linewidth = .05) +
  facet_wrap(~year, ncol = 1) +
  scale_fill_manual(values = c(HH="#d7191c", LH="#abd9e9", LL="#2c7bb6", HL="#fdae61")) +
  labs(title = "Life-expectancy spatial clusters (LISA quadrants)",
       subtitle = "HH = high country in high neighbourhood; LL = low in low; via lisa_markov()",
       fill = "quadrant") +
  theme(panel.grid = element_blank(), axis.text = element_blank())
ggsave(file.path(OUT, "world2_lisa.png"), g2, width = 8.5, height = 8, dpi = 120)

# ---- FIG 3: transition matrix heatmap -----------------------
g3 <- ggplot(heat(m$matrix),
             aes(to, factor(from, levels = rev(sort(unique(from)))), fill = v)) +
  geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", v))) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
  labs(title = "Life-expectancy class transitions (5-year steps)",
       x = "to class", y = "from class", fill = "prob")
ggsave(file.path(OUT, "world3_transition.png"), g3, width = 6.6, height = 5.4, dpi = 120)

# ---- FIG 4: spatial Markov conditional matrices -------------
mats <- c(list(Pooled = sm$pooled), sm$matrices)
smdf <- do.call(rbind, Map(function(nm, mt) cbind(panel = nm, heat(mt)),
                           names(mats), mats))
smdf$panel <- factor(smdf$panel, levels = c("Pooled", names(sm$matrices)))
g4 <- ggplot(smdf, aes(to, factor(from, levels = rev(sort(unique(from)))), fill = v)) +
  geom_tile(color = "white") + facet_wrap(~panel, nrow = 2) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
  labs(title = "Spatial Markov: does the neighbourhood's level matter?",
       x = "to class", y = "from class", fill = "prob") +
  theme(panel.grid = element_blank())
ggsave(file.path(OUT, "world4_spatial_markov.png"), g4, width = 11, height = 6, dpi = 120)

cat("Saved world figures to", normalizePath(OUT), "\n")
