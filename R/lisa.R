#' LISA Markov chain
#'
#' @description
#' A Markov chain over the four quadrants of the Moran scatterplot (Rey, 2001).
#' Each unit, in each period, is placed in a quadrant according to the sign of
#' its value (relative to the cross-sectional mean) and the sign of its spatial
#' lag:
#' \describe{
#'   \item{`HH`}{high value, high neighbourhood (quadrant 1).}
#'   \item{`LH`}{low value, high neighbourhood (quadrant 2).}
#'   \item{`LL`}{low value, low neighbourhood (quadrant 3).}
#'   \item{`HL`}{high value, low neighbourhood (quadrant 4).}
#' }
#' The resulting 4 by 4 transition matrix describes how units move between local
#' spatial-association types over time — capturing the co-evolution of a unit
#' and its neighbourhood.
#'
#' @inheritParams full_rank_markov
#' @param weights A square numeric spatial weights matrix whose row and column
#'   names match the unit ids.
#' @param row_standardize Logical; row-standardise `weights` before computing
#'   the spatial lag (default `TRUE`).
#'
#' @return An object of class `sddr_markov`: a 4 by 4 transition matrix over the
#'   quadrants `HH`, `LH`, `LL`, `HL`, plus the per-unit-per-period `quadrants`
#'   matrix.
#'
#' @references
#' Rey, S. J. (2001). Spatial empirics for economic growth and convergence.
#' *Geographical Analysis*, 33(3), 195-214.
#' \doi{10.1111/j.1538-4632.2001.tb00444.x}
#'
#' @examples
#' set.seed(1)
#' n <- 9; periods <- 8
#' W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
#' for (i in 1:n) {                       # ring contiguity
#'   W[i, i %% n + 1] <- 1
#'   W[i, (i - 2) %% n + 1] <- 1
#' }
#' df <- data.frame(id = rep(1:n, times = periods),
#'                  time = rep(1:periods, each = n),
#'                  value = rnorm(n * periods))
#' lisa_markov(df, "id", "time", "value", weights = W)
#'
#' @seealso [spatial_markov()], [markov()]
#' @export
lisa_markov <- function(data, id, time, value, weights,
                        row_standardize = TRUE) {
  Y <- .pivot_wide(data, id, time, value)
  ids <- rownames(Y)

  W <- as.matrix(weights)
  if (is.null(rownames(W)) || is.null(colnames(W))) {
    stop("`weights` must have row and column names matching unit ids.",
         call. = FALSE)
  }
  if (!all(ids %in% rownames(W))) {
    stop("`weights` is missing some unit ids present in the data.",
         call. = FALSE)
  }
  W <- W[ids, ids, drop = FALSE]
  if (row_standardize) {
    rs <- rowSums(W); rs[rs == 0] <- 1; W <- W / rs
  }

  n <- nrow(Y); periods <- ncol(Y)
  q <- matrix(0L, n, periods)
  for (t in seq_len(periods)) {
    z <- Y[, t] - mean(Y[, t])          # centre; only signs matter
    lag <- as.numeric(W %*% z)
    zp <- z > 0; lp <- lag > 0
    q[, t] <- 1L * (zp & lp) + 2L * (!zp & lp) +
              3L * (!zp & !lp) + 4L * (zp & !lp)
  }

  m <- .rank_markov(q, 4L, kind = "LISA")
  quad <- c("HH", "LH", "LL", "HL")
  dimnames(m$matrix) <- dimnames(m$transitions) <- list(quad, quad)
  names(m$steady_state) <- quad
  m$quadrants <- q
  m
}
