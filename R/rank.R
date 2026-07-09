#' Kendall's tau rank correlation
#'
#' @description
#' Kendall's \eqn{\tau} measures the rank concordance between two variables —
#' for distribution dynamics, typically a set of units' values at two points in
#' time. It is the (concordant minus discordant) pair count, normalised, and is
#' a natural summary of *positional* (exchange) mobility: \eqn{\tau = 1} means
#' ranks are perfectly preserved, \eqn{\tau = -1} a complete reversal, and
#' \eqn{\tau \approx 0} extensive rank shuffling. Ties are handled with the
#' tau-b correction.
#'
#' @param x,y Numeric vectors of equal length.
#'
#' @return An object of class `sddr_tau`: a list with the statistic `tau`, the
#'   `concordant` and `discordant` pair counts, an asymptotic (normal) two-sided
#'   `p_value`, and `n`.
#'
#' @references
#' Kendall, M. G. (1938). A new measure of rank correlation. *Biometrika*,
#' 30(1/2), 81-93. \doi{10.2307/2332226}
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(30)
#' y <- x + rnorm(30)          # correlated ranks
#' tau(x, y)
#'
#' @seealso [theta()]
#' @export
tau <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }
  n <- length(x)
  if (n < 2L) stop("need at least two observations.", call. = FALSE)

  s <- sign(outer(x, x, "-")) * sign(outer(y, y, "-"))
  upper <- upper.tri(s)
  concordant <- sum(s[upper] > 0)
  discordant <- sum(s[upper] < 0)

  n0 <- n * (n - 1) / 2
  tie <- function(v) sum(vapply(table(v), function(t) t * (t - 1) / 2, numeric(1)))
  tx <- tie(x); ty <- tie(y)
  tau_b <- (concordant - discordant) / sqrt((n0 - tx) * (n0 - ty))

  var0 <- n * (n - 1) * (2 * n + 5) / 18
  z <- (concordant - discordant) / sqrt(var0)
  p <- 2 * stats::pnorm(-abs(z))

  structure(
    list(tau = tau_b, concordant = concordant, discordant = discordant,
         p_value = p, n = n),
    class = "sddr_tau"
  )
}

#' @export
print.sddr_tau <- function(x, digits = 4, ...) {
  cat("<sddr> Kendall's tau\n")
  cat(sprintf("  tau = %.*f  (concordant %d, discordant %d, n = %d)\n",
              digits, x$tau, x$concordant, x$discordant, x$n))
  cat(sprintf("  asymptotic two-sided p-value = %.4g\n", x$p_value))
  invisible(x)
}

#' Theta rank-mobility statistic
#'
#' @description
#' The Theta statistic (Rey, 2004) measures how much of the overall rank
#' movement in a system is *within-regime* exchange. For each period-to-period
#' transition it is the sum over regimes of the absolute *net* within-regime
#' rank change, divided by the total absolute rank change. Values near zero
#' indicate that rank movement is dominated by within-regime shuffling that
#' nets out; larger values indicate coherent regime-level movement.
#'
#' @param y A numeric matrix of values, units in rows and time periods in
#'   columns (at least two columns).
#' @param regime A vector of length `nrow(y)` assigning each unit to a regime.
#'
#' @return An object of class `sddr_theta`: a list with `theta` (one value per
#'   period-to-period transition), the per-period `total` absolute rank change,
#'   and the column rank matrix `ranks`.
#'
#' @references
#' Rey, S. J. (2004). Spatial dependence in the evolution of regional income
#' distributions. In *Spatial Econometrics and Spatial Statistics*, 194-213.
#'
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(40), nrow = 10)    # 10 units, 4 periods
#' regime <- rep(c("a", "b"), each = 5)
#' theta(y, regime)
#'
#' @seealso [tau()]
#' @export
theta <- function(y, regime) {
  y <- as.matrix(y)
  if (ncol(y) < 2L) stop("`y` needs at least two time periods (columns).",
                         call. = FALSE)
  if (length(regime) != nrow(y)) {
    stop("`regime` must have one entry per row of `y`.", call. = FALSE)
  }

  ranks <- apply(y, 2, rank)                     # rank within each period
  tt <- ncol(ranks)
  rd <- ranks[, -1, drop = FALSE] - ranks[, -tt, drop = FALSE]
  total <- colSums(abs(rd))

  regs <- sort(unique(regime))
  within <- vapply(regs,
                   function(r) abs(colSums(rd[regime == r, , drop = FALSE])),
                   numeric(tt - 1))
  within <- matrix(within, nrow = tt - 1)
  th <- rowSums(within) / total

  structure(
    list(theta = th, total = total, ranks = ranks, regimes = regs),
    class = "sddr_theta"
  )
}

#' @export
print.sddr_theta <- function(x, digits = 4, ...) {
  cat("<sddr> Theta rank-mobility statistic\n")
  cat(sprintf("  %d regimes | %d period transitions\n",
              length(x$regimes), length(x$theta)))
  cat("  theta by transition:\n")
  print(round(x$theta, digits))
  invisible(x)
}
