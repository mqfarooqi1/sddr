#' Spatial Markov chain for distribution dynamics
#'
#' @description
#' Estimates a spatial Markov chain (Rey, 2001): a set of class-transition
#' probability matrices, each conditioned on the *spatial context* of a unit at
#' the start of the period. Where [markov()] asks "given a unit's own class,
#' where does it move next?", `spatial_markov()` asks "given a unit's own class
#' **and the class of its spatial neighbourhood**, where does it move next?".
#' Comparing the conditional matrices with the pooled (aspatial) matrix reveals
#' whether, and how, spatial context shapes distribution dynamics.
#'
#' The neighbourhood is summarised by the *spatial lag*: a (row-standardised)
#' weighted average of neighbours' values. Supply either a spatial weights
#' matrix via `weights` (the lag is computed for you, per period), or a
#' precomputed spatial-lag column via `lag` (keeping the estimator independent
#' of any particular spatial-data stack).
#'
#' @inheritParams markov
#' @param weights Optional spatial weights: a square numeric matrix whose row
#'   and column names match the unit ids, or an \pkg{spdep} `listw` or `nb`
#'   object (e.g. built from an `sf` layer with [spdep::poly2nb()] and
#'   [spdep::nb2listw()]). Used to compute the spatial lag of `value` within
#'   each period. Requires a balanced panel. Exactly one of `weights` or `lag`
#'   must be supplied.
#' @param lag Optional name (character scalar) of a column giving a precomputed
#'   spatial lag of `value`.
#' @param m Integer number of spatial-lag classes to condition on. Defaults to
#'   `k`.
#' @param lag_breaks Optional numeric vector of *interior* class boundaries for
#'   the spatial lag (the lag analogue of `breaks`). If supplied, these fixed
#'   cut points are used instead of data-driven quantiles. Supplying both
#'   `breaks` and `lag_breaks` makes the classification fully explicit, which is
#'   how exact agreement with `giddy` (`cutoffs` / `lag_cutoffs`) is obtained.
#' @param row_standardize Logical; row-standardise `weights` before computing
#'   the lag (default `TRUE`). Ignored when `lag` is supplied.
#'
#' @return An object of class `sddr_spatial_markov`: a list with elements
#'   \describe{
#'     \item{`matrices`}{Length-`m` list of `k` by `k` conditional transition
#'       probability matrices, one per spatial-lag class.}
#'     \item{`transitions`}{Length-`m` list of conditional count matrices.}
#'     \item{`pooled`}{The `k` by `k` pooled (aspatial) transition matrix; this
#'       equals the matrix from [markov()] with the same class definition.}
#'     \item{`pooled_transitions`}{Pooled count matrix.}
#'     \item{`steady_states`}{Length-`m` list of conditional ergodic
#'       distributions (`NA` where a conditional matrix has unvisited classes).}
#'     \item{`classes`, `lag_classes`}{`k` and `m`.}
#'     \item{`n`, `n_transitions`}{Units contributing, and total transitions.}
#'     \item{`data`}{Input augmented with `lagval`, value state and lag state.}
#'   }
#'
#' @references
#' Rey, S. J. (2001). Spatial empirics for economic growth and convergence.
#' *Geographical Analysis*, 33(3), 195-214.
#' \doi{10.1111/j.1538-4632.2001.tb00444.x}
#'
#' @examples
#' set.seed(1)
#' n <- 12; periods <- 8
#' # Ring contiguity: each unit neighbours the next and previous on a circle.
#' W <- matrix(0, n, n, dimnames = list(1:n, 1:n))
#' for (i in 1:n) {
#'   W[i, i %% n + 1] <- 1
#'   W[i, (i - 2) %% n + 1] <- 1
#' }
#' df <- do.call(rbind, lapply(seq_len(periods), function(p)
#'   data.frame(id = 1:n, time = p, value = rnorm(n))))
#' sm <- spatial_markov(df, "id", "time", "value", weights = W, k = 3)
#' sm
#'
#' @seealso [markov()]
#' @export
spatial_markov <- function(data, id, time, value,
                           weights = NULL, lag = NULL,
                           k = 5, m = k, breaks = NULL, lag_breaks = NULL,
                           fixed = TRUE, row_standardize = TRUE) {
  stopifnot(is.data.frame(data))
  if (is.null(weights) == is.null(lag)) {
    stop("Supply exactly one of `weights` or `lag`.", call. = FALSE)
  }
  for (col in c(id, time, value)) {
    if (!is.character(col) || length(col) != 1L || !col %in% names(data)) {
      stop("`id`, `time` and `value` must each name a column in `data`.",
           call. = FALSE)
    }
  }

  d <- data.frame(
    id    = data[[id]],
    time  = data[[time]],
    value = as.numeric(data[[value]]),
    stringsAsFactors = FALSE
  )
  if (anyNA(d$value)) {
    stop("`value` contains missing values; remove or impute them first.",
         call. = FALSE)
  }

  # Spatial lag of `value`.
  if (!is.null(lag)) {
    if (!lag %in% names(data)) {
      stop(sprintf("Lag column '%s' not found in `data`.", lag), call. = FALSE)
    }
    d$lagval <- as.numeric(data[[lag]])
    if (anyNA(d$lagval)) {
      stop("`lag` contains missing values.", call. = FALSE)
    }
  } else {
    d$lagval <- .spatial_lag(d, weights, row_standardize)
  }

  # Classify value (k classes) and spatial lag (m classes).
  vcls <- .classify(d$value, d$time, k = k, breaks = breaks, fixed = fixed)
  lcls <- .classify(d$lagval, d$time, k = m, breaks = lag_breaks, fixed = fixed)
  d$vstate <- as.integer(vcls); k <- attr(vcls, "k")
  d$lstate <- as.integer(lcls); m <- attr(lcls, "k")
  sv <- seq_len(k)

  # Count transitions, conditioning on the origin-period spatial-lag class.
  d <- d[order(d$id, d$time), , drop = FALSE]
  cond <- replicate(m, matrix(0L, k, k, dimnames = list(sv, sv)),
                    simplify = FALSE)
  n_units <- 0L
  for (rows in split(seq_len(nrow(d)), d$id)) {
    vs <- d$vstate[rows]; ls <- d$lstate[rows]
    nobs <- length(rows)
    if (nobs >= 2L) {
      n_units <- n_units + 1L
      for (i in seq_len(nobs - 1L)) {
        cond[[ls[i]]][vs[i], vs[i + 1L]] <- cond[[ls[i]]][vs[i], vs[i + 1L]] + 1L
      }
    }
  }
  pooled_counts <- Reduce(`+`, cond)

  norm <- function(C) {
    rs <- rowSums(C)
    P <- C / ifelse(rs == 0, 1, rs)
    P[rs == 0, ] <- 0
    P
  }
  cond_steady <- function(P) {
    if (any(rowSums(P) == 0)) {
      stats::setNames(rep(NA_real_, nrow(P)), rownames(P))
    } else {
      .steady_state(P)
    }
  }

  matrices <- lapply(cond, norm)
  nm <- paste("lag", seq_len(m))
  names(matrices) <- names(cond) <- nm

  structure(
    list(
      matrices           = matrices,
      transitions        = cond,
      pooled             = norm(pooled_counts),
      pooled_transitions = pooled_counts,
      steady_states      = stats::setNames(lapply(matrices, cond_steady), nm),
      classes            = k,
      lag_classes        = m,
      fixed              = fixed,
      n                  = n_units,
      n_transitions      = sum(pooled_counts),
      data               = d,
      call               = match.call()
    ),
    class = "sddr_spatial_markov"
  )
}

# Row-standardised spatial lag of `value`, computed within each period.
.spatial_lag <- function(d, W, row_standardize) {
  W <- .as_weights_matrix(W)
  if (nrow(W) != ncol(W)) stop("`weights` must be a square matrix.",
                               call. = FALSE)
  ids <- rownames(W)
  if (is.null(ids) || is.null(colnames(W))) {
    stop("`weights` must have row and column names matching unit ids.",
         call. = FALSE)
  }
  if (row_standardize) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }
  lagval <- numeric(nrow(d))
  for (tt in unique(d$time)) {
    idx <- which(d$time == tt)
    yt <- d$value[idx]
    names(yt) <- as.character(d$id[idx])
    if (!all(ids %in% names(yt))) {
      stop("A balanced panel is required with `weights`: period '", tt,
           "' is missing some units present in `weights`.", call. = FALSE)
    }
    lag_tt <- as.numeric(W %*% yt[ids])
    names(lag_tt) <- ids
    lagval[idx] <- lag_tt[as.character(d$id[idx])]
  }
  lagval
}
