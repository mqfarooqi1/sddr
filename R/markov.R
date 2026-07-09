#' Classic Markov chain for distribution dynamics
#'
#' @description
#' Discretises a continuous longitudinal variable into ordered classes and
#' estimates a discrete-time, first-order Markov transition probability matrix
#' together with its ergodic (steady-state) distribution. This is the classic
#' distribution-dynamics estimator of Quah (1993): the transition matrix
#' summarises how units move between parts of the distribution from one period
#' to the next, and the steady-state distribution describes where the process
#' settles in the long run.
#'
#' Input is *tidy*, long-format panel data with one row per unit-period. The
#' estimator is agnostic to what the units are (regions, firms, households,
#' \dots) and to whether space is involved; spatial conditioning is added by a
#' companion `spatial_markov()` function (available in a later release).
#'
#' @param data A data frame in long format with one row per unit and period.
#' @param id Name (character scalar) of the column identifying the unit.
#' @param time Name (character scalar) of the column giving the period. Periods
#'   are compared in sorted order within each unit.
#' @param value Name (character scalar) of the numeric column to analyse.
#' @param k Integer number of classes to discretise `value` into. Default `5`
#'   (quintiles). Ignored if `breaks` is supplied.
#' @param breaks Optional numeric vector of *interior* class boundaries. If
#'   supplied these fixed cut points are used instead of data-driven quantiles,
#'   giving `length(breaks) + 1` classes.
#' @param fixed Logical. If `TRUE` (default) a single set of quantile breaks is
#'   estimated from the pooled distribution and applied to every period
#'   (absolute mobility). If `FALSE`, breaks are re-estimated within each period
#'   (relative mobility). Ignored when `breaks` is supplied.
#'
#' @return An object of class `sddr_markov`: a list with elements
#'   \describe{
#'     \item{`matrix`}{The `k` by `k` row-stochastic transition probability
#'       matrix `P`.}
#'     \item{`transitions`}{The `k` by `k` matrix of observed transition counts.}
#'     \item{`steady_state`}{The ergodic (stationary) distribution of `P`.}
#'     \item{`classes`}{Number of classes `k`.}
#'     \item{`breaks`}{The class boundaries used (when `fixed`).}
#'     \item{`n`}{Number of units contributing at least one transition.}
#'     \item{`n_transitions`}{Total number of observed transitions.}
#'     \item{`data`}{The input data augmented with the integer `state` column.}
#'   }
#'
#' @references
#' Quah, D. (1993). Empirical cross-section dynamics in economic growth.
#' *European Economic Review*, 37(2-3), 426-434.
#' \doi{10.1016/0014-2921(93)90031-3}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   id = rep(1:50, each = 4),
#'   time = rep(2000:2003, times = 50),
#'   value = rnorm(200)
#' )
#' m <- markov(df, id = "id", time = "time", value = "value", k = 4)
#' m
#' m$steady_state
#'
#' @seealso [steady_state()]
#' @export
markov <- function(data, id, time, value, k = 5, breaks = NULL, fixed = TRUE) {
  stopifnot(is.data.frame(data))
  for (col in c(id, time, value)) {
    if (!is.character(col) || length(col) != 1L) {
      stop("`id`, `time` and `value` must each be a single column name.",
           call. = FALSE)
    }
    if (!col %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", col), call. = FALSE)
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

  cls <- .classify(d$value, d$time, k = k, breaks = breaks, fixed = fixed)
  d$state <- as.integer(cls)
  k <- attr(cls, "k")
  states <- seq_len(k)

  # Order within unit, then count consecutive-period transitions per unit.
  d <- d[order(d$id, d$time), , drop = FALSE]
  counts <- matrix(0L, k, k, dimnames = list(states, states))
  n_units <- 0L
  for (rows in split(seq_len(nrow(d)), d$id)) {
    s <- d$state[rows]
    if (length(s) >= 2L) {
      n_units <- n_units + 1L
      tab <- table(
        factor(s[-length(s)], levels = states),
        factor(s[-1L],         levels = states)
      )
      counts <- counts + as.matrix(tab)
    }
  }

  rs <- rowSums(counts)
  P <- counts / ifelse(rs == 0, 1, rs)
  P[rs == 0, ] <- 0

  structure(
    list(
      matrix        = P,
      transitions   = counts,
      steady_state  = .steady_state(P),
      classes       = k,
      breaks        = attr(cls, "breaks"),
      fixed         = fixed,
      n             = n_units,
      n_transitions = sum(counts),
      data          = d,
      call          = match.call()
    ),
    class = "sddr_markov"
  )
}

# Discretise a numeric vector into ordered integer classes.
# Returns an integer vector with attributes `k` and `breaks`.
.classify <- function(value, time, k = 5, breaks = NULL, fixed = TRUE) {
  if (!is.null(breaks)) {
    cuts <- c(-Inf, sort(unique(as.numeric(breaks))), Inf)
    state <- as.integer(cut(value, breaks = cuts, labels = FALSE,
                            include.lowest = TRUE))
    return(structure(state, k = length(cuts) - 1L, breaks = cuts))
  }

  k <- as.integer(k)
  if (k < 2L) stop("`k` must be at least 2.", call. = FALSE)
  probs <- seq(0, 1, length.out = k + 1L)

  cut_fixed <- function(v) {
    qs <- stats::quantile(v, probs = probs, names = FALSE, type = 7)
    qs[1L] <- -Inf
    qs[length(qs)] <- Inf
    if (anyDuplicated(qs)) {
      stop("Quantile breaks are not unique (too many ties in `value` for ",
           k, " classes); reduce `k` or supply `breaks`.", call. = FALSE)
    }
    qs
  }

  if (fixed) {
    qs <- cut_fixed(value)
    state <- as.integer(cut(value, breaks = qs, labels = FALSE,
                            include.lowest = TRUE))
    structure(state, k = k, breaks = qs)
  } else {
    state <- integer(length(value))
    for (tt in unique(time)) {
      idx <- which(time == tt)
      qs <- cut_fixed(value[idx])
      state[idx] <- as.integer(cut(value[idx], breaks = qs, labels = FALSE,
                                   include.lowest = TRUE))
    }
    structure(state, k = k, breaks = NA_real_)
  }
}
