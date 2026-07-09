#' Full-rank Markov chain
#'
#' @description
#' A Markov chain over *ranks* rather than discretised classes: each unit's
#' rank within the cross-section (1 = highest value, `n` = lowest) becomes its
#' state, and the `n` by `n` matrix records how ranks transition from one
#' period to the next (Rey, 2014). Unlike [markov()] no binning is needed — the
#' full rank ordering is retained.
#'
#' @param data A data frame in long format with one row per unit and period.
#'   The panel must be balanced (every unit present in every period).
#' @param id,time,value Character scalars naming the unit, period, and numeric
#'   value columns.
#'
#' @return An object of class `sddr_markov` (an `n` by `n` rank transition
#'   matrix, its counts, and ergodic distribution).
#'
#' @references
#' Rey, S. J. (2014). Fast algorithms for a space-time concordance measure.
#' *Computational Statistics*, 29(3-4), 799-811.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   id = rep(1:6, times = 8),
#'   time = rep(1:8, each = 6),
#'   value = rnorm(48)
#' )
#' full_rank_markov(df, "id", "time", "value")
#'
#' @seealso [markov()], [geo_rank_markov()]
#' @export
full_rank_markov <- function(data, id, time, value) {
  W <- .pivot_wide(data, id, time, value)
  n <- nrow(W)
  ranks_asc <- apply(W, 2, function(col) rank(col, ties.method = "first"))
  ranks_hl <- n - ranks_asc + 1                     # 1 = highest value
  .rank_markov(ranks_hl, n, kind = "full-rank")
}

#' Geographic-rank Markov chain
#'
#' @description
#' A Markov chain that tracks, for each rank position, *which unit* occupies it
#' over time (Rey, 2016). States are the geographic units, so the `n` by `n`
#' matrix records how units exchange rank positions — a spatially explicit view
#' of positional mobility.
#'
#' @inheritParams full_rank_markov
#'
#' @return An object of class `sddr_markov` (an `n` by `n` transition matrix
#'   over units, its counts, and ergodic distribution).
#'
#' @references
#' Rey, S. J. (2016). Space-time patterns of rank concordance: Local indicators
#' of mobility association. *Annals of the AAG*, 106(4), 788-803.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   id = rep(1:6, times = 8),
#'   time = rep(1:8, each = 6),
#'   value = rnorm(48)
#' )
#' geo_rank_markov(df, "id", "time", "value")
#'
#' @seealso [full_rank_markov()], [markov()]
#' @export
geo_rank_markov <- function(data, id, time, value) {
  W <- .pivot_wide(data, id, time, value)
  n <- nrow(W)
  ranks_asc <- apply(W, 2, function(col) rank(col, ties.method = "first"))
  geo <- apply(ranks_asc, 2, order)                 # unit at each rank position
  .rank_markov(geo, n, kind = "geo-rank")
}

# Reshape tidy long data to a units x periods matrix; require a balanced panel.
.pivot_wide <- function(data, id, time, value) {
  for (col in c(id, time, value)) {
    if (!is.character(col) || length(col) != 1L || !col %in% names(data)) {
      stop("`id`, `time` and `value` must each name a column in `data`.",
           call. = FALSE)
    }
  }
  ids <- sort(unique(data[[id]]))
  times <- sort(unique(data[[time]]))
  W <- matrix(NA_real_, length(ids), length(times),
              dimnames = list(as.character(ids), as.character(times)))
  W[cbind(match(data[[id]], ids), match(data[[time]], times))] <-
    as.numeric(data[[value]])
  if (anyNA(W)) {
    stop("a balanced panel is required (every unit present in every period).",
         call. = FALSE)
  }
  W
}

# Classic Markov counting over an integer state matrix (rows = series over
# time, columns = periods, entries = states in 1..k).
.rank_markov <- function(S, k, kind) {
  counts <- matrix(0L, k, k, dimnames = list(seq_len(k), seq_len(k)))
  for (r in seq_len(nrow(S))) {
    s <- S[r, ]
    from <- s[-length(s)]; to <- s[-1L]
    for (i in seq_along(from)) {
      counts[from[i], to[i]] <- counts[from[i], to[i]] + 1L
    }
  }
  rs <- rowSums(counts)
  P <- counts / ifelse(rs == 0, 1, rs)
  P[rs == 0, ] <- 0
  structure(
    list(matrix = P, transitions = counts, steady_state = .steady_state(P),
         classes = k, n = nrow(S), n_transitions = sum(counts),
         kind = kind, fixed = TRUE),
    class = "sddr_markov"
  )
}
