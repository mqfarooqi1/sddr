#' Ergodic (steady-state) distribution of a Markov chain
#'
#' @description
#' Computes the long-run stationary distribution of a row-stochastic transition
#' matrix, i.e. the probability vector \eqn{\pi} satisfying
#' \eqn{\pi P = \pi} with \eqn{\sum_i \pi_i = 1}. It is obtained as the left
#' eigenvector of `P` associated with the eigenvalue closest to 1.
#'
#' @param x A row-stochastic matrix, or an object of class `sddr_markov`.
#' @param ... Currently unused.
#'
#' @return A named numeric vector giving the stationary distribution.
#'
#' @examples
#' P <- matrix(c(0.8, 0.2,
#'               0.3, 0.7), nrow = 2, byrow = TRUE)
#' steady_state(P)
#'
#' @export
steady_state <- function(x, ...) {
  UseMethod("steady_state")
}

#' @rdname steady_state
#' @export
steady_state.default <- function(x, ...) {
  .steady_state(as.matrix(x))
}

#' @rdname steady_state
#' @export
steady_state.sddr_markov <- function(x, ...) {
  x$steady_state
}

# Left eigenvector of P for eigenvalue ~1, normalised to sum to 1.
.steady_state <- function(P) {
  P <- as.matrix(P)
  k <- nrow(P)
  if (k != ncol(P)) stop("`P` must be square.", call. = FALSE)
  ev <- eigen(t(P))
  i <- which.min(abs(ev$values - 1))
  v <- Re(ev$vectors[, i])
  v <- v / sum(v)
  stats::setNames(v, rownames(P) %||% seq_len(k))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
