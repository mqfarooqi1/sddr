#' Mean first passage times
#'
#' @description
#' The mean first passage time from class \eqn{i} to class \eqn{j} is the
#' expected number of periods for a unit currently in \eqn{i} to first reach
#' \eqn{j}. Together with the transition matrix these times summarise how
#' quickly a distribution mixes. The diagonal holds the *mean recurrence time*
#' \eqn{1/\pi_j} (the expected return time to a class), where \eqn{\pi} is the
#' ergodic distribution.
#'
#' Computed with the Kemeny-Snell fundamental matrix
#' \eqn{Z = (I - P + A)^{-1}}, where \eqn{A} has every row equal to \eqn{\pi};
#' then \eqn{M_{ij} = (\delta_{ij} - Z_{ij} + Z_{jj}) / \pi_j}.
#'
#' @param x A row-stochastic matrix, or an object of class `sddr_markov`.
#' @param ... Currently unused.
#'
#' @return A numeric matrix of mean first passage times.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*. Springer.
#'
#' @examples
#' P <- matrix(c(0.8, 0.2,
#'               0.3, 0.7), nrow = 2, byrow = TRUE)
#' mfpt(P)
#'
#' @seealso [steady_state()], [sojourn_time()]
#' @export
mfpt <- function(x, ...) {
  UseMethod("mfpt")
}

#' @rdname mfpt
#' @export
mfpt.default <- function(x, ...) .mfpt(as.matrix(x))

#' @rdname mfpt
#' @export
mfpt.sddr_markov <- function(x, ...) .mfpt(x$matrix)

.mfpt <- function(P) {
  P <- as.matrix(P)
  k <- nrow(P)
  if (k != ncol(P)) stop("`P` must be square.", call. = FALSE)
  pistat <- .steady_state(P)
  A <- matrix(pistat, k, k, byrow = TRUE)
  Z <- solve(diag(k) - P + A)
  M <- matrix(0, k, k, dimnames = dimnames(P))
  for (j in seq_len(k)) {
    for (i in seq_len(k)) {
      M[i, j] <- ((i == j) - Z[i, j] + Z[j, j]) / pistat[j]
    }
  }
  M
}

#' Sojourn times
#'
#' @description
#' The expected number of consecutive periods a unit remains in a class before
#' leaving it, \eqn{1 / (1 - p_{ii})} for each class \eqn{i}. Large values
#' indicate "sticky" classes with strong persistence.
#'
#' @param x A row-stochastic matrix, or an object of class `sddr_markov`.
#' @param ... Currently unused.
#'
#' @return A named numeric vector of sojourn times (one per class). Absorbing
#'   classes (\eqn{p_{ii} = 1}) return `Inf`.
#'
#' @examples
#' P <- matrix(c(0.8, 0.2,
#'               0.3, 0.7), nrow = 2, byrow = TRUE)
#' sojourn_time(P)
#'
#' @seealso [mfpt()], [steady_state()]
#' @export
sojourn_time <- function(x, ...) {
  UseMethod("sojourn_time")
}

#' @rdname sojourn_time
#' @export
sojourn_time.default <- function(x, ...) .sojourn(as.matrix(x))

#' @rdname sojourn_time
#' @export
sojourn_time.sddr_markov <- function(x, ...) .sojourn(x$matrix)

.sojourn <- function(P) {
  P <- as.matrix(P)
  d <- diag(P)
  st <- 1 / (1 - d)
  stats::setNames(st, rownames(P) %||% seq_len(nrow(P)))
}
