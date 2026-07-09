#' Markov mobility indices
#'
#' @description
#' Scalar summaries of how much mixing a transition matrix implies — larger
#' values mean more mobility (a faster-churning distribution), values near zero
#' mean a near-immobile system. Several classical indices are provided:
#'
#' \describe{
#'   \item{`prais`}{Prais index \eqn{(k - \mathrm{tr}(P)) / (k - 1)} — average
#'     probability of leaving a class.}
#'   \item{`determinant`}{\eqn{1 - |\det(P)|}.}
#'   \item{`L2`}{\eqn{1 - |\lambda_2|}, one minus the second-largest eigenvalue
#'     modulus (the mixing rate).}
#'   \item{`shorrock1`}{\eqn{(k - k \sum_i \pi_i p_{ii}) / (k - 1)}, weighted by
#'     an initial distribution `ini`.}
#'   \item{`shorrock2`}{\eqn{\sum_i \sum_j \pi_i p_{ij} |i - j| / (k - 1)},
#'     penalising longer jumps.}
#' }
#'
#' @param x A row-stochastic matrix, or an object of class `sddr_markov`.
#' @param measure One of `"prais"`, `"determinant"`, `"L2"`, `"shorrock1"`,
#'   `"shorrock2"`, or `"all"` (default) to return every index.
#' @param ini Optional initial distribution (length `k`) for the Shorrock
#'   indices; defaults to uniform.
#'
#' @return A named numeric vector of the requested index or indices.
#'
#' @references
#' Prais, S. J. (1955). Measuring social mobility. *JRSS A*, 118(1), 56-66.
#' Shorrocks, A. F. (1978). The measurement of mobility. *Econometrica*,
#' 46(5), 1013-1024. \doi{10.2307/1911433}
#'
#' @examples
#' P <- matrix(c(0.7, 0.2, 0.1,
#'               0.2, 0.6, 0.2,
#'               0.1, 0.3, 0.6), nrow = 3, byrow = TRUE)
#' mobility(P)
#' mobility(P, "prais")
#'
#' @export
mobility <- function(x, measure = "all", ini = NULL) {
  P <- if (inherits(x, "sddr_markov")) x$matrix else as.matrix(x)
  k <- nrow(P)
  if (k != ncol(P)) stop("transition matrix must be square.", call. = FALSE)
  measure <- match.arg(
    tolower(measure),
    c("all", "prais", "determinant", "l2", "shorrock1", "shorrock2")
  )
  w <- if (is.null(ini)) rep(1 / k, k) else as.numeric(ini)

  prais <- function() (k - sum(diag(P))) / (k - 1)
  determ <- function() 1 - abs(det(P))
  l2 <- function() {
    ev <- sort(abs(eigen(P, only.values = TRUE)$values))
    1 - ev[length(ev) - 1L]
  }
  b1 <- function() (k - k * sum(w * diag(P))) / (k - 1)
  b2 <- function() {
    idx <- seq_len(k)
    D <- abs(outer(idx, idx, "-"))
    sum((w * P) * D) / (k - 1)
  }

  switch(
    measure,
    prais       = c(prais = prais()),
    determinant = c(determinant = determ()),
    l2          = c(L2 = l2()),
    shorrock1   = c(shorrock1 = b1()),
    shorrock2   = c(shorrock2 = b2()),
    all = c(prais = prais(), determinant = determ(), L2 = l2(),
            shorrock1 = b1(), shorrock2 = b2())
  )
}
