#' @export
print.sddr_markov <- function(x, digits = 3, ...) {
  cat("<sddr> classic Markov chain\n")
  cat(sprintf("  units: %d | transitions: %d | classes: %d | breaks: %s\n",
              x$n, x$n_transitions, x$classes,
              if (isTRUE(x$fixed)) "fixed" else "per-period"))
  cat("\nTransition probability matrix (rows = from, cols = to):\n")
  print(round(x$matrix, digits))
  cat("\nErgodic (steady-state) distribution:\n")
  print(round(x$steady_state, digits))
  invisible(x)
}
