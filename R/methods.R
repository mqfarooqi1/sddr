#' @export
print.sddr_markov <- function(x, digits = 3, ...) {
  cat(sprintf("<sddr> %s Markov chain\n", x$kind %||% "classic"))
  cat(sprintf("  units: %d | transitions: %d | classes: %d | breaks: %s\n",
              x$n, x$n_transitions, x$classes,
              if (isTRUE(x$fixed)) "fixed" else "per-period"))
  cat("\nTransition probability matrix (rows = from, cols = to):\n")
  print(round(x$matrix, digits))
  cat("\nErgodic (steady-state) distribution:\n")
  print(round(x$steady_state, digits))
  invisible(x)
}

#' @export
print.sddr_spatial_markov <- function(x, digits = 3, ...) {
  cat("<sddr> spatial Markov chain\n")
  cat(sprintf(paste0("  units: %d | transitions: %d | value classes: %d | ",
                     "lag classes: %d | breaks: %s\n"),
              x$n, x$n_transitions, x$classes, x$lag_classes,
              if (isTRUE(x$fixed)) "fixed" else "per-period"))
  cat("\nPooled transition matrix (aspatial):\n")
  print(round(x$pooled, digits))
  for (l in seq_along(x$matrices)) {
    cat(sprintf("\nConditional on spatial-lag class %d (%d transitions):\n",
                l, sum(x$transitions[[l]])))
    print(round(x$matrices[[l]], digits))
  }
  invisible(x)
}
