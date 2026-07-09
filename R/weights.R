# Coerce a spatial-weights specification to a numeric matrix with id dimnames.
#
# Accepts a plain matrix, an 'spdep' listw object, or an 'nb' neighbours
# object, so the spatial estimators can be driven directly from the weights
# built for an sf layer (e.g. spdep::poly2nb() |> spdep::nb2listw()). 'spdep'
# is only needed when a listw/nb is supplied, and stays in Suggests.
.as_weights_matrix <- function(w) {
  if (inherits(w, "listw")) {
    .need_spdep("a 'listw' weights object")
    m <- spdep::listw2mat(w)
    .set_weight_ids(m, attr(w, "region.id"))
  } else if (inherits(w, "nb")) {
    .need_spdep("an 'nb' neighbours object")
    m <- spdep::nb2mat(w, style = "B", zero.policy = TRUE)
    .set_weight_ids(m, attr(w, "region.id"))
  } else if (is.matrix(w) || is.data.frame(w)) {
    as.matrix(w)
  } else {
    stop("`weights` must be a matrix, an 'spdep' 'listw', or an 'nb' object.",
         call. = FALSE)
  }
}

.need_spdep <- function(what) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required to use ", what,
         "; install it or pass a weights matrix instead.", call. = FALSE)
  }
}

.set_weight_ids <- function(m, ids) {
  if (!is.null(ids)) dimnames(m) <- list(as.character(ids), as.character(ids))
  m
}
