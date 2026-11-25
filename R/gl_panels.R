#' Composite Gauss-Legendre rule with panels split at given points
#'
#' Without split points inside the interval, this is the Gauss-Legendre rule
#' with \code{n} nodes. Otherwise the interval is divided at the split points
#' and each panel gets a number of nodes proportional to its length, but at
#' least \code{ceiling(n / 2)}, so that a function with a kink at a split
#' point is integrated accurately also when the mass of the integrand lies in
#' a short panel.
#'
#' @param lo,hi Limits of the interval.
#' @param n Number of nodes without split points.
#' @param splits Split points; those outside \code{(lo, hi)} are ignored.
#'
#' @return A list with elements \code{x} (nodes in increasing order) and
#'   \code{w} (weights).
#' @keywords internal
#' @noRd
gl_panels <- function(lo, hi, n, splits = numeric(0)) {
  s <- sort(splits[is.finite(splits) & splits > lo & splits < hi])
  if (length(s) == 0L) {
    g <- gl_nodes(n)
    return(list(x = (hi - lo) / 2 * g$x + (hi + lo) / 2,
                w = g$w * (hi - lo) / 2))
  }
  e <- c(lo, s, hi)
  m <- pmax(as.integer(ceiling(n / 2)),
            as.integer(floor(n * diff(e) / (hi - lo) + 0.5)))
  x <- numeric(0)
  w <- numeric(0)
  for (j in seq_along(m)) {
    g <- gl_nodes(m[j])
    a <- e[j]
    b <- e[j + 1L]
    x <- c(x, (b - a) / 2 * g$x + (b + a) / 2)
    w <- c(w, g$w * (b - a) / 2)
  }
  list(x = x, w = w)
}
