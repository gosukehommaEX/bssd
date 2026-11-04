#' Chebyshev-Lobatto points
#'
#' @param n Number of points (at least 2).
#'
#' @return The points \code{-cos(pi * j / (n - 1))}, \code{j = 0, ..., n - 1},
#'   on \code{[-1, 1]} in increasing order.
#' @keywords internal
#' @noRd
cheb_nodes <- function(n) {
  -cos(pi * (0:(n - 1L)) / (n - 1L))
}
