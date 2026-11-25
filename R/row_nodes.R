#' Chebyshev points of the interval of each row
#'
#' @param lo,hi Lower and upper limits, one per row.
#' @param n Number of Chebyshev points per row (at least 2).
#'
#' @return A matrix with one row per interval and \code{n} columns: the
#'   points of \code{cheb_nodes(n)} mapped to \code{[lo[i], hi[i]]}.
#' @keywords internal
#' @noRd
row_nodes <- function(lo, hi, n) {
  (hi + lo) / 2 + outer((hi - lo) / 2, cheb_nodes(n))
}
