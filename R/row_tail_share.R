#' Share of the highest Chebyshev coefficients in the density of each row
#'
#' For each row of a matrix of log densities at Chebyshev points, computes
#' the Chebyshev coefficients of the density (scaled to a maximum of 1) and
#' returns the share of the absolute coefficients of the highest quarter of
#' the degrees. The share is small when the interpolation of the row is
#' accurate; it is used to choose the coordinate of each row.
#'
#' @param L Matrix of log densities, one row per node of the first
#'   coordinate and one column per Chebyshev point.
#'
#' @return Vector with the share of each row (\code{Inf} for rows whose
#'   density is zero at all points).
#' @keywords internal
#' @noRd
row_tail_share <- function(L) {
  L[is.na(L)] <- -Inf
  m <- apply(L, 1L, max)
  f <- exp(L - m)
  f[!is.finite(m), ] <- 0
  a <- abs(f %*% t(cheb_coef_matrix(ncol(L))))
  k0 <- floor(0.75 * ncol(L))
  out <- rowSums(a[, (k0 + 1L):ncol(L), drop = FALSE]) /
    pmax(rowSums(a), .Machine$double.xmin)
  out[!is.finite(m)] <- Inf
  out
}
