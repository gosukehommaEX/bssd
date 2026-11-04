#' Matrix mapping function values to Chebyshev coefficients
#'
#' For values \code{f} of a function at the points returned by
#' \code{cheb_nodes(n)}, \code{cheb_coef_matrix(n) \%*\% f} gives the
#' coefficients \code{a_0, ..., a_(n - 1)} of the interpolating polynomial
#' \code{sum(a_k * T_k(x))}, where \code{T_k} is the Chebyshev polynomial of
#' the first kind.
#'
#' @param n Number of points (at least 2).
#'
#' @return An \code{n} by \code{n} matrix.
#' @keywords internal
#' @noRd
cheb_coef_matrix <- function(n) {
  j <- 0:(n - 1L)
  theta <- pi * (n - 1L - j) / (n - 1L)
  k <- 0:(n - 1L)
  C <- cos(outer(k, theta)) * 2 / (n - 1L)
  C[, 1L] <- C[, 1L] / 2
  C[, n] <- C[, n] / 2
  C[1L, ] <- C[1L, ] / 2
  C[n, ] <- C[n, ] / 2
  C
}
