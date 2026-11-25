#' Two-dimensional posterior integration with an interval for each row
#'
#' Integrates a posterior density given at Gauss-Legendre nodes \code{u} of
#' the first coordinate and, for each node (one row), at Chebyshev points of
#' an interval of the second coordinate that may differ from row to row
#' (\code{row_nodes(lo, hi, n_v)}). The density of each row is represented
#' by its Chebyshev series, so the posterior probability below a threshold
#' that depends on the row is obtained without new evaluations of the
#' density.
#'
#' @param L Matrix of log posterior densities with one row per node
#'   \code{u} and one column per Chebyshev point
#'   (\code{n_v = ncol(L)}). The coordinate of the second dimension may
#'   differ between rows, as long as \code{L} is the density with respect to
#'   that coordinate.
#' @param u Nodes of the first coordinate.
#' @param gw Quadrature weights of \code{u}.
#' @param lo,hi Interval of the second coordinate of each row.
#'
#' @return A list with the nodes \code{u}, the matrix \code{v} of the
#'   Chebyshev points of the rows, the matrix \code{wq} of quadrature weights
#'   times the normalized posterior density (its sum is 1), and the function
#'   \code{cdf_rows(s)}: the posterior probability that the second coordinate
#'   is below \code{s[i]} in the \code{i}-th row (\code{s} is recycled to
#'   the number of rows; values outside the interval of a row give 0 or the
#'   whole mass of the row).
#' @keywords internal
#' @noRd
post2d <- function(L, u, gw, lo, hi) {
  n_u <- length(u)
  n_v <- ncol(L)
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  L <- L + log(gw * half)
  L[is.na(L)] <- -Inf
  if (!any(is.finite(L))) {
    stop("The posterior density is zero on the whole integration box.",
         call. = FALSE)
  }
  W <- exp(L - max(L))
  ccw <- cc_weights(n_v)
  W <- W / sum(W %*% ccw)
  anti <- cheb_antideriv(cheb_coef_matrix(n_v) %*% t(W))
  list(
    u = u,
    v = row_nodes(lo, hi, n_v),
    wq = W * rep(ccw, each = n_u),
    cdf_rows = function(s) {
      s <- rep_len(s, n_u)
      p <- sum(cheb_cumint(anti, (s - mid) / half))
      min(max(p, 0), 1)
    }
  )
}
