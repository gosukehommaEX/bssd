#' Two-dimensional posterior integration on a product grid
#'
#' Evaluates the log posterior density on the product of Gauss-Legendre
#' nodes in the first coordinate and Chebyshev points in the second
#' coordinate of a box and normalizes it. Each row of the density (one node
#' of the first coordinate) is represented by its Chebyshev series in the
#' second coordinate, so the posterior probability below a threshold that
#' depends on the first coordinate is obtained without new evaluations of
#' the density.
#'
#' @param logpost Function that takes a list of two node vectors and returns
#'   the matrix of log posterior densities on their product grid (rows: first
#'   coordinate).
#' @param box List with elements \code{lo} and \code{hi} (length 2).
#' @param n_u Number of Gauss-Legendre nodes in the first coordinate.
#' @param n_v Number of Chebyshev points in the second coordinate.
#'
#' @return A list with the nodes \code{u} and \code{v}, the matrix \code{wq}
#'   of quadrature weights times the normalized posterior density (its sum is
#'   1), and two functions: \code{cdf_v(vq)}, the posterior distribution
#'   function of the second coordinate, and \code{cdf_rows(s)}, the posterior
#'   probability that the second coordinate is below \code{s[i]} at the
#'   \code{i}-th node of the first coordinate (\code{s} has length
#'   \code{n_u}).
#' @keywords internal
#' @noRd
post2d <- function(logpost, box, n_u, n_v) {
  g <- gl_nodes(n_u)
  u <- (box$hi[1L] - box$lo[1L]) / 2 * g$x + (box$hi[1L] + box$lo[1L]) / 2
  gw <- g$w * (box$hi[1L] - box$lo[1L]) / 2
  half <- (box$hi[2L] - box$lo[2L]) / 2
  mid <- (box$hi[2L] + box$lo[2L]) / 2
  v <- half * cheb_nodes(n_v) + mid
  L <- logpost(list(u, v)) + log(gw)
  L[is.na(L)] <- -Inf
  if (!any(is.finite(L))) {
    stop("The posterior density is zero on the whole integration box.",
         call. = FALSE)
  }
  W <- exp(L - max(L))
  ccw <- cc_weights(n_v) * half
  W <- W / sum(W %*% ccw)
  cm <- cheb_coef_matrix(n_v)
  anti_v <- cheb_antideriv(cm %*% colSums(W))
  anti_r <- cheb_antideriv(cm %*% t(W))
  list(
    u = u,
    v = v,
    wq = W * rep(ccw, each = n_u),
    cdf_v = function(vq) {
      # one copy of the series per point (cheb_cumint takes one limit per
      # column)
      rep_anti <- list(B = anti_v$B[, rep(1L, length(vq)), drop = FALSE],
                       base = rep(anti_v$base, length(vq)))
      pmin(pmax(cheb_cumint(rep_anti, (vq - mid) / half) * half, 0), 1)
    },
    cdf_rows = function(s) {
      p <- sum(cheb_cumint(anti_r, (s - mid) / half)) * half
      min(max(p, 0), 1)
    }
  )
}
