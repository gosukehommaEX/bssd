#' One-dimensional posterior integration at Chebyshev points
#'
#' Evaluates the log posterior density at Chebyshev points of an interval,
#' normalizes it with Clenshaw-Curtis weights, and returns the posterior
#' distribution function from the Chebyshev series of the density.
#'
#' @param logpost Function that takes a list with one node vector and returns
#'   the log posterior densities at the nodes.
#' @param box List with elements \code{lo} and \code{hi} (length 1).
#' @param n Number of Chebyshev points.
#'
#' @return A list with the nodes \code{v}, the weights \code{wq} (quadrature
#'   weights times the normalized density; their sum is 1) and the functions
#'   \code{cdf(vq)} (distribution function) and \code{dens(vq)} (density,
#'   from the Chebyshev series).
#' @keywords internal
#' @noRd
post1d <- function(logpost, box, n) {
  half <- (box$hi - box$lo) / 2
  mid <- (box$hi + box$lo) / 2
  v <- half * cheb_nodes(n) + mid
  L <- logpost(list(v))
  L[is.na(L)] <- -Inf
  if (!any(is.finite(L))) {
    stop("The posterior density is zero on the whole integration interval.",
         call. = FALSE)
  }
  f <- exp(L - max(L))
  ccw <- cc_weights(n) * half
  f <- f / sum(ccw * f)
  coef <- cheb_coef_matrix(n) %*% f
  anti <- cheb_antideriv(coef)
  list(
    v = v,
    wq = ccw * f,
    cdf = function(vq) {
      # one copy of the series per point (cheb_cumint takes one limit per
      # column)
      rep_anti <- list(B = anti$B[, rep(1L, length(vq)), drop = FALSE],
                       base = rep(anti$base, length(vq)))
      pmin(pmax(cheb_cumint(rep_anti, (vq - mid) / half) * half, 0), 1)
    },
    dens = function(vq) {
      x <- pmin(pmax((vq - mid) / half, -1), 1)
      drop(cos(outer(acos(x), 0:(n - 1L))) %*% coef)
    }
  )
}
