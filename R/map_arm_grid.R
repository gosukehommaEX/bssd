#' Quadrature grid for the likelihood of one historical arm
#'
#' Evaluates the log-likelihood of the subject-level data of one arm as a
#' function of \eqn{z = \mathrm{logit}(\pi)}{z = logit(pi)} and returns the
#' nodes and log weights of a composite Gauss-Legendre rule over the region
#' where the likelihood is not negligible (within 60 log units of its
#' maximum on \eqn{z \in [-30, 30]}{z in [-30, 30]}).
#'
#' @param y Event indicators.
#' @param t Times at risk.
#' @param n_pi Gauss-Legendre nodes per panel (8 panels).
#'
#' @return A list with \code{l1} (\eqn{\log \pi}{log(pi)}), \code{l0}
#'   (\eqn{\log(1 - \pi)}{log(1 - pi)}) and \code{lw} (log quadrature weight
#'   plus log-likelihood) at the nodes.
#' @keywords internal
#' @noRd
map_arm_grid <- function(y, t, n_pi) {
  pat <- group_patterns(y, t)
  ll <- function(z) {
    tl <- outer(plogis(-z, log.p = TRUE), pat$t)
    ev <- log(-expm1(tl))
    Y <- matrix(pat$y > 0.5, length(z), length(pat$y), byrow = TRUE)
    drop(ifelse(Y, ev, tl) %*% pat$w)
  }
  zg <- seq(-30, 30, by = 0.05)
  lg <- ll(zg)
  keep <- which(lg > max(lg) - 60)
  lo <- zg[max(min(keep) - 1L, 1L)]
  hi <- zg[min(max(keep) + 1L, length(zg))]
  n_panel <- 8L
  g <- gl_nodes(n_pi)
  edges <- seq(lo, hi, length.out = n_panel + 1L)
  half <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-(n_panel + 1L)]) / 2
  z <- as.vector(outer(g$x, half) + rep(mid, each = n_pi))
  w <- rep(g$w, n_panel) * rep(half, each = n_pi)
  list(l1 = plogis(z, log.p = TRUE), l0 = plogis(-z, log.p = TRUE),
       lw = log(w) + ll(z))
}
