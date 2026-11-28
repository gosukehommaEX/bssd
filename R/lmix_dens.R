#' Logarithm of the density of a Beta mixture on the logit scale
#'
#' Computes the log density of \eqn{\mathrm{logit}(\pi)}{logit(pi)} for a
#' mixture of Beta distributions of \eqn{\pi}{pi} (the density of
#' \eqn{\pi}{pi} times the Jacobian \eqn{\pi (1 - \pi)}{pi * (1 - pi)}),
#' from \eqn{\log \pi}{log(pi)} and \eqn{\log(1 - \pi)}{log(1 - pi)}, so that
#' the result is finite for any \eqn{x}{x}.
#'
#' @param x Vector of values of the logit.
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#'
#' @return Vector of log densities.
#' @keywords internal
#' @noRd
lmix_dens <- function(x, mix) {
  l1 <- plogis(x, log.p = TRUE)
  l0 <- plogis(-x, log.p = TRUE)
  L <- outer(l1, mix$shape1) + outer(l0, mix$shape2) +
    rep(log(mix$w) - lbeta(mix$shape1, mix$shape2), each = length(x))
  m <- apply(L, 1L, max)
  out <- m + log(rowSums(exp(L - m)))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}
