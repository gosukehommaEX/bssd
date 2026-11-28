#' Logarithm of the distribution function of a Beta mixture on the logit scale
#'
#' Computes \eqn{\log P(\mathrm{logit}(\pi) \le x)}{log P(logit(pi) <= x)}
#' (\code{lower = TRUE}) or \eqn{\log P(\mathrm{logit}(\pi) > x)}{log P(logit(pi) > x)}
#' for a mixture of Beta distributions of \eqn{\pi}{pi}. The tail
#' probability \eqn{I_p(a, b)}{I_p(a, b)} of each component, with
#' \eqn{p = \pi}{p = pi} or \eqn{p = 1 - \pi}{p = 1 - pi}, is computed with
#' \code{\link[stats]{pbeta}}. Where \eqn{p}{p} is below \code{plogis(-700)},
#' or where \code{pbeta()} returns \code{-Inf} on the log scale because the
#' probability is below the smallest double, the leading term
#' \eqn{p^a (1 - p)^b / (a B(a, b))}{p^a * (1 - p)^b / (a * B(a, b))} of the
#' series of the incomplete Beta function is used, so that the result is
#' finite and increases with \eqn{x}{x} (lower tail) or decreases with
#' \eqn{x}{x} (upper tail) for any \eqn{x}{x}. Its relative error is of the
#' order of \eqn{p (a + b) / (a + 1)}{p * (a + b) / (a + 1)}; it is used
#' only for probabilities far below those that the integration uses.
#'
#' @param x Vector of values of the logit.
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#' @param lower Logical; lower tail if \code{TRUE}, upper tail otherwise.
#'
#' @return Vector of log probabilities.
#' @keywords internal
#' @noRd
lmix_cdf <- function(x, mix, lower = TRUE) {
  nx <- length(x)
  k <- nrow(mix)
  X <- matrix(x, nx, k)
  S1 <- matrix(mix$shape1, nx, k, byrow = TRUE)
  S2 <- matrix(mix$shape2, nx, k, byrow = TRUE)
  LB <- matrix(lbeta(mix$shape1, mix$shape2), nx, k, byrow = TRUE)
  # tail probability I_p(A, B): p = pi with (A, B) = (shape1, shape2) for
  # the lower tail, p = 1 - pi with (A, B) = (shape2, shape1) for the upper
  # tail
  if (lower) {
    sx <- X
    A <- S1
    B <- S2
  } else {
    sx <- -X
    A <- S2
    B <- S1
  }
  L <- matrix(-Inf, nx, k)
  near <- sx > -700
  L[near] <- pbeta(plogis(sx[near]), A[near], B[near], log.p = TRUE)
  # leading term of the series where p is tiny or pbeta() underflows
  far <- !near | (!is.na(L) & L == -Inf)
  L[far] <- A[far] * plogis(sx[far], log.p = TRUE) +
    B[far] * plogis(-sx[far], log.p = TRUE) - log(A[far]) - LB[far]
  L <- L + matrix(log(mix$w), nx, k, byrow = TRUE)
  # row maxima (max.col() with ties.method = "first" compares exactly)
  m <- L[cbind(seq_len(nx), max.col(L, ties.method = "first"))]
  out <- m + log(rowSums(exp(L - m)))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}
