#' Quantile function of a mixture of Beta distributions
#'
#' @param p Probabilities.
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#'
#' @return Quantiles.
#' @keywords internal
#' @noRd
qmixbeta <- function(p, mix) {
  vapply(p, function(q) {
    uniroot(function(x) pmixbeta(x, mix) - q, c(0, 1), tol = 1e-12)$root
  }, numeric(1))
}
