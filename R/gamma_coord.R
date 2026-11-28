#' First coordinate of the integration for a control rate
#'
#' Returns the first coordinate of the two-dimensional integration of
#' \code{mix2_post()} for a control rate with a Gamma prior. With a shape
#' parameter of at least 1, the coordinate is
#' \eqn{u = \log \lambda}{u = log(lambda)}, with initial limits at the 1e-10
#' and 1 - 1e-10 quantiles of the prior. With a smaller shape parameter, the
#' prior density of \eqn{\log \lambda}{log(lambda)} decreases slowly as
#' \eqn{\lambda \to 0}{lambda -> 0}, the posterior, which the blinded data
#' do not bound there, inherits that tail, and the integration box becomes
#' too wide for the Gauss-Legendre rule. The coordinate is then
#' \eqn{u = \Phi^{-1}(F(\lambda))}{u = qnorm(F(lambda))}, where \eqn{F}{F}
#' is the prior distribution function, and the prior density of \eqn{u}{u}
#' is \eqn{\phi(u)}{dnorm(u)}.
#'
#' @param shape,rate Parameters of the Gamma prior.
#'
#' @return List as described in \code{beta_coord}, with \code{coordinate}
#'   \code{"log"} or \code{"normal"}.
#' @keywords internal
#' @noRd
gamma_coord <- function(shape, rate) {
  if (shape >= 1) {
    return(list(
      c_of = exp, u_of = log,
      lprior_u = function(u) {
        shape * log(rate) - lgamma(shape) + shape * u - rate * exp(u)
      },
      lo_u = log(max(qgamma(1e-10, shape, rate), .Machine$double.xmin)),
      hi_u = log(qgamma(1 - 1e-10, shape, rate)),
      coordinate = "log"))
  }
  list(
    c_of = function(u) {
      out <- u
      neg <- !is.na(u) & u <= 0
      pos <- !is.na(u) & u > 0
      out[neg] <- qgamma(pnorm(u[neg], log.p = TRUE), shape, rate,
                         log.p = TRUE)
      out[pos] <- qgamma(pnorm(-u[pos], log.p = TRUE), shape, rate,
                         lower.tail = FALSE, log.p = TRUE)
      out
    },
    u_of = function(x) {
      lf <- pgamma(x, shape, rate, log.p = TRUE)
      ls <- pgamma(x, shape, rate, lower.tail = FALSE, log.p = TRUE)
      ifelse(lf <= log(0.5), qnorm(lf, log.p = TRUE),
             -qnorm(ls, log.p = TRUE))
    },
    lprior_u = function(u) dnorm(u, log = TRUE),
    lo_u = -8.5, hi_u = 8.5, coordinate = "normal")
}
