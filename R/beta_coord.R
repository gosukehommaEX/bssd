#' First coordinate of the integration for a control probability
#'
#' Returns the first coordinate of the two-dimensional integration of
#' \code{mix2_post()} for a control probability whose prior is a mixture of
#' Beta distributions. If the total weight of the components with a shape
#' parameter below 1 is at most \code{eps}, the coordinate is
#' \eqn{u = \mathrm{logit}(\pi)}{u = logit(pi)}, with initial limits at the
#' 1e-10 and 1 - 1e-10 quantiles of the components (within [-40, 40]).
#' Otherwise the prior density of the logit decreases slowly in a tail, the
#' posterior, which the blinded data do not bound as \eqn{\pi \to 0}{pi -> 0}
#' or \eqn{\pi \to 1}{pi -> 1}, inherits that tail, and the integration box
#' becomes too wide for the Gauss-Legendre rule. The coordinate is then
#' \eqn{u = \Phi^{-1}(F_e(x))}{u = qnorm(F_e(x))}, \eqn{x = \mathrm{logit}(\pi)}{x = logit(pi)},
#' where \eqn{F_e}{F_e} is the distribution function of the envelope
#' (\code{mix_envelope}) of the mixture. The prior density of \eqn{u}{u} is
#' \eqn{f(x) / f_e(x) \phi(u)}{f(x) / f_e(x) * dnorm(u)}, where \eqn{f}{f} and
#' \eqn{f_e}{f_e} are the densities of the logit under the mixture and the
#' envelope; the ratio is bounded, so the tails of \eqn{u}{u} are like those
#' of the standard normal distribution. For a single Beta distribution, the
#' envelope is the prior itself and the prior density of \eqn{u}{u} is
#' \eqn{\phi(u)}{dnorm(u)}. The logit at the last vector of \eqn{u}{u} is
#' kept, because \code{mix2_post()} evaluates the same nodes more than once.
#'
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#' @param eps Largest total weight of the components with a shape parameter
#'   below 1 for which the logit is used.
#'
#' @return List with the functions \code{c_of} (probability of \eqn{u}{u}),
#'   \code{u_of} (its inverse) and \code{lprior_u} (log prior density of
#'   \eqn{u}{u}), the initial limits \code{lo_u} and \code{hi_u}, and
#'   \code{coordinate} (\code{"logit"} or \code{"normal"}).
#' @keywords internal
#' @noRd
beta_coord <- function(mix, eps) {
  heavy <- pmin(mix$shape1, mix$shape2) < 1
  if (sum(mix$w[heavy]) <= eps) {
    return(list(
      c_of = plogis, u_of = qlogis,
      lprior_u = function(u) lmix_dens(u, mix),
      lo_u = max(min(qlogis(qbeta(1e-10, mix$shape1, mix$shape2))), -40),
      hi_u = min(max(qlogis(qbeta(1 - 1e-10, mix$shape1, mix$shape2))), 40),
      coordinate = "logit"))
  }
  env <- mix_envelope(mix)
  x_of_z <- mix_logit_quantile(env)
  cache <- new.env(parent = emptyenv())
  x_of <- function(u) {
    if (!identical(u, cache$u)) {
      cache$u <- u
      cache$x <- x_of_z(u)
    }
    cache$x
  }
  list(
    c_of = function(u) plogis(x_of(u)),
    u_of = function(p) {
      x <- qlogis(p)
      lf <- lmix_cdf(x, env, lower = TRUE)
      ls <- lmix_cdf(x, env, lower = FALSE)
      ifelse(lf <= log(0.5), qnorm(lf, log.p = TRUE),
             -qnorm(ls, log.p = TRUE))
    },
    lprior_u = function(u) {
      x <- x_of(u)
      lmix_dens(x, mix) - lmix_dens(x, env) + dnorm(u, log = TRUE)
    },
    lo_u = -8.5, hi_u = 8.5, coordinate = "normal")
}
