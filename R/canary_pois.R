#' Posterior of the rate ratio under the Poisson model of Gould and Wang
#'
#' @param y Total number of events.
#' @param exposure Total exposure.
#' @param k Randomization ratio.
#' @param prior_c Shape and rate of the Gamma prior of the control rate.
#' @param metric_prior \code{c(meanlog, sdlog)} of the lognormal prior of the
#'   rate ratio.
#' @param crit Critical value of the rate ratio.
#' @param level Probability of the credible interval.
#' @param details Logical; compute the quantiles.
#' @param ctrl Control list.
#'
#' @return A list with \code{prob}, \code{q}, \code{mean_m}, \code{means},
#'   \code{method}, \code{nodes}, \code{box} and \code{error}.
#' @keywords internal
#' @noRd
canary_pois <- function(y, exposure, k, prior_c, metric_prior, crit, level,
                        details, ctrl) {
  a <- prior_c[[1L]]
  b <- prior_c[[2L]]
  mu <- metric_prior[[1L]]
  s <- metric_prior[[2L]]
  # expected count per unit of control rate
  m_of <- function(v) exposure * (1 + k * exp(v)) / (k + 1)
  logpost <- function(nodes) {
    v <- nodes[[1L]]
    m <- m_of(v)
    (if (y > 0) y * log(m) else 0) - (a + y) * log(b + m) +
      dnorm(v, mu, s, log = TRUE)
  }
  box <- refine_box(logpost, mu - 10 * s, mu + 10 * s, ctrl$n_coarse,
                    ctrl$eps, ctrl$max_iter, -Inf, Inf)
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  fine <- function(n_v, quantiles) {
    g <- post1d(logpost, box, n_v)
    r <- exp(g$v)
    lam <- (a + y) / (b + m_of(g$v))
    q <- rep(NA_real_, 3L)
    if (quantiles) {
      q <- exp(inv_cdf(g$cdf, probs, box$lo, box$hi))
    }
    list(prob = 1 - g$cdf(log(crit)), q = q, mean_m = sum(g$wq * r),
         means = c(lambda_c = sum(g$wq * lam), lambda_d = sum(g$wq * r * lam),
                   metric = sum(g$wq * r)))
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_v, quantiles = FALSE)
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_v, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_v <- ctrl$n_v
    res <- fine(n_v, quantiles = details)
    err <- NA_real_
  }
  c(res, list(method = "Chebyshev quadrature", nodes = c(v = n_v),
              box = c(lower = box$lo, upper = box$hi), error = err))
}
