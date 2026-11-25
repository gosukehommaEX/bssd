#' Posterior summaries of a blinded two-arm model with independent priors
#'
#' Internal engine of \code{exptime_post()}, \code{binbeta_post()} and
#' \code{mapbern_post()}. The first coordinate is a transformation of the
#' control parameter (log rate or logit probability), integrated by
#' Gauss-Legendre quadrature. The second coordinate is the prior
#' distribution function of the study drug parameter, on which the posterior
#' density is interpolated at Chebyshev points; in this coordinate the
#' density is proportional to the likelihood times the prior density of the
#' control parameter, so a vague prior of the study drug parameter does not
#' create long tails.
#'
#' @param loglik Function of the vectors of control and study drug parameters
#'   that returns the log-likelihood matrix on their product grid.
#' @param c_of Function giving the control parameter from the first
#'   coordinate.
#' @param lprior_u Log prior density of the first coordinate (including the
#'   Jacobian of the transformation).
#' @param d_of Quantile function of the prior of the study drug parameter.
#' @param F_d Distribution function of the prior of the study drug parameter.
#' @param lo_u,hi_u Initial limits of the first coordinate.
#' @param d_max Upper limit of the study drug parameter (\code{Inf} for a
#'   rate, 1 for a probability).
#' @param rate0 Critical value of the study drug parameter, or \code{NULL}.
#' @param rr0 Critical value of the ratio of the study drug parameter to the
#'   control parameter.
#' @param level Probability of the equal-tailed credible intervals.
#' @param details Logical; if \code{TRUE}, the quantiles are computed.
#' @param ctrl Control list with elements \code{n_coarse}, \code{n_u},
#'   \code{n_v}, \code{check}, \code{eps} and \code{max_iter}.
#'
#' @return A list with the posterior probabilities \code{prob_rate}
#'   (\code{NA} when \code{rate0} is \code{NULL}) and \code{prob_ratio}, the
#'   posterior means \code{means} (control, drug, ratio, difference), the
#'   quantiles \code{q_d} of the study drug parameter and \code{q_ratio} of
#'   the ratio (lower limit, median, upper limit), the integration \code{box},
#'   the numbers of \code{nodes} and the \code{error} estimate.
#' @keywords internal
#' @noRd
mix2_post <- function(loglik, c_of, lprior_u, d_of, F_d, lo_u, hi_u, d_max,
                      rate0, rr0, level, details, ctrl) {
  v_max <- 1 - 1e-15
  logpost <- function(nodes) {
    loglik(c_of(nodes[[1L]]), d_of(pmin(nodes[[2L]], v_max))) +
      lprior_u(nodes[[1L]])
  }
  box <- refine_box(logpost, c(lo_u, 0), c(hi_u, 1), ctrl$n_coarse,
                    ctrl$eps, ctrl$max_iter, lim_lo = c(-Inf, 0),
                    lim_hi = c(Inf, 1))
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  fine <- function(n_u, n_v, quantiles) {
    g <- post2d(logpost, box, n_u, n_v)
    tc <- c_of(g$u)
    td <- d_of(pmin(g$v, v_max))
    # P(theta_d / theta_c <= x): the threshold of the second coordinate
    # differs from row to row
    cdf_ratio <- function(x) g$cdf_rows(F_d(pmin(x * tc, d_max)))
    p_rate <- if (is.null(rate0)) NA_real_ else 1 - g$cdf_v(F_d(rate0))
    p_ratio <- 1 - cdf_ratio(rr0)
    ratio <- outer(1 / tc, td)
    means <- c(control = sum(g$wq * tc),
               drug = sum(g$wq * rep(td, each = n_u)),
               ratio = sum(g$wq * ratio),
               difference = sum(g$wq * outer(-tc, td, "+")))
    q_d <- rep(NA_real_, 3L)
    q_ratio <- rep(NA_real_, 3L)
    if (quantiles) {
      q_d <- d_of(pmin(inv_cdf(g$cdf_v, probs, box$lo[2L], box$hi[2L]),
                       v_max))
      pos <- ratio[g$wq > 0 & is.finite(ratio) & ratio > 0]
      start <- log(range(pos))
      q_ratio <- exp(inv_cdf(function(lx) cdf_ratio(exp(lx)), probs,
                             start[1L], start[2L], extend = TRUE))
    }
    list(prob_rate = p_rate, prob_ratio = p_ratio, means = means,
         q_d = q_d, q_ratio = q_ratio)
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_u, ctrl$n_v, quantiles = FALSE)
    n_u <- as.integer(round(1.5 * ctrl$n_u))
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_u, n_v, quantiles = details)
    err <- max(abs(c(res$prob_rate - res1$prob_rate,
                     res$prob_ratio - res1$prob_ratio)), na.rm = TRUE)
  } else {
    n_u <- ctrl$n_u
    n_v <- ctrl$n_v
    res <- fine(n_u, n_v, quantiles = details)
    err <- NA_real_
  }
  c(res, list(box = box, nodes = c(u = n_u, v = n_v), error = err))
}
