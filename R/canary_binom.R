#' Posterior of the metric under the binomial model of Gould and Wang
#'
#' @param y Total number of subjects with the event.
#' @param n Number of subjects.
#' @param k Randomization ratio.
#' @param prior_c Beta prior of the control probability.
#' @param metric \code{"ratio"}, \code{"difference"} or \code{"odds"}.
#' @param metric_prior \code{c(meanlog, sdlog)} of the lognormal prior.
#' @param crit Critical value of the metric.
#' @param level Probability of the credible interval.
#' @param details Logical; compute the quantiles.
#' @param ctrl Control list.
#'
#' @return A list with \code{prob}, \code{q}, \code{mean_m}, \code{means},
#'   \code{method}, \code{nodes}, \code{box} and \code{error}.
#' @keywords internal
#' @noRd
canary_binom <- function(y, n, k, prior_c, metric, metric_prior, crit, level,
                         details, ctrl) {
  n1 <- round(k / (k + 1) * n)
  n0 <- n - n1
  xs <- seq(max(0, y - n0), min(y, n1))
  lc_x <- lchoose(n1, xs) + lchoose(n0, y - xs)
  ac <- prior_c[[1L]]
  bc <- prior_c[[2L]]
  mu <- metric_prior[[1L]]
  s <- metric_prior[[2L]]
  xlog <- function(x, l) if (x == 0) 0 else x * l

  logpost <- function(nodes) {
    u <- nodes[[1L]]
    v <- nodes[[2L]]
    nu <- length(u)
    nv <- length(v)
    lpc <- plogis(u, log.p = TRUE)
    l1pc <- plogis(-u, log.p = TRUE)
    valid <- NULL
    if (metric == "ratio") {
      lpt <- outer(lpc, v, "+")
      valid <- lpt < 0
      l1pt <- log1p(-pmin(exp(lpt), 1))
    } else if (metric == "odds") {
      lt <- outer(u, v, "+")
      lpt <- plogis(lt, log.p = TRUE)
      l1pt <- plogis(-lt, log.p = TRUE)
    } else {
      lpt <- matrix(log(v), nu, nv, byrow = TRUE)
      l1pt <- matrix(log1p(-v), nu, nv, byrow = TRUE)
    }
    acc <- matrix(-Inf, nu, nv)
    for (i in seq_along(xs)) {
      x <- xs[i]
      term <- lc_x[i] + xlog(x, lpt) + xlog(n1 - x, l1pt) +
        (xlog(y - x, lpc) + xlog(n0 - y + x, l1pc))
      m <- pmax(acc, term)
      d <- -abs(acc - term)
      d[is.nan(d)] <- -Inf
      acc <- m + log1p(exp(d))
    }
    lp <- acc + (ac * lpc + bc * l1pc - lbeta(ac, bc))
    if (metric == "ratio") {
      lp <- lp + rep(dnorm(v, mu, s, log = TRUE), each = nu) -
        pnorm((-lpc - mu) / s, log.p = TRUE)
      lp[!valid] <- -Inf
    } else if (metric == "odds") {
      lp <- lp + rep(dnorm(v, mu, s, log = TRUE), each = nu)
    }
    lp
  }

  lo_u <- max(qlogis(qbeta(1e-10, ac, bc)), -40)
  hi_u <- min(qlogis(qbeta(1 - 1e-10, ac, bc)), 40)
  if (metric == "difference") {
    lo <- c(lo_u, 0)
    hi <- c(hi_u, 1)
    lim_lo <- c(-Inf, 0)
    lim_hi <- c(Inf, 1)
  } else {
    lo <- c(lo_u, mu - 10 * s)
    hi <- c(hi_u, mu + 10 * s)
    lim_lo <- c(-Inf, -Inf)
    lim_hi <- c(Inf, Inf)
  }
  box <- refine_box(logpost, lo, hi, ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
                    lim_lo, lim_hi)
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  fine <- function(n_u, n_v, quantiles) {
    g <- post2d(logpost, box, n_u, n_v)
    pc <- plogis(g$u)
    q <- rep(NA_real_, 3L)
    if (metric == "difference") {
      cdf_m <- function(d) g$cdf_rows(pc + d)
      prob <- 1 - cdf_m(crit)
      pt <- matrix(g$v, n_u, n_v, byrow = TRUE)
      mm <- pt - pc
      if (quantiles) {
        q <- inv_cdf(cdf_m, probs, -1, 1)
      }
    } else {
      prob <- 1 - g$cdf_v(log(crit))
      mm <- matrix(exp(g$v), n_u, n_v, byrow = TRUE)
      pt <- if (metric == "ratio") {
        outer(pc, exp(g$v))
      } else {
        plogis(outer(g$u, g$v, "+"))
      }
      if (quantiles) {
        q <- exp(inv_cdf(g$cdf_v, probs, box$lo[2L], box$hi[2L]))
      }
    }
    list(prob = prob, q = q, mean_m = sum(g$wq * mm),
         means = c(p_c = sum(g$wq * pc), p_d = sum(g$wq * pt),
                   metric = sum(g$wq * mm)))
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_u, ctrl$n_v, quantiles = FALSE)
    n_u <- as.integer(round(1.5 * ctrl$n_u))
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_u, n_v, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_u <- ctrl$n_u
    n_v <- ctrl$n_v
    res <- fine(n_u, n_v, quantiles = details)
    err <- NA_real_
  }
  c(res, list(method = "Gauss-Legendre and Chebyshev quadrature",
              nodes = c(u = n_u, v = n_v),
              box = rbind(lower = box$lo, upper = box$hi), error = err))
}
