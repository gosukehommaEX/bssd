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

  # log posterior density at the first coordinates u (one per row) and the
  # matrix V of second coordinates (one row per element of u)
  logpost <- function(u, V) {
    nu <- nrow(V)
    nv <- ncol(V)
    lpc <- plogis(u, log.p = TRUE)
    l1pc <- plogis(-u, log.p = TRUE)
    valid <- NULL
    if (metric == "ratio") {
      lpt <- lpc + V
      # the last point of a row is at p_d = 1 up to rounding
      valid <- lpt <= 1e-10
      lpt <- pmin(lpt, 0)
      l1pt <- log1p(-exp(lpt))
    } else if (metric == "odds") {
      lt <- u + V
      lpt <- plogis(lt, log.p = TRUE)
      l1pt <- plogis(-lt, log.p = TRUE)
    } else {
      lpt <- log(V)
      l1pt <- log1p(-V)
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
      lp <- lp + dnorm(V, mu, s, log = TRUE) -
        pnorm((-lpc - mu) / s, log.p = TRUE)
      lp[!valid] <- -Inf
    } else if (metric == "odds") {
      lp <- lp + dnorm(V, mu, s, log = TRUE)
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
  box <- refine_box(function(nodes) {
    logpost(nodes[[1L]], matrix(nodes[[2L]], length(nodes[[1L]]),
                                length(nodes[[2L]]), byrow = TRUE))
  }, lo, hi, ctrl$n_coarse, ctrl$eps, ctrl$max_iter, lim_lo, lim_hi)
  # limits of the second coordinate in each row: R <= 1 / p_c for the ratio
  # (the density of log(R) ends at -log(p_c)), (0, 1) for p_d
  row_lim <- function(u) {
    nu <- length(u)
    switch(metric,
           ratio = list(lo = rep(-Inf, nu), hi = -plogis(u, log.p = TRUE)),
           odds = list(lo = rep(-Inf, nu), hi = rep(Inf, nu)),
           difference = list(lo = rep(0, nu), hi = rep(1, nu)))
  }
  # the probability of the metric above crit has a kink where the threshold
  # meets a limit of the rows: p_c = 1 / crit for the ratio, p_c = 1 - crit
  # or p_c = -crit for the difference
  u_split <- numeric(0)
  if (metric == "ratio" && crit > 1) {
    u_split <- qlogis(1 / crit)
  } else if (metric == "difference" && crit != 0) {
    u_split <- if (crit > 0) qlogis(1 - crit) else qlogis(-crit)
  }
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)

  # grid with the Gauss-Legendre rule of logit(p_c) split at 'splits'
  build <- function(n_u, n_v, splits) {
    gq <- gl_panels(box$lo[1L], box$hi[1L], n_u, splits)
    u <- gq$x
    nu <- length(u)
    lim <- row_lim(u)
    rb <- row_refine(logpost, u, rep(box$lo[2L], nu), rep(box$hi[2L], nu),
                     ctrl$n_coarse, ctrl$eps, ctrl$max_iter, lim$lo, lim$hi)
    V <- row_nodes(rb$lo, rb$hi, n_v)
    if (metric == "difference") {
      # rounding could move the end points out of [0, 1]
      V <- pmin(pmax(V, 0), 1)
    }
    g <- post2d(logpost(u, V), u, gq$w, rb$lo, rb$hi)
    pc <- plogis(u)
    # P(metric <= x)
    cdf_m <- if (metric == "difference") {
      function(x) g$cdf_rows(pc + x)
    } else {
      function(x) g$cdf_rows(log(x))
    }
    list(g = g, u = u, V = V, pc = pc, cdf_m = cdf_m, lo = rb$lo,
         hi = rb$hi, nu = nu)
  }
  # value of logit(p_c) at which P(metric <= x) has a kink, or NA
  kink_of <- function(x) {
    if (metric == "ratio") {
      ifelse(is.finite(x) & x > 1, qlogis(1 / pmax(x, 1)), NA_real_)
    } else if (metric == "difference") {
      ifelse(is.finite(x) & x > 0, qlogis(1 - pmin(abs(x), 1)),
             ifelse(is.finite(x) & x < 0, qlogis(pmin(abs(x), 1)),
                    NA_real_))
    } else {
      rep(NA_real_, length(x))
    }
  }

  fine <- function(n_u, n_v, quantiles) {
    b <- build(n_u, n_v, u_split)
    g <- b$g
    pc <- b$pc
    q <- rep(NA_real_, 3L)
    prob <- 1 - b$cdf_m(crit)
    if (metric == "difference") {
      pt <- b$V
      mm <- b$V - pc
      if (quantiles) {
        q <- inv_cdf(b$cdf_m, probs, -1, 1)
      }
    } else {
      mm <- exp(b$V)
      pt <- if (metric == "ratio") pc * mm else plogis(b$u + b$V)
      if (quantiles) {
        q <- exp(inv_cdf(function(lv) b$cdf_m(exp(lv)), probs, min(b$lo),
                         max(b$hi)))
      }
    }
    if (quantiles) {
      # the quantiles are refined on a grid split at the kinks of their own
      # thresholds
      ks <- kink_of(q)
      sel <- !is.na(ks) & ks > box$lo[1L] & ks < box$hi[1L]
      if (any(sel)) {
        b2 <- build(n_u, n_v, c(u_split, ks[sel]))
        for (j in which(sel)) {
          qj <- if (metric == "difference") {
            inv_cdf(b2$cdf_m, probs[j], q[j] - 0.02, q[j] + 0.02,
                    extend = TRUE)
          } else {
            exp(inv_cdf(function(lv) b2$cdf_m(exp(lv)), probs[j],
                        log(q[j]) - 0.02, log(q[j]) + 0.02, extend = TRUE))
          }
          if (!is.na(qj)) {
            q[j] <- qj
          }
        }
      }
    }
    list(prob = prob, q = q, mean_m = sum(g$wq * mm),
         means = c(p_c = sum(g$wq * pc), p_d = sum(g$wq * pt),
                   metric = sum(g$wq * mm)), n_u = b$nu)
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_u, ctrl$n_v, quantiles = FALSE)
    n_u <- as.integer(round(1.5 * ctrl$n_u))
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_u, n_v, quantiles = details)
    err <- abs(res$prob - res1$prob)
  } else {
    n_v <- ctrl$n_v
    res <- fine(ctrl$n_u, n_v, quantiles = details)
    err <- NA_real_
  }
  nodes <- c(u = res$n_u, v = n_v)
  res$n_u <- NULL
  c(res, list(method = "Gauss-Legendre and Chebyshev quadrature",
              nodes = nodes,
              box = rbind(lower = box$lo, upper = box$hi), error = err))
}
