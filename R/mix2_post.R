#' Posterior summaries of a blinded two-arm model with independent priors
#'
#' Internal engine of \code{exptime_post()}, \code{binbeta_post()} and
#' \code{mapbern_post()}. The first coordinate \code{u} is a transformation
#' of the control parameter (log rate or logit probability), integrated by
#' Gauss-Legendre quadrature. The second coordinate describes the study drug
#' parameter \eqn{\theta_d}{theta_d} through its prior distribution function
#' \eqn{F_d}{F_d}: either \eqn{v = F_d(\theta_d)}{v = F_d(theta_d)}, in
#' which the posterior density is the likelihood times the prior density of
#' the control parameter, so a vague prior of the study drug parameter does
#' not create long tails, or
#' \eqn{z = \Phi^{-1}(F_d(\theta_d))}{z = qnorm(F_d(theta_d))}, computed on
#' the logarithmic scale, which also represents the far tails of an
#' informative prior. For each node of the first coordinate, the interval of
#' \eqn{z}{z} that holds the conditional posterior mass is located
#' separately (\code{row_refine()}), so the intervals follow the posterior
#' when it lies on a ridge of the blinded likelihood or in a tail of the
#' prior of the study drug parameter. Each row is then interpolated at
#' Chebyshev points in \eqn{v}{v} or in \eqn{z}{z}, whichever gives the
#' faster decay of the Chebyshev coefficients of its density.
#'
#' @param loglik Function of the vector of control parameters (one per row)
#'   and the matrix of study drug parameters (one row per control parameter)
#'   that returns the log-likelihood matrix.
#' @param c_of Function giving the control parameter from the first
#'   coordinate.
#' @param u_of Inverse of \code{c_of}.
#' @param lprior_u Log prior density of the first coordinate (including the
#'   Jacobian of the transformation).
#' @param pfun,qfun Distribution and quantile functions of the prior of the
#'   study drug parameter, with the arguments \code{lower.tail} and
#'   \code{log.p} of \code{stats::pbeta()} and \code{stats::qbeta()}.
#' @param lo_u,hi_u Initial limits of the first coordinate.
#' @param d_max Upper limit of the study drug parameter (\code{Inf} for a
#'   rate, 1 for a probability).
#' @param u_split Values of the first coordinate at which the Gauss-Legendre
#'   rule is split (where the probability of the ratio above \code{rr0} has a
#'   kink), or \code{numeric(0)}. The quantiles of the ratio are refined on
#'   a grid split at the kinks of their own thresholds.
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
#'   the ratio (lower limit, median, upper limit), the integration \code{box}
#'   of \code{(u, z)} located by the coarse grid, the numbers of
#'   \code{nodes} and the \code{error} estimate.
#' @keywords internal
#' @noRd
mix2_post <- function(loglik, c_of, u_of, lprior_u, pfun, qfun, lo_u, hi_u,
                      d_max, u_split, rate0, rr0, level, details, ctrl) {
  # study drug parameter at z = qnorm(F_d(theta_d)), on the log scale
  d_of_z <- function(z) {
    out <- z
    neg <- !is.na(z) & z <= 0
    pos <- !is.na(z) & z > 0
    out[neg] <- qfun(pnorm(z[neg], log.p = TRUE), log.p = TRUE)
    out[pos] <- qfun(pnorm(-z[pos], log.p = TRUE), lower.tail = FALSE,
                     log.p = TRUE)
    out
  }
  # z of a value of the study drug parameter
  z_of_d <- function(x) {
    lf <- pfun(x, log.p = TRUE)
    ls <- pfun(x, lower.tail = FALSE, log.p = TRUE)
    ifelse(lf <= log(0.5), qnorm(lf, log.p = TRUE), -qnorm(ls, log.p = TRUE))
  }
  logpost_z <- function(u, Z) {
    loglik(c_of(u), d_of_z(Z)) + lprior_u(u) + dnorm(Z, log = TRUE)
  }
  box <- refine_box(function(nodes) {
    n1 <- length(nodes[[1L]])
    n2 <- length(nodes[[2L]])
    # logpost_z() of the matrix of z; z is the same in every row, so the
    # study drug parameter is computed once per column
    loglik(c_of(nodes[[1L]]),
           matrix(d_of_z(nodes[[2L]]), n1, n2, byrow = TRUE)) +
      lprior_u(nodes[[1L]]) +
      dnorm(matrix(nodes[[2L]], n1, n2, byrow = TRUE), log = TRUE)
  }, c(lo_u, -8), c(hi_u, 8), ctrl$n_coarse, ctrl$eps, ctrl$max_iter,
  lim_lo = c(-Inf, -Inf), lim_hi = c(Inf, Inf))
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)
  # grid with the Gauss-Legendre rule of the first coordinate split at
  # 'splits'
  build <- function(n_u, n_v, splits) {
    gq <- gl_panels(box$lo[1L], box$hi[1L], n_u, splits)
    u <- gq$x
    n <- length(u)
    rb <- row_refine(logpost_z, u, rep(box$lo[2L], n), rep(box$hi[2L], n),
                     ctrl$n_coarse, ctrl$eps, ctrl$max_iter, rep(-Inf, n),
                     rep(Inf, n))
    zl <- rb$lo
    zh <- rb$hi
    vl <- pnorm(zl)
    vh <- pnorm(zh)
    tc <- c_of(u)
    # densities of each row at Chebyshev points in z and in v = pnorm(z)
    Zn <- row_nodes(zl, zh, n_v)
    # (rounding could move the end points out of [0, 1])
    Vn <- pmin(pmax(row_nodes(vl, vh, n_v), 0), 1)
    DZ <- d_of_z(Zn)
    DV <- qfun(Vn)
    Lz <- loglik(tc, DZ) + dnorm(Zn, log = TRUE)
    Lv <- loglik(tc, DV)
    # v is used only where it is not 0 or 1 in double precision
    vrow <- zh <= 8 & zl >= -37 & row_tail_share(Lv) < row_tail_share(Lz)
    L <- Lz
    L[vrow, ] <- Lv[vrow, ]
    TD <- DZ
    TD[vrow, ] <- DV[vrow, ]
    g <- post2d(L + lprior_u(u), u, gq$w, ifelse(vrow, vl, zl),
                ifelse(vrow, vh, zh))
    # coordinate of the threshold x (one per row) in each row
    s_of <- function(x) {
      x <- rep_len(x, n)
      ifelse(vrow, pfun(x), z_of_d(x))
    }
    list(g = g, tc = tc, TD = TD, s_of = s_of, zl = zl, zh = zh, n = n,
         # P(theta_d / theta_c <= x): the threshold differs from row to row
         cdf_ratio = function(x) g$cdf_rows(s_of(pmin(x * tc, d_max))))
  }

  fine <- function(n_u, n_v, quantiles) {
    b <- build(n_u, n_v, u_split)
    g <- b$g
    tc <- b$tc
    TD <- b$TD
    p_rate <- if (is.null(rate0)) NA_real_ else 1 - g$cdf_rows(b$s_of(rate0))
    p_ratio <- 1 - b$cdf_ratio(rr0)
    ratio <- TD / tc
    means <- c(control = sum(g$wq * tc),
               drug = sum(g$wq * TD),
               ratio = sum(g$wq * ratio),
               difference = sum(g$wq * (TD - tc)))
    q_d <- rep(NA_real_, 3L)
    q_ratio <- rep(NA_real_, 3L)
    if (quantiles) {
      cdf_z <- function(z) g$cdf_rows(b$s_of(d_of_z(z)))
      q_d <- d_of_z(inv_cdf(cdf_z, probs, min(b$zl), max(b$zh)))
      pos <- ratio[g$wq > 0 & is.finite(ratio) & ratio > 0]
      start <- log(range(pos))
      q_ratio <- exp(inv_cdf(function(lx) b$cdf_ratio(exp(lx)), probs,
                             start[1L], start[2L], extend = TRUE))
      # the probability of the ratio below x has a kink at the control
      # parameter d_max / x; the quantiles are refined on a grid split there
      kt <- d_max / q_ratio
      sel <- is.finite(kt) & kt > 0 & kt < 1
      if (is.finite(d_max) && any(sel)) {
        b2 <- build(n_u, n_v, c(u_split, u_of(kt[sel])))
        for (j in which(sel)) {
          lx <- log(q_ratio[j])
          qj <- exp(inv_cdf(function(l) b2$cdf_ratio(exp(l)), probs[j],
                            lx - 0.02, lx + 0.02, extend = TRUE))
          if (!is.na(qj)) {
            q_ratio[j] <- qj
          }
        }
      }
    }
    list(prob_rate = p_rate, prob_ratio = p_ratio, means = means,
         q_d = q_d, q_ratio = q_ratio, n_u = b$n)
  }

  if (isTRUE(ctrl$check)) {
    res1 <- fine(ctrl$n_u, ctrl$n_v, quantiles = FALSE)
    n_u <- as.integer(round(1.5 * ctrl$n_u))
    n_v <- as.integer(round(1.5 * (ctrl$n_v - 1L))) + 1L
    res <- fine(n_u, n_v, quantiles = details)
    err <- max(abs(c(res$prob_rate - res1$prob_rate,
                     res$prob_ratio - res1$prob_ratio)), na.rm = TRUE)
  } else {
    n_v <- ctrl$n_v
    res <- fine(ctrl$n_u, n_v, quantiles = details)
    err <- NA_real_
  }
  nodes <- c(u = res$n_u, v = n_v)
  res$n_u <- NULL
  c(res, list(box = box, nodes = nodes, error = err))
}
