#' Integral over the control probability of one study in the model of Goren
#' et al.
#'
#' For each value \eqn{w = \log \delta}{w = log(delta)}, computes the
#' logarithm of
#' \deqn{\int_0^{\min(1, 1 / \delta)} \theta^{a + y - 1} (1 - \theta)^{b - 1} (1 - c \theta)^{n - y} d\theta / B(a, b),}{integral from 0 to min(1, 1 / delta) of theta^(a + y - 1) * (1 - theta)^(b - 1) * (1 - c * theta)^(n - y) d theta / B(a, b),}
#' with \eqn{c = 1 - r + r \delta}{c = 1 - r + r * delta}, on
#' \eqn{z = \mathrm{logit}(\theta)}{z = logit(theta)}. A grid with step
#' 0.25 over an interval of length 80 that ends at the upper limit of
#' \eqn{z}{z} (or at 40) locates the region where the integrand is within
#' \eqn{e^{-40}}{exp(-40)} of its maximum; this region is integrated by a
#' composite Gauss-Legendre rule with panels of length at most 2 (at least
#' 4 panels).
#'
#' @param w Values of \eqn{\log \delta}{log(delta)}.
#' @param y,n Number of subjects with the event and number of subjects.
#' @param r Probability of study drug.
#' @param a,b Shape parameters of the Beta prior.
#' @param n_theta Number of Gauss-Legendre nodes per panel.
#'
#' @return Numeric vector of the logarithms of the integral.
#' @keywords internal
#' @noRd
multibin_inner <- function(w, y, n, r, a, b, n_theta) {
  nw <- length(w)
  cc <- 1 - r + r * exp(w)
  zmax <- rep(Inf, nw)
  pos <- w > 0
  zmax[pos] <- -log(expm1(w[pos]))
  zhi <- pmin(zmax, 40)
  zlo <- zhi - 80
  h <- function(z) {
    (a + y) * plogis(z, log.p = TRUE) + b * plogis(-z, log.p = TRUE) +
      (if (n > y) (n - y) * log1p(-pmin(cc * plogis(z), 1)) else 0)
  }
  step <- 0.25
  ng <- 321L
  Z <- zlo + outer(rep(1, nw), (seq_len(ng) - 1L) * step)
  H <- h(Z)
  H[is.na(H)] <- -Inf
  hm <- apply(H, 1L, max)
  inside <- (H > hm - 40) * 1
  first <- max.col(inside, ties.method = "first")
  last <- max.col(inside, ties.method = "last")
  rows <- seq_len(nw)
  lo <- pmax(Z[cbind(rows, first)] - step, zlo)
  hi <- pmin(Z[cbind(rows, last)] + step, zhi)
  n_panel <- max(4L, as.integer(ceiling(max(hi - lo) / 2)))
  g <- gl_nodes(n_theta)
  width <- (hi - lo) / n_panel
  Zf <- do.call(cbind, lapply(seq_len(n_panel), function(j) {
    outer(width / 2, g$x) + (lo + (j - 0.5) * width)
  }))
  Wf <- outer(width / 2, rep(g$w, n_panel))
  Hf <- h(Zf)
  Hf[is.na(Hf)] <- -Inf
  mf <- pmax(apply(Hf, 1L, max), hm)
  log(rowSums(exp(Hf - mf) * Wf)) + mf - lbeta(a, b)
}
