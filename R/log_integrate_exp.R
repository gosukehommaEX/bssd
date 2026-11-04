#' Logarithm of the integral of an exponentiated function
#'
#' Computes \code{log(integral of exp(lfun(z)) over (lower, upper))} by
#' adaptive quadrature. The integrand is scaled by its maximum, and the range
#' is split around the mode so that a narrow peak is not missed.
#'
#' @param lfun Vectorised function returning the log integrand.
#' @param lower,upper Limits of integration (may be infinite).
#'
#' @return The logarithm of the integral (\code{-Inf} if it is zero).
#' @keywords internal
#' @noRd
log_integrate_exp <- function(lfun, lower = -Inf, upper = Inf) {
  if (lower >= upper) {
    return(-Inf)
  }
  g_lo <- if (is.finite(lower)) lower else min(-100, upper - 200)
  g_hi <- if (is.finite(upper)) upper else max(100, lower + 200)
  grid <- seq(g_lo, g_hi, length.out = 801L)
  lg <- lfun(grid)
  lg[is.na(lg)] <- -Inf
  if (!any(is.finite(lg))) {
    return(-Inf)
  }
  step <- grid[2L] - grid[1L]
  z0 <- grid[which.max(lg)]
  opt <- optimize(lfun, c(max(g_lo, z0 - step), min(g_hi, z0 + step)),
                  maximum = TRUE)
  zm <- if (opt$objective > max(lg)) opt$maximum else z0
  lmax <- max(opt$objective, max(lg))
  # approximate scale from the curvature at the mode
  hh <- 1e-3
  d2 <- (lfun(zm + hh) - 2 * lfun(zm) + lfun(zm - hh)) / hh^2
  w <- if (is.finite(d2) && d2 < 0) 1 / sqrt(-d2) else 1
  find_bound <- function(direction, limit) {
    s <- w
    repeat {
      z <- zm + direction * s
      if ((direction < 0 && z <= limit) || (direction > 0 && z >= limit)) {
        return(limit)
      }
      v <- lfun(z)
      if (is.na(v) || v - lmax < -80 || s > 1e6) {
        return(z)
      }
      s <- 2 * s
    }
  }
  zl <- find_bound(-1, lower)
  zu <- find_bound(1, upper)
  breaks <- sort(unique(c(zl, zm - 5 * w, zm, zm + 5 * w, zu)))
  breaks <- breaks[breaks >= zl & breaks <= zu]
  f <- function(z) {
    v <- exp(lfun(z) - lmax)
    v[is.na(v)] <- 0
    v
  }
  total <- 0
  for (j in seq_len(length(breaks) - 1L)) {
    # the integrand is scaled to a maximum of 1, so abs.tol is relative to it
    total <- total + integrate(f, breaks[j], breaks[j + 1L],
                               rel.tol = 1e-10, abs.tol = 1e-13,
                               subdivisions = 1000L)$value
  }
  lmax + log(total)
}
