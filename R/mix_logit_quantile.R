#' Logit of a Beta mixture as a function of the normal quantile
#'
#' Returns a function that solves \eqn{\Phi^{-1}(F(x)) = z}{qnorm(F(x)) = z}
#' for \eqn{x}{x}, where \eqn{F}{F} is the distribution function of the
#' logit under a mixture of Beta distributions (\code{lmix_cdf}). The
#' equation is written with the logarithm of the smaller tail probability, so
#' that it is accurate in both tails. A table of \eqn{F}{F} on the grid
#' \eqn{x = \mathrm{sign}(t) (e^{|t|} - 1)}{x = sign(t) * (exp(|t|) - 1)},
#' \eqn{t = -60, -59.75, \dots, 60}{t = -60, -59.75, ..., 60}, is computed
#' once and brackets each root (rounding errors that make the table
#' decrease locally are removed with a cumulative maximum); Newton steps that
#' leave the bracket are replaced by bisection.
#'
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#' @param max_iter Maximum number of iterations.
#'
#' @return Function of a vector of normal quantiles that returns the values
#'   of the logit (the ends of the table for normal quantiles beyond them).
#' @keywords internal
#' @noRd
mix_logit_quantile <- function(mix, max_iter = 100L) {
  tg <- seq(-60, 60, by = 0.25)
  xg <- sign(tg) * expm1(abs(tg))
  # f(x) = log F(x) (lower tail) or -log(1 - F(x)) (upper tail), both
  # increasing in x
  f_lo <- lmix_cdf(xg, mix, lower = TRUE)
  f_up <- -lmix_cdf(xg, mix, lower = FALSE)
  # findInterval() needs non-decreasing tables without NA
  f_lo <- cummax(ifelse(is.na(f_lo), -Inf, f_lo))
  f_up <- cummax(ifelse(is.na(f_up), -Inf, f_up))
  n_g <- length(xg)
  function(z) {
    n <- length(z)
    low <- z <= 0
    lt <- pnorm(-abs(z), log.p = TRUE)
    target <- ifelse(low, lt, -lt)
    j <- integer(n)
    j[low] <- findInterval(target[low], f_lo)
    j[!low] <- findInterval(target[!low], f_up)
    x <- numeric(n)
    x[j == 0L] <- xg[1L]
    x[j == n_g] <- xg[n_g]
    act <- which(j > 0L & j < n_g)
    if (length(act) == 0L) {
      return(x)
    }
    a <- xg[j[act]]
    b <- xg[j[act] + 1L]
    xa <- (a + b) / 2
    lo <- low[act]
    ta <- target[act]
    f_of <- function(x) {
      out <- numeric(length(x))
      if (any(lo)) {
        out[lo] <- lmix_cdf(x[lo], mix, lower = TRUE)
      }
      if (any(!lo)) {
        out[!lo] <- -lmix_cdf(x[!lo], mix, lower = FALSE)
      }
      out
    }
    for (iter in seq_len(max_iter)) {
      fx <- f_of(xa)
      f <- fx - ta
      a <- ifelse(f < 0, xa, a)
      b <- ifelse(f < 0, b, xa)
      # derivative of f: density over the tail probability
      d <- exp(lmix_dens(xa, mix) - ifelse(lo, fx, -fx))
      xn <- xa - f / d
      bad <- !is.finite(xn) | xn <= a | xn >= b
      xn[bad] <- (a[bad] + b[bad]) / 2
      done <- abs(xn - xa) <= 1e-12 * (1 + abs(xa))
      xa <- xn
      if (all(done)) {
        break
      }
    }
    x[act] <- xa
    x
  }
}
