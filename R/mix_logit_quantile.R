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
#' decrease locally are removed with a cumulative maximum); the iteration
#' starts from the linear interpolation of the table in the bracket, Newton
#' steps that leave the bracket are replaced by bisection, and each root is
#' iterated until its step is at most
#' \eqn{10^{-12} (1 + |x|)}{1e-12 * (1 + |x|)}.
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
    lo <- low[act]
    ta <- target[act]
    # start from the linear interpolation of the table in the bracket
    fa <- ifelse(lo, f_lo[j[act]], f_up[j[act]])
    fb <- ifelse(lo, f_lo[j[act] + 1L], f_up[j[act] + 1L])
    xa <- a + (ta - fa) / (fb - fa) * (b - a)
    mid <- !is.finite(xa) | xa <= a | xa >= b
    xa[mid] <- (a[mid] + b[mid]) / 2
    f_of <- function(x, lo) {
      out <- numeric(length(x))
      if (any(lo)) {
        out[lo] <- lmix_cdf(x[lo], mix, lower = TRUE)
      }
      if (any(!lo)) {
        out[!lo] <- -lmix_cdf(x[!lo], mix, lower = FALSE)
      }
      out
    }
    # only the roots not yet found are iterated
    todo <- seq_along(xa)
    for (iter in seq_len(max_iter)) {
      x0 <- xa[todo]
      l0 <- lo[todo]
      fx <- f_of(x0, l0)
      f <- fx - ta[todo]
      a0 <- ifelse(f < 0, x0, a[todo])
      b0 <- ifelse(f < 0, b[todo], x0)
      # derivative of f: density over the tail probability
      d <- exp(lmix_dens(x0, mix) - ifelse(l0, fx, -fx))
      xn <- x0 - f / d
      bad <- !is.finite(xn) | xn <= a0 | xn >= b0
      xn[bad] <- (a0[bad] + b0[bad]) / 2
      done <- abs(xn - x0) <= 1e-12 * (1 + abs(x0))
      xa[todo] <- xn
      a[todo] <- a0
      b[todo] <- b0
      todo <- todo[!done]
      if (length(todo) == 0L) {
        break
      }
    }
    x[act] <- xa
    x
  }
}
