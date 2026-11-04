#' Find the region that holds the posterior mass
#'
#' Evaluates the log posterior density on a coarse grid of cell midpoints and
#' shrinks or widens the box in each coordinate until the posterior mass
#' outside the box is below \code{eps} in every marginal.
#'
#' @param fun Function that takes a list of node vectors (one per coordinate)
#'   and returns the array of log posterior densities on their product grid.
#' @param lo,hi Initial lower and upper limits of the box.
#' @param n0 Number of cells per coordinate.
#' @param eps Tail mass that may be left outside the box in each marginal.
#' @param max_iter Maximum number of iterations.
#' @param lim_lo,lim_hi Limits that the box may not cross.
#'
#' @return A list with elements \code{lo}, \code{hi} and \code{iter}.
#' @keywords internal
#' @noRd
refine_box <- function(fun, lo, hi, n0, eps, max_iter, lim_lo, lim_hi) {
  d <- length(lo)
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    h <- (hi - lo) / n0
    nodes <- lapply(seq_len(d), function(i) lo[i] + (seq_len(n0) - 0.5) * h[i])
    L <- fun(nodes)
    L[is.na(L)] <- -Inf
    if (!any(is.finite(L))) {
      stop("The posterior density is zero on the whole integration box.",
           call. = FALSE)
    }
    w <- exp(L - max(L))
    w <- w / sum(w)
    new_lo <- lo
    new_hi <- hi
    for (i in seq_len(d)) {
      m <- if (d == 1L) as.vector(w) else apply(w, i, sum)
      cm <- cumsum(m)
      cr <- rev(cumsum(rev(m)))
      i1 <- which(cm > eps)[1L]
      i2 <- max(which(cr > eps))
      width <- hi[i] - lo[i]
      if (i1 == 1L) {
        new_lo[i] <- max(lo[i] - 0.5 * width, lim_lo[i])
      } else {
        new_lo[i] <- max(lo[i] + (i1 - 2L) * h[i], lim_lo[i])
      }
      if (i2 == n0) {
        new_hi[i] <- min(hi[i] + 0.5 * width, lim_hi[i])
      } else {
        new_hi[i] <- min(lo[i] + (i2 + 1L) * h[i], lim_hi[i])
      }
    }
    change <- max(pmax(abs(new_lo - lo), abs(new_hi - hi)) / (hi - lo))
    lo <- new_lo
    hi <- new_hi
    if (change < 0.05) {
      break
    }
  }
  list(lo = lo, hi = hi, iter = iter)
}
