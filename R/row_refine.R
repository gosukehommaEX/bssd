#' Find the interval of the second coordinate that holds the mass of each row
#'
#' For each node of the first coordinate (one row), evaluates the log
#' density at the midpoints of \code{n0} cells of the interval of the row and
#' shrinks or widens the interval until the mass of the row outside it is
#' below \code{eps} on each side. This is \code{refine_box()} applied to the
#' conditional density of each row, so that the intervals follow the
#' posterior when it lies on a ridge or in the tail of a prior.
#'
#' @param logpost Function of the vector \code{u} of first coordinates and
#'   a matrix of second coordinates (one row per element of \code{u}) that
#'   returns the matrix of log densities.
#' @param u Nodes of the first coordinate.
#' @param lo,hi Initial lower and upper limits, one per row.
#' @param n0 Number of cells per row.
#' @param eps Mass that may be left outside the interval on each side.
#' @param max_iter Maximum number of iterations.
#' @param lim_lo,lim_hi Limits that the interval of each row may not cross.
#'
#' @return A list with elements \code{lo} and \code{hi} (the intervals) and
#'   \code{ok} (\code{FALSE} for rows whose density is zero on the whole
#'   interval).
#' @keywords internal
#' @noRd
row_refine <- function(logpost, u, lo, hi, n0, eps, max_iter, lim_lo,
                       lim_hi) {
  w0 <- hi - lo
  # an initial interval outside the limits is moved next to the limit
  bad_hi <- lo >= lim_hi
  bad_lo <- hi <= lim_lo
  lo1 <- pmax(lo, lim_lo)
  hi1 <- pmin(hi, lim_hi)
  lo1[bad_hi] <- pmax(lim_hi[bad_hi] - w0[bad_hi], lim_lo[bad_hi])
  hi1[bad_hi] <- lim_hi[bad_hi]
  lo1[bad_lo] <- lim_lo[bad_lo]
  hi1[bad_lo] <- pmin(lim_lo[bad_lo] + w0[bad_lo], lim_hi[bad_lo])
  lo <- lo1
  hi <- hi1
  n <- length(u)
  ok <- rep(TRUE, n)
  mids <- seq_len(n0) - 0.5
  # cumulative sums along the rows from the left and from the right
  up <- upper.tri(diag(n0), diag = TRUE) * 1
  for (iter in seq_len(max_iter)) {
    h <- (hi - lo) / n0
    L <- logpost(u, lo + outer(h, mids))
    L[is.na(L)] <- -Inf
    m <- apply(L, 1L, max)
    ok <- is.finite(m)
    if (!any(ok)) {
      stop("The posterior density is zero on the whole integration box.",
           call. = FALSE)
    }
    W <- exp(L - m)
    W[!ok, ] <- 0
    W <- W / pmax(rowSums(W), .Machine$double.xmin)
    i1 <- max.col((W %*% up > eps) * 1, ties.method = "first")
    i2 <- max.col((W %*% t(up) > eps) * 1, ties.method = "last")
    width <- hi - lo
    new_lo <- ifelse(i1 == 1L, pmax(lo - 0.5 * width, lim_lo),
                     pmax(lo + (i1 - 2L) * h, lim_lo))
    new_hi <- ifelse(i2 == n0, pmin(hi + 0.5 * width, lim_hi),
                     pmin(lo + (i2 + 1L) * h, lim_hi))
    new_lo[!ok] <- lo[!ok]
    new_hi[!ok] <- hi[!ok]
    change <- max((pmax(abs(new_lo - lo), abs(new_hi - hi)) / width)[ok])
    lo <- new_lo
    hi <- new_hi
    if (change < 0.05) {
      break
    }
  }
  list(lo = lo, hi = hi, ok = ok)
}
