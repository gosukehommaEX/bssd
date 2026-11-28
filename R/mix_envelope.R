#' Envelope of a mixture of Beta distributions
#'
#' Groups the components of a mixture of Beta distributions and replaces each
#' group by one Beta distribution with the smallest shape parameters of the
#' group and the total weight of the group. The density of each component is
#' then at most a moderate multiple of the density of its group, also in the
#' tails, so that the ratio of the density of the mixture to that of the
#' envelope is bounded. Components with a shape parameter below 1 are grouped
#' by the integer parts of the base-4 logarithms of their shape parameters
#' below 1; the others are grouped along the logit of their mean, a group
#' spanning at most the standard deviation of the logit under the group's
#' distribution. Mixtures with at most \code{max_plain} components are
#' returned unchanged.
#'
#' @param mix Data frame with columns \code{w}, \code{shape1} and
#'   \code{shape2}.
#' @param max_plain Largest number of components returned unchanged.
#'
#' @return Data frame with columns \code{w}, \code{shape1} and \code{shape2}.
#' @keywords internal
#' @noRd
mix_envelope <- function(mix, max_plain = 20L) {
  if (nrow(mix) <= max_plain) {
    return(mix)
  }
  s1 <- mix$shape1
  s2 <- mix$shape2
  key <- character(nrow(mix))
  heavy <- s1 < 1 | s2 < 1
  b1 <- floor(log(pmin(s1, 1)) / log(4))
  b2 <- floor(log(pmin(s2, 1)) / log(4))
  key[heavy] <- paste("h", b1[heavy], b2[heavy])
  idx <- which(!heavy)
  if (length(idx) > 0L) {
    la <- log(s1[idx]) - log(s2[idx])
    ord <- order(la)
    o <- idx[ord]
    la_o <- la[ord]
    grp <- integer(length(o))
    g <- 1L
    start <- la_o[1L]
    m1 <- s1[o[1L]]
    m2 <- s2[o[1L]]
    grp[1L] <- g
    for (j in seq_along(o)[-1L]) {
      n1 <- min(m1, s1[o[j]])
      n2 <- min(m2, s2[o[j]])
      if (la_o[j] - start <= sqrt(trigamma(n1) + trigamma(n2))) {
        m1 <- n1
        m2 <- n2
      } else {
        g <- g + 1L
        start <- la_o[j]
        m1 <- s1[o[j]]
        m2 <- s2[o[j]]
      }
      grp[j] <- g
    }
    key[o] <- paste("l", grp)
  }
  f <- factor(key, levels = unique(key))
  data.frame(w = as.vector(tapply(mix$w, f, sum)),
             shape1 = as.vector(tapply(s1, f, min)),
             shape2 = as.vector(tapply(s2, f, min)))
}
