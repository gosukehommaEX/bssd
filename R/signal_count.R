#' Smallest event count that gives a signal
#'
#' For each number of subjects (or exposure), finds the smallest blinded
#' event count for which the posterior probability exceeds \code{cutoff}.
#' Plotted against the number of subjects, these counts are the boundaries of
#' Ball and Schnell (2016, Figure 1) and the control charts of Gould and Wang
#' (2017, Figure 2 and Table III).
#'
#' @details
#' The posterior probability increases with the event count up to the first
#' count with a signal, but not necessarily beyond it: when the event
#' probabilities are bounded (binomial likelihood), counts near \code{n}
#' imply event probabilities near 1 in both arms, and the probability of a
#' large ratio falls again. The counts 1, 2, 4, 8, ... (up to \code{y_max})
#' are therefore evaluated until the first one with a signal, and the
#' smallest count with a signal is found by bisection between that count and
#' the one before it. A signal at count
#' \eqn{y}{y} means \code{prob > cutoff}, the rule of all functions of the
#' package; Gould and Wang (2017) use \eqn{\ge}{>=}, which differs only when
#' the probability equals \code{cutoff} exactly. The largest acceptable
#' count of Ball and Schnell (2016) is the result minus 1.
#'
#' @param n Numbers of subjects (or exposures for the Poisson likelihood of
#'   \code{\link{canary_post}}).
#' @param method \code{"canary"} (\code{\link{canary_post}}) or
#'   \code{"binbeta"} (\code{\link{binbeta_post}}).
#' @param ... Further arguments of the posterior function (priors, \code{k},
#'   critical values); \code{y}, \code{n}, \code{cutoff} and \code{details}
#'   are set here.
#' @param cutoff Posterior probability threshold for a signal.
#' @param y_max Largest count to consider; by default \code{n} (rounded
#'   down). It must be given for the Poisson likelihood.
#'
#' @return A data frame with columns \code{n}, \code{y_crit} (the smallest
#'   count with a signal, \code{NA} if none of the evaluated counts up to
#'   \code{y_max} gives one), \code{prob} (the posterior probability at
#'   \code{y_crit}, or at \code{y_max} if \code{y_crit} is \code{NA}) and
#'   \code{prob_below} (at \code{y_crit - 1}, \code{NA} if \code{y_crit} is 0).
#'
#' @references
#' Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
#' the FDA IND reporting final rule. In J. Lin et al. (eds.),
#' \emph{Statistical Applications from Clinical Trials and Personalized
#' Medicine to Finance and Business Analytics}, ICSA Book Series in
#' Statistics, 201--209. Springer, Cham.
#' \doi{10.1007/978-3-319-42568-9_17}
#'
#' Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
#' rate differences using data from blinded trials: the canary in the coal
#' mine. \emph{Statistics in Medicine}, 36(1), 92--104.
#' \doi{10.1002/sim.7129}
#'
#' @seealso \code{\link{count_oc}}
#'
#' @examples
#' # Gould and Wang (2017), Table III: control prior Beta(20, 980), signal if
#' # P(R > 1) > 0.8
#' signal_count(n = 100, method = "canary", prior_c = c(20, 980),
#'              cutoff = 0.8)
#'
#' @export
signal_count <- function(n, method = c("canary", "binbeta"), ...,
                         cutoff = 0.95, y_max = NULL) {
  method <- match.arg(method)
  check_num(n, "n", len = NULL, lower = 0, lower_open = TRUE)
  check_num(cutoff, "cutoff", lower = 0, upper = 1, lower_open = TRUE,
            upper_open = TRUE)
  args <- list(...)
  fixed <- intersect(names(args), c("y", "n", "cutoff", "details"))
  if (length(fixed) > 0L) {
    stop("These arguments are set by signal_count(): ",
         paste(fixed, collapse = ", "), ".", call. = FALSE)
  }
  post <- if (method == "canary") canary_post else binbeta_post
  if (!is.null(y_max)) {
    check_num(y_max, "y_max", len = NULL, lower = 0, integer = TRUE)
  } else if (identical(args$likelihood, "poisson")) {
    stop("'y_max' must be given for the Poisson likelihood.", call. = FALSE)
  }
  y_max <- if (is.null(y_max)) floor(n) else rep_len(round(y_max), length(n))
  prob_at <- function(y, nn) {
    do.call(post, c(list(y = y, n = nn, cutoff = cutoff, details = FALSE),
                    args))$prob
  }
  res <- lapply(seq_along(n), function(j) {
    p0 <- prob_at(0, n[j])
    if (p0 > cutoff) {
      return(data.frame(n = n[j], y_crit = 0, prob = p0,
                        prob_below = NA_real_))
    }
    # counts 1, 2, 4, ... up to the first signal (see Details)
    lo <- 0
    p_lo <- p0
    hi <- 1
    repeat {
      hi <- min(hi, y_max[j])
      p_hi <- prob_at(hi, n[j])
      if (p_hi > cutoff) {
        break
      }
      if (hi >= y_max[j]) {
        return(data.frame(n = n[j], y_crit = NA_real_, prob = p_hi,
                          prob_below = NA_real_))
      }
      lo <- hi
      p_lo <- p_hi
      hi <- 2 * hi
    }
    while (hi - lo > 1) {
      mid <- floor((lo + hi) / 2)
      pm <- prob_at(mid, n[j])
      if (pm > cutoff) {
        hi <- mid
        p_hi <- pm
      } else {
        lo <- mid
        p_lo <- pm
      }
    }
    data.frame(n = n[j], y_crit = hi, prob = p_hi, prob_below = p_lo)
  })
  do.call(rbind, res)
}
