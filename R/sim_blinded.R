#' Simulate blinded clinical trial data on an adverse event
#'
#' Generates subject-level data of a randomized two-group trial at one or more
#' analysis times, as used in the simulations of Duncan (2026).
#'
#' @details
#' Subjects enter the trial uniformly over \code{[0, accrual]}, and
#' \code{round(n * k / (k + 1))} of them, chosen at random, are assigned to the
#' study drug. The time to the first occurrence of the event follows the
#' exponential distribution (\code{dist = "exponential"}; \code{par_p} and
#' \code{par_d} are the hazard rates) or the log-logistic distribution of
#' \code{\link{ll_cdf}} (\code{dist = "loglogistic"}; \code{par_p} and
#' \code{par_d} are \code{c(mu, sigma)}). When \code{dropout > 0}, the time to
#' dropout follows the exponential distribution with rate \code{dropout}.
#'
#' At the analysis time \eqn{\tau}{tau}, the follow-up time of a subject who
#' entered at time \eqn{e}{e} is the smallest of \code{follow_up},
#' \eqn{\tau - e}{tau - e} and the time to dropout, and the subject has the
#' event if the time to the event is not larger than the follow-up time.
#' Subjects who have not entered by \eqn{\tau}{tau} are not included. The same
#' subjects and event times are used at all analysis times.
#'
#' @param n Total number of subjects.
#' @param k Randomization ratio of study drug to placebo (\code{k:1}).
#' @param accrual Length of the accrual period (0 if all subjects enter at
#'   time 0).
#' @param follow_up Planned follow-up time of each subject.
#' @param analysis_time Calendar times (from the start of accrual) of the
#'   blinded analyses. The default is the time at which all subjects have
#'   completed follow-up.
#' @param dist Distribution of the time to the event: \code{"exponential"} or
#'   \code{"loglogistic"}.
#' @param par_p,par_d Parameters of the distribution for placebo and study
#'   drug: the hazard rate for \code{"exponential"}, or \code{c(mu, sigma)}
#'   for \code{"loglogistic"}.
#' @param dropout Hazard rate of dropout (0 for no dropout).
#'
#' @return A data frame with one row per subject and analysis, and the
#'   columns \code{look} (index of the analysis), \code{analysis_time},
#'   \code{id}, \code{group} (\code{"placebo"} or \code{"drug"}; not used by
#'   the blinded analyses), \code{entry}, \code{follow_up}, \code{event}
#'   (1 if the event occurred during follow-up), \code{event_time} (time from
#'   entry to the event, \code{NA} without the event) and \code{time_at_risk}
#'   (\code{event_time} for subjects with the event, \code{follow_up}
#'   otherwise).
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{sim_oc}}, \code{\link{ll_mu}}
#'
#' @examples
#' set.seed(1)
#' # 52-week incidence of 10% (placebo) and 20% (study drug)
#' d <- sim_blinded(n = 200, k = 1, accrual = 26, follow_up = 52,
#'                  analysis_time = c(39, 78), dist = "exponential",
#'                  par_p = -log(0.9) / 52, par_d = -log(0.8) / 52)
#' table(d$look, d$event)
#'
#' @export
sim_blinded <- function(n, k = 1, accrual, follow_up,
                        analysis_time = accrual + follow_up,
                        dist = c("exponential", "loglogistic"), par_p, par_d,
                        dropout = 0) {
  check_num(n, "n", lower = 2, integer = TRUE)
  check_num(k, "k", lower = 0, lower_open = TRUE)
  check_num(accrual, "accrual", lower = 0)
  check_num(follow_up, "follow_up", lower = 0, lower_open = TRUE)
  check_num(analysis_time, "analysis_time", len = NULL, lower = 0,
            lower_open = TRUE)
  dist <- match.arg(dist)
  if (dist == "exponential") {
    check_num(par_p, "par_p", lower = 0, lower_open = TRUE)
    check_num(par_d, "par_d", lower = 0, lower_open = TRUE)
  } else {
    check_num(par_p, "par_p", len = 2L)
    check_num(par_d, "par_d", len = 2L)
    if (par_p[2L] <= 0 || par_d[2L] <= 0) {
      stop("The scale parameters in 'par_p' and 'par_d' must be positive.",
           call. = FALSE)
    }
  }
  check_num(dropout, "dropout", lower = 0)

  n <- as.integer(round(n))
  n_d <- as.integer(round(n * k / (k + 1)))
  group <- sample(rep(c("drug", "placebo"), c(n_d, n - n_d)))
  is_d <- group == "drug"
  entry <- runif(n, 0, accrual)
  if (dist == "exponential") {
    tev <- rexp(n, ifelse(is_d, par_d, par_p))
  } else {
    mu <- ifelse(is_d, par_d[1L], par_p[1L])
    sg <- ifelse(is_d, par_d[2L], par_p[2L])
    tev <- exp(mu + sg * qlogis(runif(n)))
  }
  tdrop <- if (dropout > 0) rexp(n, dropout) else rep(Inf, n)

  looks <- lapply(seq_along(analysis_time), function(j) {
    tau <- analysis_time[j]
    fu <- pmin(follow_up, tau - entry, tdrop)
    keep <- fu > 0
    ev <- as.integer(tev <= fu)
    out <- data.frame(look = j, analysis_time = tau, id = seq_len(n),
                      group = group, entry = entry, follow_up = fu,
                      event = ev,
                      event_time = ifelse(ev == 1L, tev, NA_real_),
                      time_at_risk = ifelse(ev == 1L, tev, fu),
                      stringsAsFactors = FALSE)
    out[keep, , drop = FALSE]
  })
  out <- do.call(rbind, looks)
  rownames(out) <- NULL
  out
}
