#' Total exposure time and total time at risk
#'
#' Computes the two exposure measures used by the BDRIBS models of Duncan
#' (2026) from subject-level data. BDRIBS(1) uses the total exposure time
#' (the sum of the follow-up times), and BDRIBS(2) uses the total time at risk,
#' in which the time after the first occurrence of the event is removed.
#'
#' @param follow_up Numeric vector of follow-up times (positive).
#' @param event Vector of event indicators (1 if the subject had at least one
#'   occurrence of the event, 0 otherwise).
#' @param event_time Numeric vector of times to the first occurrence of the
#'   event. Values for subjects without the event are not used and may be
#'   \code{NA}.
#'
#' @return A named numeric vector with elements \code{exposure} (total
#'   exposure time), \code{time_at_risk} (total time at risk), and
#'   \code{n_event} (number of subjects with the event).
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' @seealso \code{\link{bdribs_post}}
#'
#' @examples
#' time_at_risk(follow_up = c(52, 52, 30, 12),
#'              event = c(1, 0, 1, 0),
#'              event_time = c(10, NA, 25, NA))
#'
#' @export
time_at_risk <- function(follow_up, event, event_time) {
  check_num(follow_up, "follow_up", len = NULL, lower = 0, lower_open = TRUE)
  check_num(event, "event", len = length(follow_up), lower = 0, upper = 1,
            integer = TRUE)
  if (length(event_time) != length(follow_up)) {
    stop("'event_time' must have the same length as 'follow_up'.",
         call. = FALSE)
  }
  ev <- round(event) == 1
  et <- event_time[ev]
  if (any(!is.finite(et)) || any(et <= 0) || any(et > follow_up[ev])) {
    stop("'event_time' must be positive and not larger than 'follow_up' ",
         "for subjects with the event.", call. = FALSE)
  }
  c(exposure = sum(follow_up),
    time_at_risk = sum(follow_up[!ev]) + sum(et),
    n_event = sum(ev))
}
