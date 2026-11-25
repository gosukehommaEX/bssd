#' Print a posterior summary for safety signal detection
#'
#' Summaries that are not available (\code{NA}), such as the posterior mean
#' of \code{\link{multibin_post}} or the quantiles with
#' \code{details = FALSE}, are not shown.
#'
#' @param x An object of class \code{bssd_post}.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{binbeta_post}},
#'   \code{\link{multibin_post}}, \code{\link{canary_post}},
#'   \code{\link{expbern_post}}, \code{\link{llbern_post}},
#'   \code{\link{exptime_post}}, \code{\link{mapbern_post}}
#'
#' @examples
#' fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
#'                    pt_hp = 4000)
#' print(fit)
#'
#' @export
print.bssd_post <- function(x, digits = 4, ...) {
  f <- function(v) format(signif(v, digits))
  cat("Model:", x$model, "\n")
  if (!is.null(x$description)) {
    cat(x$description, "\n", sep = "")
  } else if (!is.null(x$t_star) && x$model == "LL_BERN") {
    cat("Relative risk: ratio of cumulative incidences at t_star =",
        f(x$t_star), "\n")
  } else {
    cat("Relative risk: rate ratio r of study drug to placebo\n")
  }
  cat("P(", x$rr_label, " > ", f(x$rr0), " | data) = ", f(x$prob), "\n",
      sep = "")
  cat("Signal (posterior probability > ", f(x$cutoff), "): ",
      if (x$signal) "yes" else "no", "\n", sep = "")
  rr_name <- if (is.null(x$rr_name)) x$rr_label else x$rr_name
  # summaries that are not available (NA) are not shown
  summ <- function(v) {
    v <- unname(v[c("mean", "median", "lower", "upper")])
    parts <- character(0)
    if (!is.na(v[1L])) {
      parts <- c(parts, paste0("mean ", f(v[1L])))
    }
    if (!is.na(v[2L])) {
      parts <- c(parts, paste0("median ", f(v[2L])))
    }
    if (!is.na(v[3L]) && !is.na(v[4L])) {
      parts <- c(parts, paste0(f(100 * x$level), "% interval (", f(v[3L]),
                               ", ", f(v[4L]), ")"))
    }
    paste(parts, collapse = ", ")
  }
  s_rr <- summ(x$rr)
  if (nzchar(s_rr)) {
    cat("Posterior of ", rr_name, ": ", s_rr, "\n", sep = "")
  }
  if (!is.null(x$rate)) {
    s_rate <- summ(x$rate)
    if (nzchar(s_rate)) {
      cat("Posterior of ", x$rate_name, ": ", s_rate, "\n", sep = "")
    }
  }
  if (!is.null(x$incidence)) {
    cat("Cumulative incidence at t_star = ", f(x$t_star), ": placebo ",
        f(x$incidence["placebo"]), ", study drug ", f(x$incidence["drug"]),
        ", difference ", f(x$incidence["difference"]), "\n", sep = "")
  }
  invisible(x)
}
