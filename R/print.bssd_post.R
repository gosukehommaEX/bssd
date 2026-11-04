#' Print a posterior summary for safety signal detection
#'
#' @param x An object of class \code{bssd_post}.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{expbern_post}},
#'   \code{\link{llbern_post}}
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
  if (!is.null(x$t_star) && x$model == "LL_BERN") {
    cat("Relative risk: ratio of cumulative incidences at t_star =",
        f(x$t_star), "\n")
  } else {
    cat("Relative risk: rate ratio r of study drug to placebo\n")
  }
  cat("P(", x$rr_label, " > ", f(x$rr0), " | data) = ", f(x$prob), "\n",
      sep = "")
  cat("Signal (posterior probability > ", f(x$cutoff), "): ",
      if (x$signal) "yes" else "no", "\n", sep = "")
  cat("Posterior of ", x$rr_label, ": mean ", f(x$rr["mean"]),
      ", median ", f(x$rr["median"]), ", ", f(100 * x$level),
      "% interval (", f(x$rr["lower"]), ", ", f(x$rr["upper"]), ")\n",
      sep = "")
  if (!is.null(x$incidence)) {
    cat("Cumulative incidence at t_star = ", f(x$t_star), ": placebo ",
        f(x$incidence["placebo"]), ", study drug ", f(x$incidence["drug"]),
        ", difference ", f(x$incidence["difference"]), "\n", sep = "")
  }
  invisible(x)
}
