#' Summarize a posterior for safety signal detection as a data frame
#'
#' @param object An object of class \code{bssd_post}.
#' @param ... Not used.
#'
#' @return A data frame with one row and the columns \code{model},
#'   \code{prob}, \code{signal}, \code{rr_mean}, \code{rr_median},
#'   \code{rr_lower}, \code{rr_upper}, \code{ip_placebo}, \code{ip_drug},
#'   \code{ip_difference} and \code{error}. The incidence columns are
#'   \code{NA} when \code{t_star} was not given.
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{expbern_post}},
#'   \code{\link{llbern_post}}
#'
#' @examples
#' fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
#'                    pt_hp = 4000, t_star = 1)
#' summary(fit)
#'
#' @export
summary.bssd_post <- function(object, ...) {
  inc <- object$incidence
  if (is.null(inc)) {
    inc <- c(placebo = NA_real_, drug = NA_real_, difference = NA_real_)
  }
  data.frame(model = object$model,
             prob = object$prob,
             signal = object$signal,
             rr_mean = unname(object$rr["mean"]),
             rr_median = unname(object$rr["median"]),
             rr_lower = unname(object$rr["lower"]),
             rr_upper = unname(object$rr["upper"]),
             ip_placebo = unname(inc["placebo"]),
             ip_drug = unname(inc["drug"]),
             ip_difference = unname(inc["difference"]),
             error = object$integration$error,
             stringsAsFactors = FALSE)
}
