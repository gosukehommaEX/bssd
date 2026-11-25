#' Result object of the two-arm models with independent priors
#'
#' @param res Result of \code{mix2_post()}.
#' @param model Model name.
#' @param rate0 Critical value of the study drug parameter, or \code{NULL}.
#' @param rr0 Critical value of the ratio.
#' @param cutoff Posterior probability threshold for a signal.
#' @param level Probability of the credible intervals.
#' @param labels Named character vector with the names of the control
#'   parameter (\code{control}), the study drug parameter (\code{drug}) and
#'   the ratio (\code{ratio}), and the descriptions of the tested quantities
#'   (\code{drug_text}, \code{ratio_text}).
#' @param data,prior Lists with the data summary and the prior parameters.
#'
#' @return An object of class \code{bssd_post}.
#' @keywords internal
#' @noRd
mix2_output <- function(res, model, rate0, rr0, cutoff, level, labels, data,
                        prior) {
  tested_rate <- !is.null(rate0)
  prob <- if (tested_rate) res$prob_rate else res$prob_ratio
  est <- res$means
  names(est) <- c(labels[["control"]], labels[["drug"]], labels[["ratio"]],
                  "difference")
  out <- list(
    model = model,
    prob = prob,
    signal = prob > cutoff,
    rr = c(mean = unname(res$means["ratio"]), median = res$q_ratio[2L],
           lower = res$q_ratio[1L], upper = res$q_ratio[3L]),
    rate = c(mean = unname(res$means["drug"]), median = res$q_d[2L],
             lower = res$q_d[1L], upper = res$q_d[3L]),
    probs = c(drug = res$prob_rate, ratio = res$prob_ratio),
    estimates = est,
    incidence = NULL,
    rr_label = if (tested_rate) labels[["drug"]] else labels[["ratio"]],
    rr0 = if (tested_rate) rate0 else rr0,
    rr_name = labels[["ratio"]],
    rate_name = labels[["drug"]],
    description = paste("Tested quantity:",
                        if (tested_rate) labels[["drug_text"]] else
                          labels[["ratio_text"]]),
    cutoff = cutoff,
    level = level,
    t_star = NULL,
    data = data,
    prior = prior,
    integration = list(method = "Gauss-Legendre and Chebyshev quadrature",
                       nodes = res$nodes,
                       box = rbind(lower = res$box$lo, upper = res$box$hi),
                       error = res$error)
  )
  class(out) <- "bssd_post"
  out
}
