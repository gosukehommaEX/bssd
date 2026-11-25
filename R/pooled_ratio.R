#' Randomization ratio of blinded data pooled from several studies
#'
#' Computes the effective randomization ratio of study drug to control for
#' blinded data pooled from studies with different randomization ratios, as
#' described by Waterhouse et al. (2022, Appendix A). The pooled numbers of
#' events and the pooled exposure can then be analyzed with
#' \code{\link{bdribs_post}} and this ratio as \code{k}.
#'
#' @details
#' Study \eqn{j}{j} with exposure \eqn{E_j}{E_j} and randomization ratio
#' \eqn{k_j}{k_j} contributes the expected exposures
#' \eqn{E_j k_j / (k_j + 1)}{E_j * k_j / (k_j + 1)} on study drug and
#' \eqn{E_j / (k_j + 1)}{E_j / (k_j + 1)} on control. The pooled ratio is the
#' ratio of the totals,
#' \deqn{k = \sum_j E_j k_j / (k_j + 1) \Big/ \sum_j E_j / (k_j + 1).}{k = sum(E_j * k_j / (k_j + 1)) / sum(E_j / (k_j + 1)).}
#'
#' @param exposure Exposures of the studies (positive).
#' @param k Randomization ratios of the studies (study drug to control).
#'
#' @return The pooled randomization ratio.
#'
#' @references
#' Waterhouse, B., Hartford, A., Mukhopadhyay, S., Ferguson, R. and
#' Hendrickson, B. A. (2022). Using the Bayesian detection of potential risk
#' using inference on blinded safety data (BDRIBS) method to support the
#' decision to refer an event for unblinded evaluation. \emph{Pharmaceutical
#' Statistics}, 21(2), 372--385. \doi{10.1002/pst.2175}
#'
#' @seealso \code{\link{bdribs_post}}, \code{\link{bdribs_grid}}
#'
#' @examples
#' # Example 3 of Waterhouse et al. (2022): 250 patient-years at 1:1 and
#' # 80 patient-years at 2:1
#' k <- pooled_ratio(c(250, 80), c(1, 2))
#' k
#' bdribs_post(y = 16, exposure = 330, k = k, lambda_hp = 0.031, pt_hp = Inf)
#'
#' @export
pooled_ratio <- function(exposure, k) {
  check_num(exposure, "exposure", len = NULL, lower = 0, lower_open = TRUE)
  check_num(k, "k", len = length(exposure), lower = 0, lower_open = TRUE)
  sum(exposure * k / (k + 1)) / sum(exposure / (k + 1))
}
