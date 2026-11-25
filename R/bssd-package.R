#' bssd: Bayesian Safety Signal Detection for Blinded Clinical Trial Data
#'
#' @description
#' The package computes the posterior probability of a safety signal from
#' blinded clinical trial data with published Bayesian methods:
#' \itemize{
#'   \item the Poisson-based model BDRIBS of Mukhopadhyay et al. (2018) with
#'     the total exposure time or the total time at risk, as specified by
#'     Duncan (2026) (\code{\link{bdribs_post}}), and the tools of Waterhouse
#'     et al. (2022) for it (\code{\link{bdribs_grid}},
#'     \code{\link{pooled_ratio}});
#'   \item the binomial model of Ball and Schnell (2016)
#'     (\code{\link{binbeta_post}}) and the model of Goren et al. (2020) for
#'     several studies (\code{\link{multibin_post}});
#'   \item the method of Gould and Wang (2017) (\code{\link{canary_post}});
#'   \item the exponential Bernoulli and log-logistic Bernoulli mixture models
#'     of Duncan (2026) (\code{\link{expbern_post}},
#'     \code{\link{llbern_post}});
#'   \item the exposure-time Poisson model of Schnell and Ball (2016)
#'     (\code{\link{exptime_post}});
#'   \item the meta-analytic predictive prior and mixture likelihood of Brock
#'     et al. (2023) (\code{\link{map_beta}}, \code{\link{mapbern_post}}).
#' }
#' All posterior quantities are obtained by numerical integration, so no
#' Markov chain Monte Carlo sampling is needed.
#'
#' @references
#' Duncan, B. (2026). Robust safety monitoring and signal detection using
#' alternatives to the standard Poisson distribution. \emph{Journal of
#' Biopharmaceutical Statistics}, 36(2), 227--244.
#' \doi{10.1080/10543406.2024.2395532}
#'
#' Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
#' detection of potential risk using inference on blinded safety data.
#' \emph{Pharmaceutical Statistics}, 17(6), 823--834. \doi{10.1002/pst.1898}
#'
#' @keywords internal
#' @useDynLib bssd, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dbinom dnorm integrate optimize pbeta pbinom pgamma
#' @importFrom stats plnorm plogis pnorm ppois qbeta qgamma qlogis qnorm
#' @importFrom stats rexp runif uniroot
"_PACKAGE"
