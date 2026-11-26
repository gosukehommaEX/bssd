# Simulation of Table 1 of
#   Schnell PM, Ball G (2016). A Bayesian exposure-time method for clinical
#   trial safety monitoring with blinded data. Therapeutic Innovation &
#   Regulatory Science 50(6):833-838. doi:10.1177/2168479016656702
# with exptime_post() of the bssd package, for the control rate 0.0130.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Schnell2016_simulation.R")
# With n_cores > 1 the trials run in parallel, which needs the installed
# package (run devtools::install() first). n_cores and n_sim (trials per
# setting, default 1000) can be set before sourcing this script, for example
# n_cores <- parallel::detectCores() - 2.
# Output: inst/reproduce/results/Schnell2016_simulation.rds (read by
#   reproduce_Schnell2016_table1.R). Do not rerun unless the model code changes.
# Run time: about 20 hours on one core (the first two treatment rates took
#   11778 and 13620 seconds on 2026-11-25; exptime_post() takes 2 to 7
#   seconds per analysis with 1500 subjects), or about one hour with 20
#   cores.
#
# The design follows the authors' script example/2-simulation-study.R: 1500
# subjects accrued uniformly over 6 months, 18 months of follow-up, interim
# analyses at months 4, 8, 12 and 16, randomization probability 2/3, priors
# Gamma(340, 26064) and Gamma(340 / 100 + 4, 26064 / 100 + 167) (the script
# uses 340 events, the text 339), critical rate 1.2 * 0.013, and 1000 trials
# per setting with seeds 1, ..., 1000 (the same data as the authors'). The
# seed of a trial is set inside schnell_generate(), so the results do not
# depend on n_cores. The threshold is the smallest p for which at most half
# of the trials at the critical rate signal at some analysis, and a signal is
# a posterior probability greater than p at any analysis.

source("inst/reproduce/schnell2016_data.R")
if (!exists("n_cores")) n_cores <- 1L
if (!exists("n_sim")) n_sim <- 1000
rate_c <- 0.013
rate_d <- c(0.9, 1, 1.1, 1.2, 1.3) * rate_c
interim <- c(4, 8, 12, 16)
prior_c <- c(340, 26064)
prior_d <- c(340 / 100 + 4, 26064 / 100 + 167)
critical <- 1.2 * rate_c

# posterior probabilities at the analyses of one trial (NA if a fit failed,
# so that one failure does not stop the parallel run)
run_trial <- function(s, rate_d_j) {
  dat <- schnell_generate(s, 1500, 2 / 3, 6, 18, rate_c, rate_d_j)
  vapply(interim, function(a) {
    d <- schnell_data_at(dat, a)
    tryCatch(
      bssd::exptime_post(d$y, d$t, k = 2, prior_c = prior_c,
                         prior_d = prior_d, rate0 = critical,
                         details = FALSE,
                         control = list(check = FALSE))$prob,
      error = function(e) NA_real_)
  }, numeric(1))
}

if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, library(bssd)))
  parallel::clusterExport(cl, c("schnell_generate", "schnell_data_at",
                                "rate_c", "interim", "prior_c", "prior_d",
                                "critical"))
}
max_prob <- matrix(NA_real_, n_sim, length(rate_d))
probs <- array(NA_real_, c(n_sim, length(interim), length(rate_d)))
t0 <- proc.time()[["elapsed"]]
for (j in seq_along(rate_d)) {
  out <- if (n_cores > 1L) {
    parallel::parLapplyLB(cl, seq_len(n_sim), run_trial,
                          rate_d_j = rate_d[j])
  } else {
    lapply(seq_len(n_sim), run_trial, rate_d_j = rate_d[j])
  }
  probs[, , j] <- do.call(rbind, out)
  max_prob[, j] <- apply(probs[, , j, drop = FALSE], 1L, max)
  cat(sprintf("treatment rate %.4f done, %.0f s\n", rate_d[j],
              proc.time()[["elapsed"]] - t0))
}
if (n_cores > 1L) {
  parallel::stopCluster(cl)
}
n_failed <- sum(is.na(probs))
# quantile(maxprobs, 1 - signal.rate, type = 1) of the authors' code
j_crit <- which.min(abs(rate_d - critical))
threshold <- if (n_failed == 0L) {
  as.numeric(stats::quantile(max_prob[, j_crit], 0.5, type = 1))
} else {
  NA_real_
}
signal_rate <- colMeans(max_prob > 0.62)
res <- list(rate_c = rate_c, rate_d = rate_d, n_sim = n_sim,
            interim = interim, probs = probs, max_prob = max_prob,
            threshold = threshold, signal_rate = signal_rate,
            signal_rate_own_threshold = colMeans(max_prob > threshold),
            n_failed = n_failed,
            run_time = proc.time()[["elapsed"]] - t0, n_cores = n_cores)
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "inst/reproduce/results/Schnell2016_simulation.rds")
print(data.frame(rate_d = rate_d, signal_rate = signal_rate))
cat("threshold:", threshold, "\n")
if (n_failed > 0L) {
  stop(n_failed, " posterior probabilities could not be computed (saved as ",
       "NA in Schnell2016_simulation.rds).", call. = FALSE)
}
