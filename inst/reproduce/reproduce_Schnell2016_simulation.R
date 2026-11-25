# Simulation of Table 1 of
#   Schnell PM, Ball G (2016). A Bayesian exposure-time method for clinical
#   trial safety monitoring with blinded data. Therapeutic Innovation &
#   Regulatory Science 50(6):833-838. doi:10.1177/2168479016656702
# with exptime_post() of the bssd package, for the control rate 0.0130.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Schnell2016_simulation.R")
# Output: inst/reproduce/results/Schnell2016_simulation.rds (read by
#   reproduce_Schnell2016_table1.R). Do not rerun unless the model code changes.
# Run time: about one to two hours on one core.
#
# The design follows the authors' script example/2-simulation-study.R: 1500
# subjects accrued uniformly over 6 months, 18 months of follow-up, interim
# analyses at months 4, 8, 12 and 16, randomization probability 2/3, priors
# Gamma(340, 26064) and Gamma(340 / 100 + 4, 26064 / 100 + 167) (the script
# uses 340 events, the text 339), critical rate 1.2 * 0.013, and 1000 trials
# per setting with seeds 1, ..., 1000 (the same data as the authors'). The
# threshold is the smallest p for which at most half of the trials at the
# critical rate signal at some analysis, and a signal is a posterior
# probability greater than p at any analysis.

source("inst/reproduce/schnell2016_data.R")
n_sim <- 1000
rate_c <- 0.013
rate_d <- c(0.9, 1, 1.1, 1.2, 1.3) * rate_c
interim <- c(4, 8, 12, 16)
prior_c <- c(340, 26064)
prior_d <- c(340 / 100 + 4, 26064 / 100 + 167)
critical <- 1.2 * rate_c

max_prob <- matrix(NA_real_, n_sim, length(rate_d))
probs <- array(NA_real_, c(n_sim, length(interim), length(rate_d)))
t0 <- proc.time()[["elapsed"]]
for (j in seq_along(rate_d)) {
  for (s in seq_len(n_sim)) {
    dat <- schnell_generate(s, 1500, 2 / 3, 6, 18, rate_c, rate_d[j])
    for (a in seq_along(interim)) {
      d <- schnell_data_at(dat, interim[a])
      probs[s, a, j] <- exptime_post(d$y, d$t, k = 2, prior_c = prior_c,
                                     prior_d = prior_d, rate0 = critical,
                                     details = FALSE,
                                     control = list(check = FALSE))$prob
    }
    max_prob[s, j] <- max(probs[s, , j])
  }
  cat(sprintf("treatment rate %.4f done, %.0f s\n", rate_d[j],
              proc.time()[["elapsed"]] - t0))
}
# quantile(maxprobs, 1 - signal.rate, type = 1) of the authors' code
j_crit <- which.min(abs(rate_d - critical))
threshold <- as.numeric(stats::quantile(max_prob[, j_crit], 0.5, type = 1))
signal_rate <- colMeans(max_prob > 0.62)
res <- list(rate_c = rate_c, rate_d = rate_d, n_sim = n_sim,
            interim = interim, probs = probs, max_prob = max_prob,
            threshold = threshold, signal_rate = signal_rate,
            signal_rate_own_threshold = colMeans(max_prob > threshold))
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "inst/reproduce/results/Schnell2016_simulation.rds")
print(data.frame(rate_d = rate_d, signal_rate = signal_rate))
cat("threshold:", threshold, "\n")
