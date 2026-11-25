# Reproduction of Table 1 and of the probability threshold of
#   Schnell PM, Ball G (2016). A Bayesian exposure-time method for clinical
#   trial safety monitoring with blinded data. Therapeutic Innovation &
#   Regulatory Science 50(6):833-838. doi:10.1177/2168479016656702
# from the simulation results of reproduce_Schnell2016_simulation.R.
#
# Run from the package root after
# inst/reproduce/reproduce_Schnell2016_simulation.R has saved
# inst/reproduce/results/Schnell2016_simulation.rds:
#   source("inst/reproduce/reproduce_Schnell2016_table1.R")
# Output: inst/reproduce/results/reproduce_Schnell2016_table1.csv
# Run time: seconds.
#
# The published signal rates and threshold are from 1000 simulated trials
# analyzed with the authors' Gibbs sampler (1000 draws at each analysis);
# bssd uses the same data (seeds 1 to 1000) and the exact posterior of the
# subject-level model. The tolerance of a signal rate is three standard
# errors of the difference of two proportions from 1000 trials, plus 0.0005
# for rounding; that of the threshold is 0.03.

source("inst/reproduce/judge.R")
explained <- read_explained()

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

sim_file <- "inst/reproduce/results/Schnell2016_simulation.rds"
if (!file.exists(sim_file)) {
  stop("Run inst/reproduce/reproduce_Schnell2016_simulation.R first.")
}
sim <- readRDS(sim_file)
add("S_T1_P", "Simulated Case Study",
    "probability threshold: smallest p with at most 50% signals at the critical rate",
    0.62, sim$threshold, 0.03,
    "the published threshold is from 1000 trials analyzed by Gibbs sampling")
t1_pub <- c(0.077, 0.144, 0.275, 0.497, 0.758)
for (j in seq_along(t1_pub)) {
  p <- sim$signal_rate[j]
  se <- sqrt(2 * p * (1 - p) / sim$n_sim)
  add(sprintf("S_T1_%d", j), "Table 1",
      sprintf("early signal, control rate 0.0130, treatment rate %.4f",
              sim$rate_d[j]), t1_pub[j], p, 3 * se + 0.0005,
      sprintf("threshold 0.62; tolerance three standard errors of the difference of two proportions from %d trials",
              sim$n_sim))
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Schnell2016_table1")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab,
                 "inst/reproduce/results/reproduce_Schnell2016_table1.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
