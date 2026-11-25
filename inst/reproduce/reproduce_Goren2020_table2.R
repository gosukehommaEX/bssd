# Reproduction of Table 2 of
#   Goren E, Lin L-A, Blaustein RO, Ball G (2020). Bayesian meta-analysis of
#   safety outcomes using blinded clinical trial data. Therapeutic
#   Innovation & Regulatory Science 54:1557-1565.
#   doi:10.1007/s43441-020-00183-7
# from the simulation results of reproduce_Goren2020_simulation.R.
#
# Run from the package root after
# inst/reproduce/reproduce_Goren2020_simulation.R has saved
# inst/reproduce/results/Goren2020_simulation.rds:
#   source("inst/reproduce/reproduce_Goren2020_table2.R")
# Output: inst/reproduce/results/reproduce_Goren2020_table2.csv
# Run time: seconds.
#
# Proportions of the data sets with the true rate ratio below, within and
# above the equal-tailed 95% credible interval; tolerance 0.005 for rounding
# plus three standard errors of the difference of two proportions from 1000
# data sets.

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

sim_file <- "inst/reproduce/results/Goren2020_simulation.rds"
if (!file.exists(sim_file)) {
  stop("Run inst/reproduce/reproduce_Goren2020_simulation.R first.")
}
sim <- readRDS(sim_file)
d <- sim$design
t2 <- list(
  "100" = list(
    below = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0.01, 0.01, 0, 0),
    within = c(0.33, 0.31, 0.44, 0.57, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               0.99, 0.99, 1, 1),
    above = c(0.67, 0.69, 0.56, 0.43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0)),
  "400" = list(
    below = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0.36, 0.5, 0.59, 0.67),
    within = c(0.03, 0.04, 0.03, 0.01, 1, 1, 1, 1, 0.99, 1, 1, 1,
               0.95, 0.99, 1, 1, 0.64, 0.5, 0.41, 0.33),
    above = c(0.97, 0.96, 0.97, 0.99, 0, 0, 0, 0, 0.01, 0, 0, 0,
              0.05, 0.01, 0, 0, 0, 0, 0, 0)))
for (pn in c("100", "400")) {
  for (i in seq_len(20)) {
    sc <- (i - 1) %/% 4 + 1
    dl <- (i - 1) %% 4 + 1
    sel <- d$prior_n == as.numeric(pn) & d$scenario == sc & d$delta == dl
    for (what in c("below", "within", "above")) {
      p_pub <- t2[[pn]][[what]][i]
      p_rep <- d[[what]][sel]
      pbar <- min(max((p_pub + p_rep) / 2, 1 / sim$n_rep), 1 - 1 / sim$n_rep)
      tol <- 0.005 + 3 * sqrt(2 * pbar * (1 - pbar) / sim$n_rep)
      add(sprintf("GO_T2_N%s_S%d_D%d_%s", pn, sc, dl, what), "Table 2",
          sprintf("proportion of the true ratio %s the 95%% CI, prior N = %s, scenario %d, delta = %d",
                  what, pn, sc, dl), p_pub, p_rep, tol,
          sprintf("failed fits: %d", d$n_failed[sel]))
    }
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Goren2020_table2")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Goren2020_table2.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
