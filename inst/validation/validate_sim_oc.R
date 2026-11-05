# Validation of the simulations of sim_oc() against an independent
# implementation in Python (python/ref_sim_oc.py). The Python script
# simulates the trials with NumPy from the same design and computes the
# posterior quantities of BDRIBS and EXP_BERN without the code of the
# package. The bssd values are the results of
# inst/reproduce/reproduce_Duncan2026_simulation.R
# (inst/reproduce/results/Duncan2026_simulation.rds, 1000 simulated trials
# per scenario).
#
# Run from the package root (the package itself is not needed):
#   source("inst/validation/validate_sim_oc.R")
# Output: inst/validation/results/validate_sim_oc.csv
# Run time: seconds.
#
# Items: the power of BDRIBS(1) and BDRIBS(2) in Figures 3 and 4 and of
# BDRIBS(2) in Figures 5 and 9 (with the threshold of each figure), and the
# power (threshold 0.95) and the mean of the posterior means of the risk
# difference of EXP_BERN in scenario 5.2 (Figure 8, looks 1 to 3). An item is
# PASS when the difference is within four standard errors of the difference
# of two independent Monte Carlo estimates. For the power, the standard
# errors are computed with (x + 0.5) / (n + 1) in place of the proportion
# x / n, so that a proportion of 0 has a positive standard error.

source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")
sim <- readRDS("inst/reproduce/results/Duncan2026_simulation.rds")
out_by_id <- stats::setNames(sim$out, vapply(sim$out, function(o) o$id, ""))
py <- utils::read.csv("inst/validation/python/ref_sim_oc.csv",
                      stringsAsFactors = FALSE)

task_id <- function(figure, point) {
  switch(as.character(figure),
         "3" = sprintf("F3_S2_P%d", point),
         "4" = sprintf("F4_S3_P%d", point),
         "5" = sprintf("F5_S4_P%d", point),
         "8" = "F8_S5.2_P1",
         "9" = sprintf("F9_S6_P%d", point))
}
se_prop <- function(x, n) {
  q <- (x + 0.5) / (n + 1)
  sqrt(q * (1 - q) / n)
}

rows <- list()
add <- function(id, source, quantity, reference, bssd, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = reference,
    reproduced = bssd, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(py))) {
  q <- py[i, ]
  r <- out_by_id[[task_id(q$figure, q$point)]]$results
  r <- r[r$model == q$model & r$look == q$look, ]
  if (q$quantity == "power") {
    x <- r$prob[!is.na(r$prob)]
    n_b <- length(x)
    v_b <- mean(x > q$cutoff)
    se_b <- se_prop(sum(x > q$cutoff), n_b)
    se_p <- se_prop(round(q$value * q$n_sim), q$n_sim)
    lab <- sprintf("power of %s (threshold %g)", q$model, q$cutoff)
  } else {
    x <- r$ip_difference[!is.na(r$ip_difference)]
    n_b <- length(x)
    v_b <- mean(x)
    se_b <- stats::sd(x) / sqrt(n_b)
    se_p <- q$se
    lab <- sprintf("mean posterior mean of the risk difference of %s",
                   q$model)
  }
  add(sprintf("S_F%s_P%d_L%d_%s_%s", q$figure, q$point, q$look, q$model,
              q$quantity),
      sprintf("Figure %s, point %d, look %d", q$figure, q$point, q$look),
      lab, q$value, v_b, 4 * sqrt(se_b^2 + se_p^2),
      sprintf("Python %d trials, bssd %d trials", q$n_sim, n_b))
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_sim_oc")
self_check(tab)
names(tab)[names(tab) == "published"] <- "reference"
names(tab)[names(tab) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/validation/results/validate_sim_oc.csv",
                 row.names = FALSE)
