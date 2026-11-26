# Simulation of Table 2 of
#   Goren E, Lin L-A, Blaustein RO, Ball G (2020). Bayesian meta-analysis of
#   safety outcomes using blinded clinical trial data. Therapeutic
#   Innovation & Regulatory Science 54:1557-1565.
#   doi:10.1007/s43441-020-00183-7
# with multibin_post() of the bssd package.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Goren2020_simulation.R")
# With n_cores > 1 the settings run in parallel, which needs the installed
# package (run devtools::install() first). n_cores and n_rep (data sets per
# setting, default 1000) can be set before sourcing this script, for
# example n_cores <- parallel::detectCores() - 2.
# Output: inst/reproduce/results/Goren2020_simulation.rds (read by
#   reproduce_Goren2020_table2.R). Do not rerun unless the model code changes.
# Run time: about two to three hours on one core (multibin_post() takes
#   about 0.2 seconds per data set, measured on 2026-11-26), or about ten
#   minutes with 20 cores.
#
# Design: four studies of 200, 400, 800 and 1600 subjects with
# randomization 1:1 (the number on control is Binomial(N, 1/2), Table 1 of
# the paper); five scenarios of the true control rates; true rate ratios 1,
# 2, 3 and 4; 1000 data sets per setting; homogeneous control priors with
# mean 0.02 and prior sample size 100 (Beta(2, 98)) or 400 (Beta(8, 392));
# log(delta) ~ N(0, 5^2). For each data set the equal-tailed 95% credible
# interval of delta is computed, and the proportions of data sets with the
# true rate ratio below, within and above the interval are saved. The random
# numbers of setting i start from set.seed(20201124 + i), so the results do
# not depend on n_cores.

if (!exists("n_cores")) n_cores <- 1L
if (!exists("n_rep")) n_rep <- 1000
sizes <- c(200, 400, 800, 1600)
scen <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01, 0.03, 0.03),
             c(0.02, 0.02, 0.02, 0.02), c(0.03, 0.01, 0.03, 0.01),
             c(0.03, 0.03, 0.03, 0.03))
priors <- list("100" = c(2, 98), "400" = c(8, 392))
design <- expand.grid(delta = 1:4, scenario = 1:5, prior_n = c(100, 400))

# proportions below, within and above the interval in setting i, and the
# number of failed fits (left out of the proportions)
run_setting <- function(i) {
  set.seed(20201124L + i)
  th <- scen[[design$scenario[i]]]
  dl <- design$delta[i]
  pc <- priors[[as.character(design$prior_n[i])]]
  res <- matrix(NA_real_, n_rep, 2)
  for (s in seq_len(n_rep)) {
    n0 <- rbinom(4, sizes, 0.5)
    y <- rbinom(4, n0, th) + rbinom(4, sizes - n0, dl * th)
    fit <- tryCatch(bssd::multibin_post(y, sizes, k = 1, prior_c = pc,
                                        control = list(check = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      res[s, ] <- fit$rr[c("lower", "upper")]
    }
  }
  ok <- !is.na(res[, 1])
  below <- mean(dl < res[ok, 1])
  above <- mean(dl > res[ok, 2])
  c(below = below, within = 1 - below - above, above = above,
    n_failed = sum(!ok))
}

t0 <- proc.time()[["elapsed"]]
out <- if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, library(bssd)))
  parallel::clusterExport(cl, c("design", "scen", "priors", "sizes",
                                "n_rep"))
  res_par <- parallel::parLapplyLB(cl, seq_len(nrow(design)), run_setting)
  parallel::stopCluster(cl)
  res_par
} else {
  lapply(seq_len(nrow(design)), function(i) {
    r <- run_setting(i)
    cat(sprintf("%d of %d done (%.0f s)\n", i, nrow(design),
                proc.time()[["elapsed"]] - t0))
    r
  })
}
out <- do.call(rbind, out)
design$below <- out[, "below"]
design$within <- out[, "within"]
design$above <- out[, "above"]
design$n_failed <- as.integer(out[, "n_failed"])
run_time <- proc.time()[["elapsed"]] - t0
cat(sprintf("all %d settings done (%.0f s)\n", nrow(design), run_time))
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(design = design, n_rep = n_rep, run_time = run_time,
             n_cores = n_cores),
        "inst/reproduce/results/Goren2020_simulation.rds")
print(design)
