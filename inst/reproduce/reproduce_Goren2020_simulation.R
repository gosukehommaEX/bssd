# Simulation of Table 2 of
#   Goren E, Lin L-A, Blaustein RO, Ball G (2020). Bayesian meta-analysis of
#   safety outcomes using blinded clinical trial data. Therapeutic
#   Innovation & Regulatory Science 54:1557-1565.
#   doi:10.1007/s43441-020-00183-7
# with multibin_post() of the bssd package.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Goren2020_simulation.R")
# Output: inst/reproduce/results/Goren2020_simulation.rds (read by
#   reproduce_Goren2020_table2.R). Do not rerun unless the model code changes.
# Run time: about one to two hours on one core.
#
# Design: four studies of 200, 400, 800 and 1600 subjects with
# randomization 1:1 (the number on control is Binomial(N, 1/2), Table 1 of
# the paper); five scenarios of the true control rates; true rate ratios 1,
# 2, 3 and 4; 1000 data sets per setting; homogeneous control priors with
# mean 0.02 and prior sample size 100 (Beta(2, 98)) or 400 (Beta(8, 392));
# log(delta) ~ N(0, 5^2). For each data set the equal-tailed 95% credible
# interval of delta is computed, and the proportions of data sets with the
# true rate ratio below, within and above the interval are saved.

n_rep <- 1000
set.seed(20201124)
sizes <- c(200, 400, 800, 1600)
scen <- list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01, 0.03, 0.03),
             c(0.02, 0.02, 0.02, 0.02), c(0.03, 0.01, 0.03, 0.01),
             c(0.03, 0.03, 0.03, 0.03))
priors <- list("100" = c(2, 98), "400" = c(8, 392))
design <- expand.grid(delta = 1:4, scenario = 1:5, prior_n = c(100, 400))
design$below <- NA_real_
design$within <- NA_real_
design$above <- NA_real_
design$n_failed <- NA_integer_
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(design))) {
  th <- scen[[design$scenario[i]]]
  dl <- design$delta[i]
  pc <- priors[[as.character(design$prior_n[i])]]
  res <- matrix(NA_real_, n_rep, 2)
  for (s in seq_len(n_rep)) {
    n0 <- rbinom(4, sizes, 0.5)
    y <- rbinom(4, n0, th) + rbinom(4, sizes - n0, dl * th)
    fit <- tryCatch(multibin_post(y, sizes, k = 1, prior_c = pc,
                                  control = list(check = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      res[s, ] <- fit$rr[c("lower", "upper")]
    }
  }
  ok <- !is.na(res[, 1])
  design$below[i] <- mean(dl < res[ok, 1])
  design$above[i] <- mean(dl > res[ok, 2])
  design$within[i] <- 1 - design$below[i] - design$above[i]
  design$n_failed[i] <- sum(!ok)
  cat(sprintf("%d of %d done (%.0f s)\n", i, nrow(design),
              proc.time()[["elapsed"]] - t0))
}
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(design = design, n_rep = n_rep),
        "inst/reproduce/results/Goren2020_simulation.rds")
print(design)
