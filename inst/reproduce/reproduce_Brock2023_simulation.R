# Simulation of Tables 1 and 3 of
#   Brock K, Chen C, Ho S, Fuller G, Woolfolk J, McShea C, Penard N (2023).
#   A Bayesian method for safety signal detection in ongoing blinded
#   randomised controlled trials. Pharmaceutical Statistics 22(2):378-395.
#   doi:10.1002/pst.2278
# with map_beta() and mapbern_post() of the bssd package, for the columns
# "None" and "High" of across-study heterogeneity.
#
# Run from the package root:
#   devtools::load_all()            # or library(bssd)
#   source("inst/reproduce/reproduce_Brock2023_simulation.R")
# With n_cores > 1 the settings run in parallel, which needs the installed
# package (run devtools::install() first). n_cores and n_rep (trials per
# setting, default 1000) can be set before sourcing this script, for example
# n_cores <- parallel::detectCores() - 2. For a short run that measures the
# run time, set n_rep to a small number and out_file to a temporary file.
# Output: inst/reproduce/results/Brock2023_simulation.rds (read by
#   reproduce_Brock2023.R; out_file). Do not rerun unless the model code
#   changes.
# Run time: about 1.5 hours with 20 cores (96 settings of 1000 trials,
#   about one second per trial; measured on 2026-11-27).
#
# Design (Section 2.4): constant hazards with complete follow-up of one year
# in the historical and the ongoing studies; control hazard 0.1 (common
# events, Table 1) or 0.02 (rare events, Table 3); hazard ratio 1 or 2;
# historical arms of 40 x 3 (N = 120), 100 x 3 (N = 300) or
# 50, 50, 50, 100, 100, 100, 200, 200, 200 (N = 1050) subjects with hazards
# equal to the control hazard ("None") or 0.5, 1 and 1.5 times it ("High",
# repeated for each arm size when N = 1050); ongoing study of n = 50, 100,
# 200 or 1000 subjects with randomization 3:1 (study drug to control; simple
# randomization is used here, as the paper does not say otherwise). The
# hyperpriors are a ~ Beta(1/3, 1/3) and b ~ Beta(5, 5), and the prior of the
# study drug probability is Beta(1/3, 1/3). The quantity is
# gamma = P(pi_d > pi_c | data) in each trial; its mean and Monte Carlo
# standard error over the trials are saved. The random numbers of setting i
# start from set.seed(20231124 + i), so the results do not depend on
# n_cores.

if (!exists("n_cores")) n_cores <- 1L
if (!exists("n_rep")) n_rep <- 1000
if (!exists("out_file")) {
  out_file <- "inst/reproduce/results/Brock2023_simulation.rds"
}
design <- expand.grid(hr = c(1, 2), n = c(50, 100, 200, 1000),
                      N = c(120, 300, 1050), het = c("None", "High"),
                      rate = c("Common", "Rare"),
                      stringsAsFactors = FALSE)
arm_sizes <- list("120" = rep(40, 3), "300" = rep(100, 3),
                  "1050" = rep(c(50, 100, 200), each = 3))
mult <- list(None = c(1, 1, 1), High = c(0.5, 1, 1.5))

# mean and standard error of gamma over the trials of setting i, the
# number of trials whose fit failed (left out of the mean, so that one
# failure does not stop the parallel run) and the first error message
run_setting <- function(i) {
  set.seed(20231124L + i)
  lam_c <- if (design$rate[i] == "Common") 0.1 else 0.02
  sizes <- arm_sizes[[as.character(design$N[i])]]
  lam_h <- lam_c * rep(mult[[design$het[i]]], length.out = length(sizes))
  gam <- rep(NA_real_, n_rep)
  first_error <- NA_character_
  for (s in seq_len(n_rep)) {
    hist_events <- rbinom(length(sizes), sizes, 1 - exp(-lam_h))
    drug <- rbinom(design$n[i], 1, 0.75)
    p <- 1 - exp(-ifelse(drug == 1, design$hr[i] * lam_c, lam_c))
    y <- rbinom(design$n[i], 1, p)
    gam[s] <- tryCatch({
      prior <- bssd::map_beta(events = hist_events, n = sizes,
                              control = list(check = FALSE))
      bssd::mapbern_post(y, rep(1, design$n[i]), k = 3, prior_c = prior,
                         details = FALSE, control = list(check = FALSE))$prob
    }, error = function(e) {
      if (is.na(first_error)) {
        first_error <<- sprintf("trial %d: %s", s, conditionMessage(e))
      }
      NA_real_
    })
  }
  ok <- !is.na(gam)
  list(stats = c(mean_gamma = mean(gam[ok]),
                 se_gamma = stats::sd(gam[ok]) / sqrt(sum(ok)),
                 n_failed = sum(!ok)),
       first_error = first_error)
}

t0 <- proc.time()[["elapsed"]]
out <- if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, library(bssd)))
  parallel::clusterExport(cl, c("design", "arm_sizes", "mult", "n_rep"))
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
design$first_error <- vapply(out, function(o) o$first_error, character(1))
out <- do.call(rbind, lapply(out, function(o) o$stats))
design$mean_gamma <- out[, "mean_gamma"]
design$se_gamma <- out[, "se_gamma"]
design$n_failed <- as.integer(out[, "n_failed"])
run_time <- proc.time()[["elapsed"]] - t0
cat(sprintf("all %d settings done (%.0f s)\n", nrow(design), run_time))
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(design = design, n_rep = n_rep, run_time = run_time,
             n_cores = n_cores),
        out_file)
print(design)
