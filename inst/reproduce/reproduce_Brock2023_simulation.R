# Simulation of Tables 1 and 3 of
#   Brock K, Chen C, Ho S, Fuller G, Woolfolk J, McShea C, Penard N (2023).
#   A Bayesian method for safety signal detection in ongoing blinded
#   randomised controlled trials. Pharmaceutical Statistics 22(2):378-395.
#   doi:10.1002/pst.2278
# with map_beta() and mapbern_post() of the bssd package, for the columns
# "None" and "High" of across-study heterogeneity.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Brock2023_simulation.R")
# Output: inst/reproduce/results/Brock2023_simulation.rds (read by
#   reproduce_Brock2023.R). Do not rerun unless the model code changes.
# Run time: about three to five hours on one core (96 settings of 1000
#   trials; set n_rep below to a smaller value for a quicker check).
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
# standard error over the trials are saved.

n_rep <- 1000
set.seed(20231124)
design <- expand.grid(hr = c(1, 2), n = c(50, 100, 200, 1000),
                      N = c(120, 300, 1050), het = c("None", "High"),
                      rate = c("Common", "Rare"),
                      stringsAsFactors = FALSE)
arm_sizes <- list("120" = rep(40, 3), "300" = rep(100, 3),
                  "1050" = rep(c(50, 100, 200), each = 3))
mult <- list(None = c(1, 1, 1), High = c(0.5, 1, 1.5))
design$mean_gamma <- NA_real_
design$se_gamma <- NA_real_
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(design))) {
  lam_c <- if (design$rate[i] == "Common") 0.1 else 0.02
  sizes <- arm_sizes[[as.character(design$N[i])]]
  lam_h <- lam_c * rep(mult[[design$het[i]]], length.out = length(sizes))
  gam <- numeric(n_rep)
  for (s in seq_len(n_rep)) {
    hist_events <- rbinom(length(sizes), sizes, 1 - exp(-lam_h))
    prior <- map_beta(events = hist_events, n = sizes,
                      control = list(check = FALSE))
    drug <- rbinom(design$n[i], 1, 0.75)
    p <- 1 - exp(-ifelse(drug == 1, design$hr[i] * lam_c, lam_c))
    y <- rbinom(design$n[i], 1, p)
    gam[s] <- mapbern_post(y, rep(1, design$n[i]), k = 3, prior_c = prior,
                           details = FALSE,
                           control = list(check = FALSE))$prob
  }
  design$mean_gamma[i] <- mean(gam)
  design$se_gamma[i] <- stats::sd(gam) / sqrt(n_rep)
  cat(sprintf("%d of %d done (%.0f s)\n", i, nrow(design),
              proc.time()[["elapsed"]] - t0))
}
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(design = design, n_rep = n_rep),
        "inst/reproduce/results/Brock2023_simulation.rds")
print(design)
