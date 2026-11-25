# Reproduction of Tables 1 and 3 (columns "None" and "High") of
#   Brock K, Chen C, Ho S, Fuller G, Woolfolk J, McShea C, Penard N (2023).
#   A Bayesian method for safety signal detection in ongoing blinded
#   randomised controlled trials. Pharmaceutical Statistics 22(2):378-395.
#   doi:10.1002/pst.2278
# from the simulation results of reproduce_Brock2023_simulation.R.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Brock2023.R")
# Output: inst/reproduce/results/reproduce_Brock2023.csv
# Run time: seconds.
#
# The published values are means (and standard errors) of gamma over 1000
# simulated trials, with two decimals. The tolerance is 0.005 for rounding
# plus three standard errors of the difference of the two Monte Carlo means
# (a published standard error of 0.00 is taken as 0.005). Brock et al.
# fitted the model by MCMC (Stan) and approximated the MAP prior by a
# mixture fitted to the draws (RBesT); bssd integrates both numerically.

source("inst/reproduce/judge.R")
explained <- read_explained()
sim_file <- "inst/reproduce/results/Brock2023_simulation.rds"
if (!file.exists(sim_file)) {
  stop("Run inst/reproduce/reproduce_Brock2023_simulation.R first.")
}
sim <- readRDS(sim_file)$design

# Tables 1 (common) and 3 (rare): rows in the order N (120, 300, 1050),
# n (50, 100, 200, 1000), HR (1, 2); columns None and High
pub <- list(
  Common = list(
    None = c(0.43, 0.72, 0.45, 0.76, 0.46, 0.80, 0.47, 0.82,
             0.44, 0.77, 0.44, 0.83, 0.46, 0.87, 0.48, 0.92,
             0.43, 0.79, 0.47, 0.89, 0.48, 0.95, 0.48, 0.99),
    None_se = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.00,
                0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.00,
                0.01, 0.01, 0.01, 0.00, 0.01, 0.00, 0.01, 0.00),
    High = c(0.41, 0.70, 0.42, 0.74, 0.42, 0.77, 0.43, 0.79,
             0.41, 0.74, 0.40, 0.79, 0.41, 0.84, 0.42, 0.89,
             0.37, 0.74, 0.39, 0.85, 0.38, 0.91, 0.33, 0.98),
    High_se = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.00,
                0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.00,
                0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.00)),
  Rare = list(
    None = c(0.55, 0.65, 0.56, 0.69, 0.57, 0.70, 0.59, 0.73,
             0.45, 0.60, 0.44, 0.62, 0.45, 0.67, 0.45, 0.72,
             0.42, 0.58, 0.41, 0.64, 0.44, 0.73, 0.47, 0.86),
    None_se = c(rep(0.01, 23), 0.00),
    High = c(0.51, 0.62, 0.52, 0.66, 0.53, 0.67, 0.55, 0.70,
             0.43, 0.58, 0.41, 0.59, 0.41, 0.64, 0.41, 0.69,
             0.38, 0.55, 0.37, 0.60, 0.38, 0.68, 0.38, 0.80),
    High_se = rep(0.01, 24)))
grid <- expand.grid(hr = c(1, 2), n = c(50, 100, 200, 1000),
                    N = c(120, 300, 1050))

rows <- list()
for (rate in c("Common", "Rare")) {
  for (het in c("None", "High")) {
    for (j in seq_len(nrow(grid))) {
      sel <- sim$rate == rate & sim$het == het & sim$N == grid$N[j] &
        sim$n == grid$n[j] & sim$hr == grid$hr[j]
      if (sum(sel) != 1L) {
        stop("The simulation results do not contain the setting.")
      }
      se_pub <- max(pub[[rate]][[paste0(het, "_se")]][j], 0.005)
      tol <- 0.005 + 3 * sqrt(se_pub^2 + sim$se_gamma[sel]^2)
      rows[[length(rows) + 1L]] <- data.frame(
        id = sprintf("B_%s_%s_N%d_n%d_HR%d", substr(rate, 1, 1),
                     substr(het, 1, 1), grid$N[j], grid$n[j], grid$hr[j]),
        source = if (rate == "Common") "Table 1" else "Table 3",
        quantity = sprintf(
          "mean gamma, %s events, %s heterogeneity, N = %d, n = %d, HR = %d",
          tolower(rate), tolower(het), grid$N[j], grid$n[j], grid$hr[j]),
        published = pub[[rate]][[het]][j],
        reproduced = sim$mean_gamma[sel], tolerance = tol,
        note = sprintf("bssd Monte Carlo SE %.4f", sim$se_gamma[sel]),
        stringsAsFactors = FALSE)
    }
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "reproduce_Brock2023")
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Brock2023.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
