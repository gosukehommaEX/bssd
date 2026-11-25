# Validation of the functions added in bssd 0.2.0 (exptime_post(),
# binbeta_post(), canary_post(), mapbern_post(), map_beta() and
# multibin_post()) against independent Python computations
# (inst/validation/python/ref_*.py): composite Gauss-Legendre rules split
# exactly at the thresholds, an exact finite mixture for binbeta_post(), and
# nested adaptive quadrature for multibin_post(). The reference values are
# collected in python/ref_methods.csv and the data in
# python/data_methods.csv.
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_methods_python.R")
# Output: inst/validation/results/validate_methods_python.csv
# Run time: about one minute.
#
# Tolerances: 2e-6 for posterior probabilities and the MAP prior moments,
# 1e-5 relative (at least 1e-7) for posterior means, quantiles and the MAP
# prior distribution function.

source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")
ref <- utils::read.csv("inst/validation/python/ref_methods.csv",
                       stringsAsFactors = FALSE)
dat <- utils::read.csv("inst/validation/python/data_methods.csv",
                       stringsAsFactors = FALSE)
get_data <- function(name) dat[dat$dataset == name, ]

rows <- list()
add <- function(id, value) {
  r <- ref[ref$id == id, ]
  if (nrow(r) != 1L) {
    stop("Reference value not found: ", id)
  }
  tol <- if (grepl("_prob|_mean$|_sd$", id)) 2e-6 else
    max(1e-5 * abs(r$reference), 1e-7)
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = "Python", quantity = r$quantity,
    published = r$reference, reproduced = value, tolerance = tol,
    note = "", stringsAsFactors = FALSE)
}
add_q <- function(prefix, v) {
  add(paste0(prefix, "_lower"), v[["lower"]])
  add(paste0(prefix, "_median"), v[["median"]])
  add(paste0(prefix, "_upper"), v[["upper"]])
}

# EXPTIME ---------------------------------------------------------------
d <- get_data("E1")
f <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500), prior_d = c(2, 70),
                  rate0 = 0.03, rr0 = 1.2)
add("PY_EXPTIME_E1_prob_rate", f$probs[["drug"]])
add("PY_EXPTIME_E1_prob_ratio", f$probs[["ratio"]])
add("PY_EXPTIME_E1_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_E1_mean_lambda_d", f$estimates[["lambda_d"]])
add("PY_EXPTIME_E1_mean_r", f$estimates[["r"]])
add_q("PY_EXPTIME_E1_q_lambda_d", f$rate)
add_q("PY_EXPTIME_E1_q_r", f$rr)
d <- get_data("E2")
f <- exptime_post(d$y, d$t, k = 1, prior_c = c(50, 300), prior_d = c(1, 5))
add("PY_EXPTIME_E2_prob_ratio", f$prob)
add("PY_EXPTIME_E2_mean_lambda_c", f$estimates[["lambda_c"]])
add("PY_EXPTIME_E2_mean_lambda_d", f$estimates[["lambda_d"]])
add("PY_EXPTIME_E2_mean_r", f$estimates[["r"]])
add_q("PY_EXPTIME_E2_q_lambda_d", f$rate)
add_q("PY_EXPTIME_E2_q_r", f$rr)

# BINBETA -----------------------------------------------------------------
f <- binbeta_post(40, 300, k = 1, prior_c = c(31, 271), prior_d = c(1, 1),
                  p0 = 0.15)
add("PY_BINBETA_B1_prob_rate", f$probs[["drug"]])
add("PY_BINBETA_B1_prob_ratio", f$probs[["ratio"]])
add("PY_BINBETA_B1_mean_theta_c", f$estimates[["theta_c"]])
add("PY_BINBETA_B1_mean_theta_d", f$estimates[["theta_d"]])
add("PY_BINBETA_B1_mean_rr", f$estimates[["rr"]])
add_q("PY_BINBETA_B1_q_theta_d", f$rate)
add_q("PY_BINBETA_B1_q_rr", f$rr)
f <- binbeta_post(12, 80, k = 2, prior_c = c(5, 45), prior_d = c(0.5, 0.5),
                  rr0 = 1.5)
add("PY_BINBETA_B2_prob_ratio", f$prob)
add("PY_BINBETA_B2_mean_theta_c", f$estimates[["theta_c"]])
add("PY_BINBETA_B2_mean_theta_d", f$estimates[["theta_d"]])
add("PY_BINBETA_B2_mean_rr", f$estimates[["rr"]])
add_q("PY_BINBETA_B2_q_theta_d", f$rate)
add_q("PY_BINBETA_B2_q_rr", f$rr)

# CANARY ------------------------------------------------------------------
cases <- list(
  C1 = list(y = 12, n = 100, k = 1, prior_c = c(1, 49), metric = "ratio",
            crit = 1),
  C2 = list(y = 12, n = 100, k = 1, prior_c = c(1, 49),
            metric = "difference", crit = 0.02),
  C3 = list(y = 20, n = 120, k = 2, prior_c = c(5, 45), metric = "odds",
            crit = 1.2),
  C4 = list(y = 250, n = 20000, k = 1, prior_c = c(75.71, 8500),
            metric = "ratio", crit = 2, likelihood = "poisson"))
for (id in names(cases)) {
  f <- do.call(canary_post, cases[[id]])
  add(sprintf("PY_CANARY_%s_prob", id), f$prob)
  add(sprintf("PY_CANARY_%s_mean_metric", id), f$rr[["mean"]])
  if (id == "C4") {
    add("PY_CANARY_C4_mean_lambda_c", f$estimates[["lambda_c"]])
  } else {
    add(sprintf("PY_CANARY_%s_mean_p_c", id), f$estimates[["p_c"]])
  }
  add_q(sprintf("PY_CANARY_%s_q", id), f$rr)
}

# MAPBERN and MAP prior ---------------------------------------------------
d <- get_data("M1")
mix <- data.frame(w = c(0.7, 0.3), shape1 = c(20, 2), shape2 = c(180, 18))
for (cs in list(c("M1", 1), c("M2", 1.5))) {
  f <- mapbern_post(d$y, d$t, k = 3, prior_c = mix,
                    rr0 = as.numeric(cs[2]))
  add(sprintf("PY_MAPBERN_%s_prob_ratio", cs[1]), f$prob)
  add(sprintf("PY_MAPBERN_%s_mean_pi_c", cs[1]), f$estimates[["pi_c"]])
  add(sprintf("PY_MAPBERN_%s_mean_pi_d", cs[1]), f$estimates[["pi_d"]])
  add_q(sprintf("PY_MAPBERN_%s_q_pi_d", cs[1]), f$rate)
  add_q(sprintf("PY_MAPBERN_%s_q_rr", cs[1]), f$rr)
}
for (cs in list(list(id = "P1", e = c(9, 11, 10), n = c(100, 100, 100)),
                list(id = "P2", e = c(3, 9, 1, 6), n = c(40, 60, 30, 80)))) {
  p <- map_beta(events = cs$e, n = cs$n)
  add(sprintf("PY_MAP_%s_mean", cs$id), attr(p, "moments")[["mean"]])
  add(sprintf("PY_MAP_%s_sd", cs$id), attr(p, "moments")[["sd"]])
  for (x in c(0.05, 0.08, 0.1, 0.12, 0.15)) {
    add(sprintf("PY_MAP_%s_cdf_%g", cs$id, x), pmixbeta(x, p))
  }
}
d <- get_data("MAP_SUBJ")
p <- map_beta(y = d$y, t = d$t, study = d$study)
# the reference for subject-level data uses a coarser grid (relative error
# of about 1e-6), so a relative tolerance of 1e-5 is used for both moments
add("PY_MAP_SUBJ_mean", attr(p, "moments")[["mean"]])
add("PY_MAP_SUBJ_sd", attr(p, "moments")[["sd"]])
rows[[length(rows) - 1L]]$tolerance <- 1e-5 * rows[[length(rows) - 1L]]$published
rows[[length(rows)]]$tolerance <- 1e-5 * rows[[length(rows)]]$published

# MULTIBIN ----------------------------------------------------------------
f <- multibin_post(c(6, 10, 25, 47), c(200, 400, 800, 1600), k = 1,
                   prior_c = c(2, 98))
add("PY_MULTIBIN_G1_prob", f$prob)
add("PY_MULTIBIN_G1_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G1_q", f$rr)
f <- multibin_post(c(3, 9), c(150, 300), k = c(1, 2),
                   prior_c = rbind(c(2, 98), c(4, 196)), lrr_prior = c(0, 2),
                   rr0 = 1.5)
add("PY_MULTIBIN_G2_prob", f$prob)
add("PY_MULTIBIN_G2_mean_log_delta", f$estimates[["log_delta"]])
add_q("PY_MULTIBIN_G2_q", f$rr)

missing <- setdiff(ref$id, vapply(rows, function(r) r$id, character(1)))
if (length(missing) > 0L) {
  stop("Reference values not compared: ", paste(missing, collapse = ", "))
}
tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_methods_python")
out <- tab
names(out)[names(out) == "published"] <- "reference"
names(out)[names(out) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "inst/validation/results/validate_methods_python.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
