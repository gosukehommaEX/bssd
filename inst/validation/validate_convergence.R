# Convergence of the numerical integration of expbern_post() and
# llbern_post(): the result with the default control is compared with the
# result on a finer grid, and the difference must not exceed the integration
# error estimate reported with the default control. For LL_BERN, the finer
# grid has 24 Gauss-Legendre nodes and 721 Chebyshev points on the same box
# as the default (a grid with 241 Chebyshev points was used before, but its
# own error was about half of the error of the default for the data set L1).
# For L1, whose posterior of sigma_d has much mass near zero, the default
# result and the finer grid are also compared with an independent reference
# computed in Python (python/ref_llbern_breakpoints.py).
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_convergence.R")
# Output: inst/validation/results/validate_convergence.csv
# Run time: a few minutes; the finer LL_BERN grids need about 1 GB of memory.
#
# An item is PASS when |default - finer| <= max(reported error, 1e-6) for the
# posterior probabilities; the tolerance for the incidences is 1e-6 (EXP_BERN)
# and max(reported error, 1e-4) (LL_BERN). The comparisons of the default
# result with the Python reference use the same tolerances, and the finer
# grid must be within a quarter of the reported error of the Python
# reference.

source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")

rows <- list()
add <- function(id, source, quantity, reference, bssd, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = reference,
    reproduced = bssd, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

exp_cases <- list(
  list(id = "E1", seed = 101, k = 1, n = 150, ip_p = 0.10, ip_d = 0.20,
       tau = 52, a = 0.5, b = 0.5, rr0 = 1),
  list(id = "E2", seed = 102, k = 1, n = 300, ip_p = 0.30, ip_d = 0.30,
       tau = 78, a = 0.5, b = 0.5, rr0 = 1),
  list(id = "E3", seed = 103, k = 2, n = 200, ip_p = 0.05, ip_d = 0.10,
       tau = 40, a = 1, b = 1, rr0 = 1.5),
  list(id = "E4", seed = 104, k = 1, n = 60, ip_p = 0.02, ip_d = 0.02,
       tau = 30, a = 0.5, b = 0.5, rr0 = 1)
)
for (cs in exp_cases) {
  set.seed(cs$seed)
  d <- sim_blinded(n = cs$n, k = cs$k, accrual = 26, follow_up = 52,
                   analysis_time = cs$tau, par_p = -log(1 - cs$ip_p) / 52,
                   par_d = -log(1 - cs$ip_d) / 52)
  args <- list(y = d$event, t = d$follow_up, k = cs$k,
               lambda_hp = -log(1 - cs$ip_p) / 52, pt_hp = cs$n / 2 * 52,
               df = 2, p_prior = c(cs$a, cs$b), rr0 = cs$rr0, t_star = 52,
               details = FALSE)
  f0 <- do.call(expbern_post, args)
  f1 <- do.call(expbern_post, c(args, list(control = list(
    n_coarse = 81, n_lambda = 120, n_p = 121, check = FALSE))))
  lab <- sprintf("%s: n = %d, events = %d", cs$id, length(d$event),
                 sum(d$event))
  tol <- max(f0$integration$error, 1e-6)
  add(sprintf("C_%s_PROB", cs$id), "EXP_BERN",
      sprintf("P(r > %g | data), finer grid vs default, %s", cs$rr0, lab),
      f1$prob, f0$prob, tol,
      sprintf("reported error %.2e", f0$integration$error))
  add(sprintf("C_%s_IPD", cs$id), "EXP_BERN",
      sprintf("study drug incidence, finer grid vs default, %s", lab),
      f1$incidence[["drug"]], f0$incidence[["drug"]], 1e-6)
}

ll_cases <- list(
  list(id = "L1", seed = 201, n = 150, par_p = c(7.3, 0.8),
       par_d = c(6.1, 0.8), follow_up = 104, tau = 117),
  list(id = "L2", seed = 202, n = 200, par_p = c(6.17, 1.6),
       par_d = c(6.17, 1.6), follow_up = 52, tau = 104 / 3),
  list(id = "L3", seed = 203, n = 100, par_p = c(7.3, 0.8),
       par_d = c(7.3, 0.8), follow_up = 104, tau = 78),
  list(id = "L4", seed = 204, n = 375, par_p = c(9.85, 2.0),
       par_d = c(8.35, 2.0), follow_up = 52, tau = 52)
)
for (cs in ll_cases) {
  set.seed(cs$seed)
  d <- sim_blinded(n = cs$n, k = 1, accrual = 26, follow_up = cs$follow_up,
                   analysis_time = cs$tau, dist = "loglogistic",
                   par_p = cs$par_p, par_d = cs$par_d)
  args <- list(y = d$event, t = d$follow_up, k = 1, mu_hp = cs$par_p[1],
               sigma_hp = cs$par_p[2], pt_hp = cs$n / 2 * cs$follow_up / 52,
               df_p = 2, t_star = cs$follow_up, details = FALSE)
  f0 <- do.call(llbern_post, args)
  f1 <- do.call(llbern_post, c(args, list(control = list(
    n_outer = 24, n_inner = 721, check = FALSE))))
  lab <- sprintf("%s: n = %d, events = %d", cs$id, length(d$event),
                 sum(d$event))
  tol <- max(f0$integration$error, 1e-6)
  add(sprintf("C_%s_PROB", cs$id), "LL_BERN",
      sprintf("P(RR(t_star) > 1 | data), finer grid vs default, %s", lab),
      f1$prob, f0$prob, tol,
      sprintf("reported error %.2e", f0$integration$error))
  add(sprintf("C_%s_IPD", cs$id), "LL_BERN",
      sprintf("study drug incidence, finer grid vs default, %s", lab),
      f1$incidence[["drug"]], f0$incidence[["drug"]], max(tol, 1e-4))
  if (cs$id == "L1") {
    # the data must be those used by python/ref_llbern_breakpoints.py
    ref_d <- utils::read.csv("inst/validation/python/data_llbern_L1.csv")
    if (!identical(as.numeric(d$event), as.numeric(ref_d$y)) ||
        !isTRUE(all.equal(d$follow_up, ref_d$t, tolerance = 1e-12))) {
      stop("The data set L1 differs from python/data_llbern_L1.csv.")
    }
    py <- utils::read.csv("inst/validation/python/ref_llbern_breakpoints.csv",
                          stringsAsFactors = FALSE)
    py <- py[py$case == "L1" & py$setting == "reference", ]
    stopifnot(nrow(py) == 1L)
    add("C_L1_PY_PROB", "LL_BERN",
        sprintf("P(RR(t_star) > 1 | data), Python reference vs default, %s",
                lab), py$prob, f0$prob, tol,
        sprintf("reported error %.2e", f0$integration$error))
    add("C_L1_PY_IPD", "LL_BERN",
        sprintf("study drug incidence, Python reference vs default, %s", lab),
        py$ip_drug, f0$incidence[["drug"]], max(tol, 1e-4))
    add("C_L1_PY_FINE", "LL_BERN",
        sprintf("P(RR(t_star) > 1 | data), Python reference vs finer grid, %s",
                lab), py$prob, f1$prob, tol / 4,
        "tolerance: a quarter of the reported error of the default")
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_convergence")
self_check(tab)
names(tab)[names(tab) == "published"] <- "reference"
names(tab)[names(tab) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/validation/results/validate_convergence.csv",
                 row.names = FALSE)
