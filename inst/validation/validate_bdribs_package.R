# Validation of bdribs_post() against the bdribs package (Mukhopadhyay et al.,
# 2018), which samples the same model by MCMC with JAGS.
#
# Requirements: JAGS 4.x (https://mcmc-jags.sourceforge.io/) and the R
# packages rjags, coda and bdribs. The rjags 4-x package on CRAN does not work
# with JAGS 5.
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/validation/validate_bdribs_package.R")
# Output: inst/validation/results/validate_bdribs_package.csv
# Run time: a few minutes.
#
# bdribs() with inf.type = 1 fixes the background rate at bg.events / bg.pyr,
# which is bdribs_post() with pt_hp = Inf. With inf.type = 2 it uses
# del0 ~ Gamma(bg.events + 0.001, bg.pyr + 0.001), which is bdribs_post()
# with lambda_hp = (bg.events + 0.001) / (bg.pyr + 0.001),
# pt_hp = bg.pyr + 0.001 and df = 1. adj.k = FALSE (the default) uses
# r = p / (k (1 - p)) as in bssd.
#
# An item is PASS when the bssd value is within four Monte Carlo standard
# errors of the MCMC value; the standard errors use the effective sample
# size of the draws (coda::effectiveSize).

for (pkg in c("rjags", "coda", "bdribs")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed (bdribs and rjags also need JAGS).")
  }
}
source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")

cases <- data.frame(
  id = c("B1", "B2", "B3", "B4", "B5", "B6", "B7"),
  y = c(5, 5, 5, 9, 9, 7, 30),
  pyr = c(500, 500, 500, 540, 540, 800, 2000),
  k = c(2, 2, 2, 2, 2, 1, 1),
  bg_events = c(18, 18, 18, 33, 33, 18, 20),
  bg_pyr = c(4000, 4000, 4000, 3000, 3000, 4000, 2000),
  inf_type = c(1, 2, 1, 1, 2, 2, 2),
  a = c(1, 1, 0.5, 1, 1, 1, 0.5),
  b = c(1, 1, 0.5, 1, 1, 1, 0.5),
  stringsAsFactors = FALSE)
thresholds <- c(1, 1.2, 1.5)
n_iter <- 100000
n_chain <- 2  # bdribs() supports two chains only

rows <- list()
add <- function(id, source, quantity, reference, bssd, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = reference,
    reproduced = bssd, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

set.seed(20261104)
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, ]
  draws <- bdribs::bdribs(y = cs$y, pyr = cs$pyr, bg.events = cs$bg_events,
                          bg.pyr = cs$bg_pyr, k = cs$k,
                          p.params = list(a = cs$a, b = cs$b),
                          mc.params = list(burn = 5000, iter = n_iter,
                                           nc = n_chain),
                          inf.type = cs$inf_type, plots = FALSE,
                          prnt = FALSE)
  r <- draws$r
  if (cs$inf_type == 1) {
    lam <- cs$bg_events / cs$bg_pyr
    pt <- Inf
  } else {
    lam <- (cs$bg_events + 0.001) / (cs$bg_pyr + 0.001)
    pt <- cs$bg_pyr + 0.001
  }
  label <- sprintf("y = %g, E = %g, k = %g, inf.type = %g, Beta(%g, %g)",
                   cs$y, cs$pyr, cs$k, cs$inf_type, cs$a, cs$b)
  for (cc in thresholds) {
    fit <- bdribs_post(y = cs$y, exposure = cs$pyr, k = cs$k,
                       lambda_hp = lam, pt_hp = pt, df = 1,
                       p_prior = c(cs$a, cs$b), rr0 = cc)
    ind <- as.numeric(r > cc)
    ess <- coda::effectiveSize(ind)
    p_mc <- mean(ind)
    se <- sqrt(p_mc * (1 - p_mc) / ess)
    add(sprintf("V_%s_P%g", cs$id, cc), "bdribs package",
        sprintf("P(r > %g | y), %s", cc, label), p_mc, fit$prob,
        4 * se, sprintf("MCMC SE %.4f, ESS %.0f", se, ess))
  }
  fit <- bdribs_post(y = cs$y, exposure = cs$pyr, k = cs$k, lambda_hp = lam,
                     pt_hp = pt, df = 1, p_prior = c(cs$a, cs$b))
  # the bssd median is checked by the MCMC proportion of draws below it
  ind <- as.numeric(r <= fit$rr[["median"]])
  ess <- coda::effectiveSize(ind)
  se <- sqrt(0.25 / ess)
  add(sprintf("V_%s_MED", cs$id), "bdribs package",
      sprintf("MCMC proportion of r below the bssd median, %s", label),
      0.5, mean(ind), 4 * se, sprintf("MCMC SE %.4f", se))
  if (cs$b + cs$bg_events > 1) {
    ess <- coda::effectiveSize(r)
    se <- stats::sd(r) / sqrt(ess)
    add(sprintf("V_%s_MEAN", cs$id), "bdribs package",
        sprintf("posterior mean of r, %s", label), mean(r),
        fit$rr[["mean"]], 4 * se, sprintf("MCMC SE %.4f", se))
  }
}

tab <- judge(do.call(rbind, rows), explained)
print_counts(tab, "validate_bdribs_package")
self_check(tab)
names(tab)[names(tab) == "published"] <- "reference"
names(tab)[names(tab) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/validation/results/validate_bdribs_package.csv",
                 row.names = FALSE)
