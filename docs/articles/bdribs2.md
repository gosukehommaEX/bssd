# BDRIBS(2): Poisson model with the total time at risk

``` r

library(bssd)
```

## Model

BDRIBS(2) is the modification of BDRIBS(1) examined by Duncan (2026,
Section 3.2.3): the total exposure time is replaced by the total time at
risk, in which the time after the first occurrence of the event is
discarded. For subject \\i\\ with the follow-up time \\t_i\\ and the
time to the first event \\s_i\\ (\\s_i = \infty\\ for a subject without
the event), the two exposures are \\E_1 = \sum_i t_i \quad
\text{(BDRIBS(1))}, \qquad E_2 = \sum_i \min(t_i, s_i) \quad
\text{(BDRIBS(2))}.\\ The rest of the model is that of BDRIBS(1)
([`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md)):
\\Y \sim \mathrm{Poisson}\\E_2 \delta (r k + 1) / (k + 1)\\, \quad
\delta \sim \mathrm{Gamma}(\lambda\_{HP} PT\_{HP} / DF, PT\_{HP} / DF),
\quad r = \frac{p}{k (1 - p)}, \quad p \sim \mathrm{Beta}(a, b).\\

The number of subjects with the event, \\Y\\, counts first events only,
and a subject cannot have a first event after it. With the total
exposure time, this time is still counted as exposure, so that the
exposure is too large for the observed number of events. Duncan (2026,
Section 3.2.3) notes that the distinction has no impact for small
incidences and a noticeable impact for larger ones, and BDRIBS(1) has
lower power than BDRIBS(2) at medium and high incidences (Duncan, 2026,
Section 4.1).

BDRIBS(2) uses the same function as BDRIBS(1),
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
with the total time at risk as `exposure`. The other arguments are
described in
[`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md).
[`time_at_risk()`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md)
computes both exposures from subject-level data.

## Example

The following data are simulated with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)
with incidences of the high-incidence scenarios of Duncan (2026, Table
1): 200 subjects randomized 1:1 with uniform entry over 26 weeks, 52
weeks of follow-up, and exponential times to the event with cumulative
incidences of 30% on placebo and 45% on the study drug at week 52. The
blinded analysis is at week 78, when all subjects have completed
follow-up.

``` r

set.seed(2026)
rate_of <- function(ip, t) -log(1 - ip) / t
d <- sim_blinded(n = 200, k = 1, accrual = 26, follow_up = 52,
                 analysis_time = 78, dist = "exponential",
                 par_p = rate_of(0.30, 52), par_d = rate_of(0.45, 52))
tar <- time_at_risk(d$follow_up, d$event, d$event_time)
tar
#>     exposure time_at_risk      n_event 
#>    10400.000     7837.746       84.000
```

The total time at risk is 75% of the total exposure time. The prior of
\\\delta\\ is based on a historical placebo group of 100 subjects
followed for 52 weeks with the incidence of 30%, discounted by \\DF =
2\\.

``` r

fits <- list(
  BDRIBS1 = bdribs_post(y = tar[["n_event"]], exposure = tar[["exposure"]],
                        k = 1, lambda_hp = rate_of(0.30, 52),
                        pt_hp = 100 * 52, df = 2, t_star = 52),
  BDRIBS2 = bdribs_post(y = tar[["n_event"]],
                        exposure = tar[["time_at_risk"]], k = 1,
                        lambda_hp = rate_of(0.30, 52), pt_hp = 100 * 52,
                        df = 2, t_star = 52)
)
res <- do.call(rbind, lapply(fits, summary))
res$model <- names(fits)
knitr::kable(res[, c("model", "prob", "signal", "rr_median", "ip_placebo",
                     "ip_drug", "ip_difference")], digits = 3,
             row.names = FALSE)
```

| model   |  prob | signal | rr_median | ip_placebo | ip_drug | ip_difference |
|:--------|------:|:-------|----------:|-----------:|--------:|--------------:|
| BDRIBS1 | 0.668 | FALSE  |     1.241 |      0.311 |   0.359 |         0.048 |
| BDRIBS2 | 0.948 | FALSE  |     2.012 |      0.308 |   0.514 |         0.205 |

The last three columns are the posterior means of the cumulative
incidences at week 52 and of their difference; the true values are 0.30,
0.45 and 0.15. The posterior means of the placebo incidence are 0.311
with the total exposure time and 0.308 with the total time at risk.

## Operating characteristics

The folder `reproduce` of the installed package contains the simulations
of Figures 2 to 10 of Duncan (2026) with 1,000 simulated trials per
scenario. For the medium (Figure 3) and high (Figure 4) incidences, the
proportions of trials with a signal (power) of the two models are

| Figure | Placebo | Study drug |   N | Analysis (week) | Threshold | BDRIBS(1) | BDRIBS(2) |
|:-------|:--------|:-----------|----:|----------------:|:----------|----------:|----------:|
| 3      | 10%     | 10%        | 350 |            69.3 | 0.925     |     0.000 |     0.000 |
| 3      | 10%     | 15%        | 350 |            69.3 | 0.925     |     0.002 |     0.016 |
| 3      | 10%     | 20%        | 350 |            69.3 | 0.925     |     0.039 |     0.151 |
| 3      | 10%     | 25%        | 350 |            69.3 | 0.925     |     0.381 |     0.659 |
| 4      | 30%     | 30%        | 200 |            78.0 | 0.95      |     0.000 |     0.000 |
| 4      | 30%     | 40%        | 200 |            78.0 | 0.95      |     0.000 |     0.009 |
| 4      | 30%     | 45%        | 200 |            78.0 | 0.95      |     0.000 |     0.068 |
| 4      | 30%     | 50%        | 200 |            78.0 | 0.95      |     0.000 |     0.182 |

The incidences are at week 52. BDRIBS(1) has much lower power than
BDRIBS(2), as stated by Duncan (2026, Section 4.1); the corresponding
statements of the article are judged PASS (Figure 3) and PASS (Figure 4)
in `inst/reproduce/reproduce_Duncan2026_figures.R`. The article does not
report all settings of the simulations, so the published powers are not
reproduced exactly; `results/SUMMARY.md` in the folder `reproduce` gives
the details.

## Numerical integration and comparison with MCMC

The posterior computation of BDRIBS(2) is that of BDRIBS(1) with another
exposure. The numerical integration is described in
[`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md),
and the comparison with MCMC is repeated here.

The bdribs package of Mukhopadhyay et al. (2018) samples the same model
by MCMC with JAGS. The script
`inst/validation/validate_bdribs_package.R` compares
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with [`bdribs::bdribs()`](https://rdrr.io/pkg/bdribs/man/bdribs.html)
for the settings below, in which E is the exposure given to both
functions. The results in this section are read from the output of that
script (`validate_bdribs_package.csv` in the folder `validation/results`
of the installed package); JAGS is not run when this article is built.

| Data | y | E | k | Historical events | Historical patient-years | Background rate | Prior of p |
|:---|---:|---:|---:|---:|---:|:---|:---|
| B1 | 5 | 500 | 2 | 18 | 4000 | fixed (inf.type = 1) | Beta(1, 1) |
| B2 | 5 | 500 | 2 | 18 | 4000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B3 | 5 | 500 | 2 | 18 | 4000 | fixed (inf.type = 1) | Beta(0.5, 0.5) |
| B4 | 9 | 540 | 2 | 33 | 3000 | fixed (inf.type = 1) | Beta(1, 1) |
| B5 | 9 | 540 | 2 | 33 | 3000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B6 | 7 | 800 | 1 | 18 | 4000 | Gamma prior (inf.type = 2) | Beta(1, 1) |
| B7 | 30 | 2000 | 1 | 20 | 2000 | Gamma prior (inf.type = 2) | Beta(0.5, 0.5) |

[`bdribs::bdribs()`](https://rdrr.io/pkg/bdribs/man/bdribs.html) was run
with 2 chains, 5,000 burn-in iterations and 100,000 iterations. The
Monte Carlo standard error (SE) of each MCMC estimate uses the effective
sample size of the draws
([`coda::effectiveSize()`](https://rdrr.io/pkg/coda/man/effectiveSize.html));
for the posterior probabilities, the effective sample sizes are between
46,649 and 115,901. For the median, the table shows the MCMC proportion
of draws below the posterior median computed by bssd, which should be
0.5. The column Ratio is the absolute difference divided by the SE.

| Data | Quantity             |   bssd |   MCMC | Difference |      SE | Ratio |
|:-----|:---------------------|-------:|-------:|-----------:|--------:|------:|
| B1   | P(r \> 1)            | 0.8154 | 0.8163 |   -8.9e-04 | 1.6e-03 |  0.55 |
| B1   | P(r \> 1.2)          | 0.7523 | 0.7526 |   -2.8e-04 | 1.7e-03 |  0.16 |
| B1   | P(r \> 1.5)          | 0.6520 | 0.6518 |    1.9e-04 | 1.8e-03 |  0.11 |
| B1   | P(r \<= bssd median) | 0.5000 | 0.5006 |   -5.5e-04 | 1.7e-03 |  0.33 |
| B1   | mean of r            | 2.1834 | 2.1842 |   -7.7e-04 | 4.4e-03 |  0.17 |
| B2   | P(r \> 1)            | 0.7850 | 0.7859 |   -9.3e-04 | 1.7e-03 |  0.55 |
| B2   | P(r \> 1.2)          | 0.7195 | 0.7202 |   -7.4e-04 | 1.8e-03 |  0.41 |
| B2   | P(r \> 1.5)          | 0.6196 | 0.6209 |   -1.3e-03 | 1.8e-03 |  0.69 |
| B2   | P(r \<= bssd median) | 0.5000 | 0.4985 |    1.5e-03 | 1.8e-03 |  0.85 |
| B2   | mean of r            | 2.1907 | 2.1929 |   -2.2e-03 | 5.5e-03 |  0.40 |
| B3   | P(r \> 1)            | 0.8461 | 0.8465 |   -4.6e-04 | 1.6e-03 |  0.28 |
| B3   | P(r \> 1.2)          | 0.7933 | 0.7930 |    3.0e-04 | 1.8e-03 |  0.17 |
| B3   | P(r \> 1.5)          | 0.7056 | 0.7052 |    3.7e-04 | 1.9e-03 |  0.20 |
| B3   | P(r \<= bssd median) | 0.5000 | 0.5005 |   -5.4e-04 | 1.7e-03 |  0.31 |
| B3   | mean of r            | 2.3919 | 2.3897 |    2.2e-03 | 5.0e-03 |  0.44 |
| B4   | P(r \> 1)            | 0.7530 | 0.7550 |   -2.0e-03 | 1.5e-03 |  1.30 |
| B4   | P(r \> 1.2)          | 0.6392 | 0.6411 |   -1.9e-03 | 1.6e-03 |  1.21 |
| B4   | P(r \> 1.5)          | 0.4647 | 0.4670 |   -2.3e-03 | 1.5e-03 |  1.57 |
| B4   | P(r \<= bssd median) | 0.5000 | 0.4975 |    2.5e-03 | 1.5e-03 |  1.67 |
| B4   | mean of r            | 1.5218 | 1.5259 |   -4.1e-03 | 2.2e-03 |  1.84 |
| B5   | P(r \> 1)            | 0.7192 | 0.7181 |    1.1e-03 | 1.7e-03 |  0.66 |
| B5   | P(r \> 1.2)          | 0.6089 | 0.6072 |    1.8e-03 | 1.7e-03 |  1.04 |
| B5   | P(r \> 1.5)          | 0.4485 | 0.4464 |    2.0e-03 | 1.6e-03 |  1.27 |
| B5   | P(r \<= bssd median) | 0.5000 | 0.5015 |   -1.5e-03 | 1.7e-03 |  0.87 |
| B5   | mean of r            | 1.5231 | 1.5205 |    2.6e-03 | 2.9e-03 |  0.90 |
| B6   | P(r \> 1)            | 0.8152 | 0.8130 |    2.2e-03 | 1.5e-03 |  1.46 |
| B6   | P(r \> 1.2)          | 0.7580 | 0.7558 |    2.3e-03 | 1.6e-03 |  1.40 |
| B6   | P(r \> 1.5)          | 0.6679 | 0.6644 |    3.6e-03 | 1.7e-03 |  2.05 |
| B6   | P(r \<= bssd median) | 0.5000 | 0.5028 |   -2.8e-03 | 1.7e-03 |  1.61 |
| B6   | mean of r            | 2.3795 | 2.3646 |    1.5e-02 | 5.8e-03 |  2.57 |
| B7   | P(r \> 1)            | 0.8777 | 0.8752 |    2.5e-03 | 1.5e-03 |  1.68 |
| B7   | P(r \> 1.2)          | 0.7999 | 0.7971 |    2.9e-03 | 1.8e-03 |  1.57 |
| B7   | P(r \> 1.5)          | 0.6582 | 0.6548 |    3.5e-03 | 2.2e-03 |  1.58 |
| B7   | P(r \<= bssd median) | 0.5000 | 0.5017 |   -1.7e-03 | 2.3e-03 |  0.74 |
| B7   | mean of r            | 1.9489 | 1.9450 |    3.9e-03 | 4.8e-03 |  0.82 |

In all 35 items, the bssd value is within 2.6 standard errors of the
MCMC value, and 33 items are within two standard errors. The largest
absolute difference in a posterior probability is 0.0036, while the SEs
of the MCMC probabilities are between 0.0015 and 0.0023. Differences of
this size are expected from the Monte Carlo error of MCMC alone, so the
difference between the numerical integration of bssd and MCMC is
negligible in practice. The validation script judges an item PASS when
the difference is at most four standard errors; 35 of the 35 items are
PASS.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227-244.
<https://doi.org/10.1080/10543406.2024.2395532>

Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
detection of potential risk using inference on blinded safety data.
*Pharmaceutical Statistics*, 17(6), 823-834.
<https://doi.org/10.1002/pst.1898>
