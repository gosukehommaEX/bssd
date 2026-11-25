# Bayesian safety signal detection with bssd

``` r

library(bssd)
```

## Introduction

The safety of a study drug is often monitored while the treatment
assignment of an ongoing trial is still blinded. The blinded data show
the adverse events of special interest, but not the treatment groups of
the subjects. bssd computes the posterior probability that the risk of
the event is higher on the study drug than on placebo (control) from
such data, with published Bayesian methods. Several of these methods
were computed by Markov chain Monte Carlo (MCMC) sampling. bssd computes
all posterior distributions by numerical integration, so the results do
not depend on a random seed, and they are fast enough for simulation
studies of the operating characteristics
([`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)).

## Models and functions

| Model | Function | Blinded data | Tested quantity | Article |
|----|----|----|----|----|
| BDRIBS(1) (Mukhopadhyay et al., 2018; Duncan, 2026) | [`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md) | number of subjects with the event and total exposure time | rate ratio | [`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md) |
| BDRIBS(2) (Mukhopadhyay et al., 2018; Duncan, 2026) | [`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md) | number of subjects with the event and total time at risk | rate ratio | [`vignette("bdribs2")`](https://gosukehommaex.github.io/bssd/articles/bdribs2.md) |
| BINBETA (Ball and Schnell, 2016) | [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md) | numbers of subjects and of subjects with the event | study drug proportion or ratio | [`vignette("binbeta")`](https://gosukehommaex.github.io/bssd/articles/binbeta.md) |
| MULTIBIN (Goren et al., 2020) | [`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md) | the same, for each of several studies | common rate ratio | [`vignette("binbeta")`](https://gosukehommaex.github.io/bssd/articles/binbeta.md) |
| CANARY (Gould and Wang, 2017) | [`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md) | total event count, with external control data | ratio, difference or odds ratio | [`vignette("canary")`](https://gosukehommaex.github.io/bssd/articles/canary.md) |
| EXP_BERN (Duncan, 2026) | [`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md) | event indicator and follow-up time of each subject | rate ratio | [`vignette("expbern")`](https://gosukehommaex.github.io/bssd/articles/expbern.md) |
| LL_BERN (Duncan, 2026) | [`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md) | event indicator and follow-up time of each subject | ratio of the cumulative incidences at `t_star` | [`vignette("llbern")`](https://gosukehommaex.github.io/bssd/articles/llbern.md) |
| EXPTIME (Schnell and Ball, 2016) | [`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md) | number of events and follow-up time of each subject | study drug rate or rate ratio | [`vignette("exptime")`](https://gosukehommaex.github.io/bssd/articles/exptime.md) |
| MAPBERN (Brock et al., 2023) | [`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md), [`map_beta()`](https://gosukehommaex.github.io/bssd/reference/map_beta.md) | event indicator and time at risk of each subject, historical control arms | ratio of the event probabilities | [`vignette("mapbern")`](https://gosukehommaex.github.io/bssd/articles/mapbern.md) |
| BDRIBS monitoring (Waterhouse et al., 2022) | [`bdribs_grid()`](https://gosukehommaex.github.io/bssd/reference/bdribs_grid.md), [`pooled_ratio()`](https://gosukehommaex.github.io/bssd/reference/pooled_ratio.md) | as BDRIBS, for grids of data and pooled trials | rate ratio | [`vignette("bdribs-monitoring")`](https://gosukehommaex.github.io/bssd/articles/bdribs-monitoring.md) |

Each article describes the model and its numerical integration, with
examples and comparisons with MCMC or with the published results.

The study drug and placebo are randomized in the ratio \\k:1\\. The
models differ in what they assume about the treatment assignment. In
EXP_BERN, LL_BERN, BINBETA, MULTIBIN, EXPTIME and MAPBERN, each subject
is on study drug with probability \\k / (k + 1)\\ (simple
randomization); in BDRIBS the exposure is divided between the arms in
the ratio \\k:1\\, and in CANARY the arm sizes are fixed at their
expected values (Gould and Wang, 2017). The prior of the placebo
parameters must be informative in all models (for example, from
historical placebo data), because the blinded data alone do not separate
the two arms.

A signal is declared when the posterior probability that the tested
quantity exceeds its critical value is larger than a threshold, for
example \\P(RR \> rr_0 \mid \mathrm{data}) \> \mathrm{cutoff}\\ for a
relative risk \\RR\\. The defaults are \\rr_0 = 1\\ and
\\\mathrm{cutoff} = 0.95\\. BINBETA, EXPTIME and MAPBERN can instead
test the study drug proportion or rate against a critical value.

[`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
gives the smallest event count with a signal for a number of subjects
(the boundaries of Ball and Schnell, 2016, and the control charts of
Gould and Wang, 2017),
[`count_oc()`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)
the exact probability that the count reaches it, and
[`elicit_prior()`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md)
the parameters of a prior from two quantiles.

## Example with aggregate data

Hypothetical example 2 of Mukhopadhyay et al. (2018, Section 4.4) has a
2:1 randomization, a background event rate of 0.011 per patient-year on
placebo, and 540 patient-years of blinded exposure. The BDRIBS model
sets the prior of the rate ratio \\r\\ through the probability \\p = r k
/ (r k + 1)\\ that a subject with the event is on the study drug, with
\\p \sim \mathrm{Beta}(a, b)\\ (`p_prior`;
[`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md)).
The screening analysis of the article treats the background rate as
known, which corresponds to `pt_hp = Inf`. With the uniform prior for
\\p\\, the posterior probabilities for 8 and 9 events are

``` r

ex2 <- expand.grid(rr0 = c(1, 1.3), y = c(8, 9))
ex2$bssd <- mapply(function(y, rr0) {
  bdribs_post(y = y, exposure = 540, k = 2, lambda_hp = 0.011, pt_hp = Inf,
              p_prior = c(1, 1), rr0 = rr0, details = FALSE)$prob
}, ex2$y, ex2$rr0)
ex2$published <- c(0.622, 0.434, 0.756, 0.584)
knitr::kable(ex2[, c("y", "rr0", "bssd", "published")], digits = 3,
             col.names = c("Events", "rr0", "bssd", "Mukhopadhyay et al."))
```

| Events | rr0 |  bssd | Mukhopadhyay et al. |
|-------:|----:|------:|--------------------:|
|      8 | 1.0 | 0.619 |               0.622 |
|      8 | 1.3 | 0.433 |               0.434 |
|      9 | 1.0 | 0.753 |               0.756 |
|      9 | 1.3 | 0.580 |               0.584 |

The published values were computed by MCMC with the bdribs package, and
they differ from the values of bssd by at most 0.004. The print method
summarizes the posterior with 9 events:

``` r

fit9 <- bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011,
                    pt_hp = Inf, p_prior = c(1, 1))
fit9
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.753
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.522, median 1.438, 95% interval (0.379, 3.143)
fit90 <- bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011,
                     pt_hp = Inf, p_prior = c(1, 1), level = 0.90)
```

The posterior mean of \\r\\ is 1.52, and the article reports 1.53. The
article gives (0.52, 2.83) as the 95% interval of \\r\\. This agrees
with the equal-tailed 90% interval, (0.51, 2.82), which `level = 0.90`
gives, and the equal-tailed 95% interval is (0.38, 3.14).

## Example with subject-level data

The following data are simulated with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)
under the placebo distribution and the largest study drug distribution
of scenario 6 in Table 1 of Duncan (2026). The 150 subjects are
randomized 1:1, they enter the trial uniformly over 26 weeks, and each
subject is followed for up to 104 weeks. The time to the event has a
log-logistic distribution with \\\sigma = 0.8\\, for which the hazard
first increases and then decreases. The cumulative incidence at week 104
is 3.5% on placebo and 18.1% on the study drug. The blinded analysis is
at week 117, the third interim analysis in the simulations of the folder
`reproduce` of the installed package.

``` r

set.seed(2026)
d <- sim_blinded(n = 150, k = 1, accrual = 26, follow_up = 104,
                 analysis_time = 117, dist = "loglogistic",
                 par_p = c(7.3, 0.8), par_d = c(5.85, 0.8))
head(d[, c("id", "entry", "follow_up", "event", "event_time")])
#>   id     entry follow_up event event_time
#> 1  1 14.446685  102.5533     0         NA
#> 2  2  3.150173  104.0000     0         NA
#> 3  3  3.683752  104.0000     0         NA
#> 4  4  3.806313  104.0000     0         NA
#> 5  5  9.020966  104.0000     0         NA
#> 6  6  1.833814  104.0000     0         NA
tar <- time_at_risk(d$follow_up, d$event, d$event_time)
tar
#>     exposure time_at_risk      n_event 
#>     15096.61     14605.76        14.00
```

The column `group` of `d` holds the true treatment group, which the
blinded analyses do not use.

The data are analyzed with the four models of Duncan (2026). Their
priors of the placebo parameters are based on historical placebo data
with the event rate \\\lambda\_{HP}\\ and the total patient time
\\PT\_{HP}\\, discounted by the factor \\DF\\. The priors follow the
simulations in the folder `reproduce`: a historical placebo group of 75
subjects, the size of the current placebo group, followed for 104 weeks
and discounted by \\DF = 2\\. For BDRIBS and EXP_BERN, \\\lambda\_{HP}\\
is the constant hazard per week that gives the placebo incidence at week
104, and \\PT\_{HP}\\ is in patient-weeks. For LL_BERN, the Gamma priors
of \\\mu_p\\ and \\\sigma_p\\ have the means `mu_hp` and `sigma_hp`, and
`pt_hp` divided by `df_p` sets their precision. The simulations use the
historical patient time in patient-years for it. For LL_BERN, the
relative risk is the ratio of the cumulative incidences at the time
\\t^\* = 104\\; Duncan (2026) does not state which relative risk is used
for this model, and
[`?llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
describes this choice and the priors.

``` r

lambda_hp <- -log(1 - ll_cdf(104, mu = 7.3, sigma = 0.8)) / 104
pt_weeks <- 75 * 104
fits <- list(
  BDRIBS1 = bdribs_post(y = tar[["n_event"]], exposure = tar[["exposure"]],
                        k = 1, lambda_hp = lambda_hp, pt_hp = pt_weeks,
                        df = 2, t_star = 104),
  BDRIBS2 = bdribs_post(y = tar[["n_event"]],
                        exposure = tar[["time_at_risk"]], k = 1,
                        lambda_hp = lambda_hp, pt_hp = pt_weeks, df = 2,
                        t_star = 104),
  EXP_BERN = expbern_post(y = d$event, t = d$follow_up, k = 1,
                          lambda_hp = lambda_hp, pt_hp = pt_weeks, df = 2,
                          t_star = 104),
  LL_BERN = llbern_post(y = d$event, t = d$follow_up, k = 1, mu_hp = 7.3,
                        sigma_hp = 0.8, pt_hp = pt_weeks / 52, df_p = 2,
                        t_star = 104)
)
res <- do.call(rbind, lapply(fits, summary))
res$model <- names(fits)
knitr::kable(res[, c("model", "prob", "signal", "rr_median", "ip_placebo",
                     "ip_drug", "ip_difference")], digits = 3,
             row.names = FALSE)
```

| model    |  prob | signal | rr_median | ip_placebo | ip_drug | ip_difference |
|:---------|------:|:-------|----------:|-----------:|--------:|--------------:|
| BDRIBS1  | 0.791 | FALSE  |     2.879 |      0.053 |   0.118 |         0.065 |
| BDRIBS2  | 0.805 | FALSE  |     3.031 |      0.053 |   0.124 |         0.071 |
| EXP_BERN | 0.820 | FALSE  |     3.243 |      0.052 |   0.131 |         0.079 |
| LL_BERN  | 0.935 | FALSE  |     4.217 |      0.041 |   0.152 |         0.111 |

``` r

sig <- res$model[res$signal]
```

The column `prob` is \\P(RR \> 1 \mid \mathrm{data})\\. Models with a
signal at the threshold 0.95: none. The relative risk is the rate ratio
\\r\\ for the BDRIBS and EXP_BERN models and \\RR(104)\\ for LL_BERN.
The table shows its posterior median, because the posterior mean of
\\r\\ is infinite for EXP_BERN under the default Beta(0.5, 0.5) prior of
\\p\\ (see
[`?expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)).
The last three columns are the posterior means of the cumulative
incidences at week 104 and of their difference. The true values are
0.035, 0.181 and 0.146.

## Accuracy of the numerical integration

[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
integrates over one dimension with
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html). The
other functions use Gauss-Legendre and Chebyshev quadrature on a region
located by a coarse grid
([`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
in four dimensions,
[`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)
and the Poisson likelihood of
[`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
in one, and the others in two). With `control$check = TRUE`, the
default, they repeat the computation on a finer grid, return the finer
result, and report an estimate of the integration error of `prob`, for
example for the data above:

``` r

vapply(fits[c("EXP_BERN", "LL_BERN")], function(f) f$integration$error,
       numeric(1))
#>     EXP_BERN      LL_BERN 
#> 1.814905e-07 3.733315e-04
```

The number of nodes can be changed with `control`.
[`?llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
explains why the error estimate of LL_BERN is twice the difference
between the two grids.
[`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
turns the check off by default to save time.

## Operating characteristics

[`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
simulates trials with
[`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md),
applies the selected models at each analysis, and returns the proportion
of trials with a signal and the posterior mean of the difference in the
cumulative incidence. The example below uses the BDRIBS models with 200
simulated trials to keep the run time short; the column `se` is the
Monte Carlo standard error of `prob_signal`. The cumulative incidence at
week 52 is 10% on placebo and 20% on the study drug, 350 subjects enter
over 26 weeks, and the analyses are at weeks 26, 52 and 78.

``` r

rate_of <- function(ip, t) -log(1 - ip) / t
design <- list(n = 350, k = 1, accrual = 26, follow_up = 52,
               analysis_time = c(26, 52, 78), dist = "exponential",
               par_p = rate_of(0.10, 52), par_d = rate_of(0.20, 52))
oc <- sim_oc(n_sim = 200, design = design,
             models = c("BDRIBS1", "BDRIBS2"),
             priors = list(bdribs = list(lambda_hp = rate_of(0.10, 52),
                                         pt_hp = 175 * 52, df = 2)),
             t_star = 52, seed = 1)
knitr::kable(oc$summary[, c("analysis_time", "model", "prob_signal", "se",
                            "mean_ip_difference",
                            "mean_observed_difference")],
             digits = 3, row.names = FALSE)
```

| analysis_time | model | prob_signal | se | mean_ip_difference | mean_observed_difference |
|---:|:---|---:|---:|---:|---:|
| 26 | BDRIBS1 | 0.095 | 0.021 | 0.060 | 0.026 |
| 26 | BDRIBS2 | 0.115 | 0.023 | 0.067 | 0.026 |
| 52 | BDRIBS1 | 0.050 | 0.015 | 0.055 | 0.071 |
| 52 | BDRIBS2 | 0.130 | 0.024 | 0.073 | 0.071 |
| 78 | BDRIBS1 | 0.045 | 0.015 | 0.054 | 0.093 |
| 78 | BDRIBS2 | 0.120 | 0.023 | 0.076 | 0.093 |

The true difference in the cumulative incidence at week 52 is 0.1. The
column `mean_observed_difference` is the mean difference between the
proportions of subjects with the event in the two groups at the
analysis, computed from the true group assignment. The scripts in the
folder `reproduce` (`system.file("reproduce", package = "bssd")`) run
the simulations of Duncan (2026) for all four models with 1000 trials
per scenario.

## Reproduction and validation

The installed package contains the scripts and the results of the
reproduction of published results (folder `reproduce`) and of the
comparison with independent computations (folder `validation`). Each
item is judged as PASS (agreement within the tolerance), EXPLAINED (a
difference whose reason is registered), FAIL (a difference without a
registered reason) or INFO (no value to compare with). The numbers of
items are

| file                                | PASS | EXPLAINED | FAIL | INFO |
|:------------------------------------|-----:|----------:|-----:|-----:|
| reproduce_Duncan2026.csv            |   15 |         6 |    0 |    7 |
| reproduce_Duncan2026_prior.csv      |   23 |         6 |    0 |    0 |
| reproduce_Duncan2026_simulation.csv |   14 |         4 |    0 |  113 |
| reproduce_Goren2020.csv             |   12 |         0 |    0 |    0 |
| reproduce_Gould2017.csv             |   71 |        17 |    0 |    0 |
| reproduce_Mukhopadhyay2018.csv      |   82 |         4 |    0 |    1 |
| reproduce_Schnell2016.csv           |   27 |         3 |    1 |   16 |
| reproduce_Waterhouse2022.csv        |   11 |         0 |    1 |    2 |
| validate_bdribs_package.csv         |   35 |         0 |    0 |    0 |
| validate_convergence.csv            |   19 |         0 |    0 |    0 |
| validate_gould2017_author_grid.csv  |   28 |         0 |    0 |    0 |
| validate_jags.csv                   |   24 |         0 |    0 |    3 |
| validate_jags_methods.csv           |   16 |         0 |    0 |    1 |
| validate_methods_python.csv         |  174 |         0 |    0 |    0 |
| validate_sim_oc.csv                 |   29 |         0 |    0 |    0 |

The reproduction covers Table 1, the hypothetical examples and the
operating characteristics of Mukhopadhyay et al. (2018), Table 1, Figure
1 and the simulations of Figures 2 to 10 of Duncan (2026), and the
published results of Schnell and Ball (2016), Gould and Wang (2017),
Brock et al. (2023), Goren et al. (2020) and Waterhouse et al. (2022).
Many simulation results are INFO, because the articles do not report all
settings of the simulations and some published values are read from the
figures. The validation compares
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with the bdribs package; the other posterior functions with MCMC in JAGS
(for
[`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md),
with the Poisson likelihood) and with independent Python computations
(an exact finite mixture for
[`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md));
the default integration grids with finer grids; and
[`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
with an independent Python simulation. The number of FAIL items is 2.
The files `README.md` and `results/SUMMARY.md` in each folder describe
the scripts, the tolerances and the reasons for the EXPLAINED items.

## References

Ball, G. and Schnell, P. M. (2016). Blinded safety signal monitoring for
the FDA IND reporting final rule. In J. Lin et al. (eds.), *Statistical
Applications from Clinical Trials and Personalized Medicine to Finance
and Business Analytics*, ICSA Book Series in Statistics, 201-209.
Springer, Cham. <https://doi.org/10.1007/978-3-319-42568-9_17>

Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
Penard, N. (2023). A Bayesian method for safety signal detection in
ongoing blinded randomised controlled trials. *Pharmaceutical
Statistics*, 22(2), 378-395. <https://doi.org/10.1002/pst.2278>

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227-244.
<https://doi.org/10.1080/10543406.2024.2395532>

Goren, E., Lin, L.-A., Blaustein, R. O. and Ball, G. (2020). Bayesian
meta-analysis of safety outcomes using blinded clinical trial data.
*Therapeutic Innovation & Regulatory Science*, 54, 1557-1565.
<https://doi.org/10.1007/s43441-020-00183-7>

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92-104.
<https://doi.org/10.1002/sim.7129>

Mukhopadhyay, S., Waterhouse, B., and Hartford, A. (2018). Bayesian
detection of potential risk using inference on blinded safety data.
*Pharmaceutical Statistics*, 17(6), 823-834.
<https://doi.org/10.1002/pst.1898>

Schnell, P. M. and Ball, G. (2016). A Bayesian exposure-time method for
clinical trial safety monitoring with blinded data. *Therapeutic
Innovation & Regulatory Science*, 50(6), 833-838.
<https://doi.org/10.1177/2168479016656702>

Waterhouse, B., Hartford, A., Mukhopadhyay, S., Ferguson, R. and
Hendrickson, B. A. (2022). Using the Bayesian detection of potential
risk using inference on blinded safety data (BDRIBS) method to support
the decision to refer an event for unblinded evaluation. *Pharmaceutical
Statistics*, 21(2), 372-385. <https://doi.org/10.1002/pst.2175>
