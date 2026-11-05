# Operating characteristics of the signal detection models by simulation

Simulates blinded trials with
[`sim_blinded`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md),
applies the selected models at each analysis, and summarizes the
probability of declaring a signal (power, or the false positive rate
when the study drug and placebo have the same distribution) and the
posterior mean of the difference in cumulative incidence.

## Usage

``` r
sim_oc(
  n_sim,
  design,
  models = c("BDRIBS1", "BDRIBS2", "EXP_BERN", "LL_BERN"),
  priors,
  rr0 = 1,
  cutoff = 0.95,
  t_star = NULL,
  control = list(),
  seed = NULL,
  keep = FALSE
)
```

## Arguments

- n_sim:

  Number of simulated trials.

- design:

  List of arguments of
  [`sim_blinded`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md).

- models:

  Character vector of models to apply.

- priors:

  List of prior settings: `bdribs` (arguments `lambda_hp`, `pt_hp` and
  optionally `df` and `p_prior` of
  [`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
  used by `"BDRIBS1"` and `"BDRIBS2"`), `expbern` (the same arguments
  for
  [`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md);
  `bdribs` is used when it is not given) and `llbern` (arguments
  `mu_hp`, `sigma_hp`, `pt_hp` and optionally `df_p`, `mu_d_prior` and
  `sigma_d_prior` of
  [`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)).

- rr0:

  Critical value of the relative risk.

- cutoff:

  Posterior probability threshold for a signal.

- t_star:

  Time at which the relative risk of `"LL_BERN"` and the cumulative
  incidences are evaluated. Required for `"LL_BERN"`.

- control:

  List with elements `EXP_BERN` and `LL_BERN`, the `control` arguments
  of
  [`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
  and
  [`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md).
  By default the convergence check is turned off.

- seed:

  Optional seed passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html).

- keep:

  Logical; if `TRUE`, the results of each simulated trial are returned.

## Value

A list with elements

- `summary`:

  Data frame with one row per analysis and model and the columns `look`,
  `analysis_time`, `model`, `n_ok` (successful fits), `n_failed`,
  `prob_signal` (proportion of trials with a signal), `se` (its Monte
  Carlo standard error), `mean_prob` (mean posterior probability),
  `mean_ip_difference` (mean of the posterior means of the difference in
  cumulative incidence at `t_star`), `mean_observed_difference` (mean of
  the difference between the proportions of subjects with the event in
  the study drug and placebo groups at the analysis, computed from the
  true group assignment), `mean_n` and `mean_events`.

- `truth`:

  True cumulative incidences at `t_star` for placebo and study drug and
  their difference.

- `results`:

  Data frame of the results of each trial (only when `keep = TRUE`).

## Details

The models are `"BDRIBS1"`
([`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with the total exposure time), `"BDRIBS2"`
([`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
with the total time at risk), `"EXP_BERN"`
([`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md))
and `"LL_BERN"`
([`llbern_post`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)).
The posterior functions are called with `details = FALSE`. A fit that
fails is recorded as `NA` and counted in `n_failed`, and a warning with
the first error message is given.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`sim_blinded`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)

## Examples

``` r
design <- list(n = 100, k = 1, accrual = 0, follow_up = 52,
               dist = "exponential", par_p = -log(0.9) / 52,
               par_d = -log(0.75) / 52)
oc <- sim_oc(n_sim = 20, design = design, models = c("BDRIBS1", "BDRIBS2"),
             priors = list(bdribs = list(lambda_hp = -log(0.9) / 52,
                                         pt_hp = 2000)),
             t_star = 52, seed = 1)
oc$summary
#>   look analysis_time   model n_ok n_failed prob_signal         se mean_prob
#> 1    1            52 BDRIBS1   20        0         0.0 0.00000000 0.6903407
#> 2    1            52 BDRIBS2   20        0         0.1 0.06708204 0.7447935
#>   mean_ip_difference mean_observed_difference mean_n mean_events
#> 1         0.06858147                    0.157    100       17.75
#> 2         0.09781049                    0.157    100       17.75
```
