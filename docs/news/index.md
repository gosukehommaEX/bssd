# Changelog

## bssd 0.2.0

- New functions for other published methods for blinded safety data, all
  computed by numerical integration and returning class `bssd_post`:
  - [`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md):
    the exposure-time Poisson model of Schnell and Ball (2016) with
    subject-level event counts and follow-up times.
  - [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md):
    the binomial model of Ball and Schnell (2016) with the numbers of
    enrolled subjects and of subjects with the event.
  - [`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md):
    the method of Gould and Wang (2017) for the total event count, with
    the ratio, the difference or the odds ratio as the metric and a
    binomial or Poisson likelihood.
  - [`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)
    and
    [`map_beta()`](https://gosukehommaex.github.io/bssd/reference/map_beta.md):
    the meta-analytic predictive prior from historical control arms and
    the mixture likelihood of Brock et al. (2023).
  - [`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md):
    the meta-analysis of several blinded studies of Goren et al. (2020)
    with a common rate ratio.
- [`bdribs_grid()`](https://gosukehommaex.github.io/bssd/reference/bdribs_grid.md)
  and
  [`pooled_ratio()`](https://gosukehommaex.github.io/bssd/reference/pooled_ratio.md):
  the contour plots and the pooled randomization ratio of Waterhouse et
  al. (2022) for BDRIBS.
- [`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
  gives the smallest event count with a signal for a number of subjects
  ([`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
  or
  [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)),
  [`count_oc()`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)
  the exact probability that the event count reaches it, and
  [`elicit_prior()`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md)
  the parameters of a Gamma, Beta or lognormal prior from two quantiles.
- [`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
  also runs the models `"EXPTIME"` and `"MAPBERN"`.
- The two-dimensional integrals of
  [`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md),
  [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md),
  [`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)
  and
  [`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
  (binomial likelihood) locate, for each node of the control parameter,
  the interval of the second coordinate that holds the conditional
  posterior mass, and
  [`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)
  locates the posterior of the log rate ratio on its standardized scale.
  The results stay accurate when the data conflict with an informative
  prior of the study drug parameter, and, in
  [`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md),
  when the event probabilities are large (all subjects with the event,
  for example).
- [`print()`](https://rdrr.io/r/base/print.html) shows the tested
  quantity of the new models and, for
  [`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md),
  [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)
  and
  [`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md),
  the posterior of the study drug rate or proportion. Summaries that are
  not available, such as the posterior mean of
  [`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md),
  are not shown.
- The README,
  [`vignette("bssd")`](https://gosukehommaex.github.io/bssd/articles/bssd.md),
  the reference index and the articles of the website list all models
  together.
- Vignettes
  [`vignette("exptime")`](https://gosukehommaex.github.io/bssd/articles/exptime.md),
  [`vignette("binbeta")`](https://gosukehommaex.github.io/bssd/articles/binbeta.md),
  [`vignette("canary")`](https://gosukehommaex.github.io/bssd/articles/canary.md),
  [`vignette("mapbern")`](https://gosukehommaex.github.io/bssd/articles/mapbern.md)
  and
  [`vignette("bdribs-monitoring")`](https://gosukehommaex.github.io/bssd/articles/bdribs-monitoring.md);
  [`vignette("bssd")`](https://gosukehommaex.github.io/bssd/articles/bssd.md)
  lists the methods, and
  [`vignette("llbern")`](https://gosukehommaex.github.io/bssd/articles/llbern.md)
  explains that the lower limit of the credible interval of the relative
  risk can be unstable.
- `inst/reproduce` reproduces results of Schnell and Ball (2016), Gould
  and Wang (2017), Brock et al. (2023), Goren et al. (2020) and
  Waterhouse et al. (2022), and `inst/validation` compares the new
  functions with independent Python computations and with MCMC in JAGS.

## bssd 0.1.0

- First version.
- [`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
  computes the posterior distribution of the rate ratio under the BDRIBS
  model of Mukhopadhyay et al. (2018) as specified by Duncan (2026),
  with the total exposure time (BDRIBS(1)) or the total time at risk
  (BDRIBS(2)) as `exposure`. The background rate has a Gamma prior from
  historical placebo data (`lambda_hp`, `pt_hp`, `df`), or it is known
  when `pt_hp = Inf` (the screening analysis of Mukhopadhyay et al.,
  2018).
- [`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
  and
  [`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
  compute the posterior distributions under the exponential Bernoulli
  (EXP_BERN) and the log-logistic Bernoulli (LL_BERN) mixture models of
  Duncan (2026) from subject-level data (`y`, `t`). For LL_BERN, the
  relative risk is the ratio of the cumulative incidences at `t_star`,
  and the priors of the study drug parameters can be set with
  `mu_d_prior` and `sigma_d_prior`.
- All posterior quantities are computed by numerical integration instead
  of Markov chain Monte Carlo sampling.
  [`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
  and
  [`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
  repeat the computation on a finer grid and report an estimate of the
  integration error when `control$check` is `TRUE` (the default).
- [`time_at_risk()`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md)
  computes the total exposure time and the total time at risk from
  subject-level data.
  [`ll_cdf()`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)
  and
  [`ll_mu()`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md)
  give the cumulative incidence and the location parameter of the
  log-logistic distribution.
- [`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)
  simulates blinded trial data with exponential or log-logistic times to
  the event, and
  [`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
  estimates the operating characteristics of the models by simulation.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for the
  results (class `bssd_post`).
- Vignettes: an overview
  ([`vignette("bssd")`](https://gosukehommaex.github.io/bssd/articles/bssd.md))
  and one article for each model
  ([`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md),
  [`vignette("bdribs2")`](https://gosukehommaex.github.io/bssd/articles/bdribs2.md),
  [`vignette("expbern")`](https://gosukehommaex.github.io/bssd/articles/expbern.md)
  and
  [`vignette("llbern")`](https://gosukehommaex.github.io/bssd/articles/llbern.md))
  with examples, the accuracy of the numerical integration and the
  comparison with MCMC.
- The folder `inst/reproduce` reproduces the numerical results of
  Mukhopadhyay et al. (2018) and Duncan (2026), and `inst/validation`
  compares the package with the bdribs package, with MCMC in JAGS and
  with independent Python implementations.
