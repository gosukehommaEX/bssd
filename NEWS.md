# bssd 0.2.0

* New functions for other published methods for blinded safety data, all
  computed by numerical integration and returning class `bssd_post`:
  * `exptime_post()`: the exposure-time Poisson model of Schnell and
    Ball (2016) with subject-level event counts and follow-up times.
  * `binbeta_post()`: the binomial model of Ball and Schnell (2016) with the
    numbers of enrolled subjects and of subjects with the event.
  * `canary_post()`: the method of Gould and Wang (2017) for the total event
    count, with the ratio, the difference or the odds ratio as the metric
    and a binomial or Poisson likelihood.
  * `mapbern_post()` and `map_beta()`: the meta-analytic predictive prior
    from historical control arms and the mixture likelihood of Brock et al.
    (2023).
  * `multibin_post()`: the meta-analysis of several blinded studies of Goren
    et al. (2020) with a common rate ratio.
* `bdribs_grid()` and `pooled_ratio()`: the contour plots and the pooled
  randomization ratio of Waterhouse et al. (2022) for BDRIBS.
* `signal_count()` gives the smallest event count with a signal for a number
  of subjects (`canary_post()` or `binbeta_post()`), `count_oc()` the exact
  probability that the event count reaches it, and `elicit_prior()` the
  parameters of a Gamma, Beta or lognormal prior from two quantiles.
* `sim_oc()` also runs the models `"EXPTIME"` and `"MAPBERN"`.
* `print()` shows the tested quantity of the new models and, for
  `exptime_post()`, `binbeta_post()` and `mapbern_post()`, the posterior of
  the study drug rate or proportion.
* Vignettes `vignette("exptime")`, `vignette("binbeta")`,
  `vignette("canary")`, `vignette("mapbern")` and
  `vignette("bdribs-monitoring")`; `vignette("bssd")` lists the methods, and
  `vignette("llbern")` explains that the lower limit of the credible
  interval of the relative risk can be unstable.
* `inst/reproduce` reproduces results of Schnell and Ball (2016), Gould and
  Wang (2017), Brock et al. (2023), Goren et al. (2020) and Waterhouse et
  al. (2022), and `inst/validation` compares the new functions with
  independent Python computations and with MCMC in JAGS.

# bssd 0.1.0

* First version.
* `bdribs_post()` computes the posterior distribution of the rate ratio under
  the BDRIBS model of Mukhopadhyay et al. (2018) as specified by Duncan
  (2026), with the total exposure time (BDRIBS(1)) or the total time at risk
  (BDRIBS(2)) as `exposure`. The background rate has a Gamma prior from
  historical placebo data (`lambda_hp`, `pt_hp`, `df`), or it is known when
  `pt_hp = Inf` (the screening analysis of Mukhopadhyay et al., 2018).
* `expbern_post()` and `llbern_post()` compute the posterior distributions
  under the exponential Bernoulli (EXP_BERN) and the log-logistic Bernoulli
  (LL_BERN) mixture models of Duncan (2026) from subject-level data (`y`,
  `t`). For LL_BERN, the relative risk is the ratio of the cumulative
  incidences at `t_star`, and the priors of the study drug parameters can be
  set with `mu_d_prior` and `sigma_d_prior`.
* All posterior quantities are computed by numerical integration instead of
  Markov chain Monte Carlo sampling. `expbern_post()` and `llbern_post()`
  repeat the computation on a finer grid and report an estimate of the
  integration error when `control$check` is `TRUE` (the default).
* `time_at_risk()` computes the total exposure time and the total time at
  risk from subject-level data. `ll_cdf()` and `ll_mu()` give the cumulative
  incidence and the location parameter of the log-logistic distribution.
* `sim_blinded()` simulates blinded trial data with exponential or
  log-logistic times to the event, and `sim_oc()` estimates the operating
  characteristics of the models by simulation.
* `print()` and `summary()` methods for the results (class `bssd_post`).
* Vignettes: an overview (`vignette("bssd")`) and one article for each model
  (`vignette("bdribs1")`, `vignette("bdribs2")`, `vignette("expbern")` and
  `vignette("llbern")`) with examples, the accuracy of the numerical
  integration and the comparison with MCMC.
* The folder `inst/reproduce` reproduces the numerical results of
  Mukhopadhyay et al. (2018) and Duncan (2026), and `inst/validation`
  compares the package with the bdribs package, with MCMC in JAGS and with
  independent Python implementations.
