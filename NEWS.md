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
* The two-dimensional integrals of `exptime_post()`, `binbeta_post()`,
  `mapbern_post()` and `canary_post()` (binomial likelihood) locate, for
  each node of the control parameter, the interval of the second coordinate
  that holds the conditional posterior mass, and `multibin_post()` locates
  the posterior of the log rate ratio on its standardized scale. The results
  stay accurate when the data conflict with an informative prior of the
  study drug parameter, and, in `canary_post()`, when the event
  probabilities are large (all subjects with the event, for example).
* `expbern_post()` integrates in the same way: for each Gauss-Legendre node
  of the log placebo rate, the interval of the second coordinate is located
  separately, and each row is interpolated in the prior distribution
  function of `p` or in its normal quantile, whichever gives the faster
  decay of the Chebyshev coefficients of its density. Before, one interval
  on the prior distribution function scale of `p` was used for all nodes,
  which lost accuracy when the data conflicted with an informative prior of
  `p` (cases X1 to X3 of `inst/validation/validate_methods_python.R`).
* `mapbern_post()` integrates over the normal quantile of an envelope of
  the prior of `pi_c` when the prior has components with a shape parameter
  below 1 (total weight above `control$eps`), as the MAP prior from few
  historical events has. Before, the integration over the logit of `pi_c`
  missed the posterior when the logit densities of such components were nearly
  flat over a wide range, and `prob` could be wrong with a small error
  estimate (cases MS1 to MS7 of `inst/validation/validate_methods_python.R`).
  `binbeta_post()` and `exptime_post()` integrate over the normal quantile
  of the prior of the control parameter when its shape parameter is below 1.
  The coordinate is shown in `integration$coordinate`; with other priors
  the results are unchanged.
* The posterior mean of the ratio of `exptime_post()`, `binbeta_post()` and
  `mapbern_post()` is `Inf` when the shape parameter of the control prior
  (the smallest `shape1` of a mixture) is at most 1, because the posterior
  density of the control parameter near 0 is then proportional to that of
  the prior.
* The documentation of `llbern_post()` says that `control$n_outer` limits
  the accuracy when there are many events, and that the error estimate
  refers to `prob`.
* Faster computation with the same results up to rounding errors: the
  likelihoods of `expbern_post()` and `llbern_post()` are computed once for
  each group of subjects with the same event indicator and follow-up time
  (as in `exptime_post()` and `mapbern_post()`), and the C++ code of
  `llbern_post()` accumulates the product over the subjects in four parts.
  The coarse grid that locates the integration region of `expbern_post()`,
  `exptime_post()`, `binbeta_post()` and `mapbern_post()` computes the
  second coordinate once per column, the summaries of `llbern_post()`
  compute the cumulative incidence of the study drug once for each pair of
  its parameters, and the normal quantile coordinate of `mapbern_post()` and
  `binbeta_post()` needs fewer Newton steps. On Windows, where the
  exponential functions of the C library are slow, the C++ likelihoods use
  their own exp and expm1 (errors of about two units in the last place); on
  other systems they use those of the C library.
* `print()` shows the tested quantity of the new models and, for
  `exptime_post()`, `binbeta_post()` and `mapbern_post()`, the posterior of
  the study drug rate or proportion. Summaries that are not available, such
  as the posterior mean of `multibin_post()`, are not shown.
* The README, `vignette("bssd")`, the reference index and the articles of
  the website list all models together.
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
