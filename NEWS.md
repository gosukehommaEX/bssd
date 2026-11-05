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
* The folder `inst/reproduce` reproduces the numerical results of
  Mukhopadhyay et al. (2018) and Duncan (2026), and `inst/validation`
  compares the package with the bdribs package, with MCMC in JAGS and with
  independent Python implementations.
