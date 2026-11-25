# Package index

## Posterior probability of a safety signal

Posterior distributions under the BDRIBS, EXP_BERN and LL_BERN models,
and the methods for their results.

- [`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
  : Posterior probability of a safety signal under the BDRIBS model
- [`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md)
  : Posterior probability of a safety signal under the EXP_BERN model
- [`llbern_post()`](https://gosukehommaex.github.io/bssd/reference/llbern_post.md)
  : Posterior probability of a safety signal under the LL_BERN model
- [`print(`*`<bssd_post>`*`)`](https://gosukehommaex.github.io/bssd/reference/print.bssd_post.md)
  : Print a posterior summary for safety signal detection
- [`summary(`*`<bssd_post>`*`)`](https://gosukehommaex.github.io/bssd/reference/summary.bssd_post.md)
  : Summarize a posterior for safety signal detection as a data frame

## Further methods for blinded data

Posterior distributions under the methods of Schnell and Ball (2016),
Ball and Schnell (2016), Goren et al. (2020), Gould and Wang (2017) and
Brock et al. (2023).

- [`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md)
  : Posterior of the event rates under the exposure-time Poisson model
- [`binbeta_post()`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md)
  : Posterior of the event proportions under the blinded binomial model
- [`multibin_post()`](https://gosukehommaex.github.io/bssd/reference/multibin_post.md)
  : Posterior of the rate ratio from blinded event counts of several
  studies
- [`canary_post()`](https://gosukehommaex.github.io/bssd/reference/canary_post.md)
  : Posterior of the risk elevation from the total event count
- [`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)
  : Posterior probability of a higher event probability on study drug
  with a meta-analytic predictive prior
- [`map_beta()`](https://gosukehommaex.github.io/bssd/reference/map_beta.md)
  : Meta-analytic predictive prior of an event probability

## Monitoring tools

Grids of BDRIBS results, pooled trials, decision boundaries for event
counts and their exact operating characteristics, and elicitation of
priors.

- [`bdribs_grid()`](https://gosukehommaex.github.io/bssd/reference/bdribs_grid.md)
  : Posterior probabilities of the BDRIBS model on a grid of data
- [`pooled_ratio()`](https://gosukehommaex.github.io/bssd/reference/pooled_ratio.md)
  : Randomization ratio of blinded data pooled from several studies
- [`signal_count()`](https://gosukehommaex.github.io/bssd/reference/signal_count.md)
  : Smallest event count that gives a signal
- [`count_oc()`](https://gosukehommaex.github.io/bssd/reference/count_oc.md)
  : Probability that the blinded event count reaches a critical value
- [`elicit_prior()`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md)
  : Prior parameters from two quantiles

## Exposure and the log-logistic distribution

- [`time_at_risk()`](https://gosukehommaex.github.io/bssd/reference/time_at_risk.md)
  : Total exposure time and total time at risk
- [`ll_cdf()`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)
  : Cumulative incidence of the log-logistic distribution
- [`ll_mu()`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md) :
  Location parameter of the log-logistic distribution for a given
  incidence

## Simulation

- [`sim_blinded()`](https://gosukehommaex.github.io/bssd/reference/sim_blinded.md)
  : Simulate blinded clinical trial data on an adverse event
- [`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
  : Operating characteristics of the signal detection models by
  simulation
