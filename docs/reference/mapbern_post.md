# Posterior probability of a higher event probability on study drug with a meta-analytic predictive prior

Computes the posterior probability that the event probability per unit
time on study drug exceeds that on control, from blinded subject-level
data, with the mixture likelihood and the meta-analytic predictive (MAP)
prior of Brock et al. (2023).

## Usage

``` r
mapbern_post(
  y,
  t,
  k = 1,
  prior_c,
  prior_d = c(1/3, 1/3),
  rate0 = NULL,
  rr0 = 1,
  cutoff = 0.95,
  level = 0.95,
  details = TRUE,
  control = list()
)
```

## Arguments

- y:

  Event indicators of the subjects (1 if the subject had at least one
  event).

- t:

  Times at risk of the subjects (positive), in the unit of time of the
  probabilities.

- k:

  Randomization ratio of study drug to control (`k:1`); Brock et
  al. (2023) use `k = 3`.

- prior_c:

  Prior of the control probability: `c(shape1, shape2)` of a Beta
  distribution, or a data frame with columns `w`, `shape1` and `shape2`
  of a mixture of Beta distributions, such as the result of
  [`map_beta`](https://gosukehommaex.github.io/bssd/reference/map_beta.md).

- prior_d:

  Shape parameters of the Beta prior of the study drug probability.

- rate0:

  Optional critical value of the study drug probability; if given, it
  defines the signal instead of `rr0`.

- rr0:

  Critical value of the ratio \\\pi_d / \pi_c\\.

- cutoff:

  Posterior probability threshold for a signal.

- level:

  Probability of the equal-tailed credible intervals.

- details:

  Logical; if `FALSE`, the quantiles are not computed (they are `NA`).

- control:

  List of settings of the numerical integration, as in
  [`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md).

## Value

An object of class `bssd_post` as described in
[`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md),
with model `"MAPBERN"`, the probabilities `pi_c` and `pi_d`, and their
ratio `rr`. `prob` is \\\gamma\\.

## Details

Let \\\pi_c\\ and \\\pi_d\\ be the probabilities of at least one event
in one unit of time on control and study drug. With a constant hazard,
the probability of an event of subject \\i\\ with time at risk \\t_i\\
is \\1 - (1 - \pi)^{t_i}\\, and the likelihood of the blinded data is
the mixture \$\$\prod_i \[p_B(y_i \mid 1 - (1 - \pi_c)^{t_i}) /
(k + 1) + k p_B(y_i \mid 1 - (1 - \pi_d)^{t_i}) / (k + 1)\],\$\$ where
\\p_B\\ is the Bernoulli probability. The prior of \\\pi_c\\ is a
mixture of Beta distributions, typically the MAP prior from historical
control data
([`map_beta`](https://gosukehommaex.github.io/bssd/reference/map_beta.md));
the prior of \\\pi_d\\ is a Beta distribution, by default the vague
Beta(1/3, 1/3) used by Brock et al. (2023). The quantity of interest is
\\\gamma = P(\pi_d \> rr_0 \pi_c \mid \mathrm{data})\\; Brock et al.
(2023) use \\rr_0 = 1\\ and mention larger values.

The likelihood is that of the exponential Bernoulli mixture model
([`expbern_post`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md))
with hazards \\-\log(1 - \pi)\\; only the priors differ. The posterior
is integrated numerically on the coordinates \\\mathrm{logit}(\pi_c)\\
and the prior distribution function of \\\pi_d\\, as in
[`exptime_post`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md)
and
[`binbeta_post`](https://gosukehommaex.github.io/bssd/reference/binbeta_post.md).
A MAP prior from few historical events has components with a shape
parameter far below 1, whose densities of \\\mathrm{logit}(\pi_c)\\ are
nearly flat over a wide range. When the total weight of the components
with a shape parameter below 1 exceeds `control$eps`, the first
coordinate is therefore \\z_c = \Phi^{-1}(F_e(\pi_c))\\, where \\F_e\\
is the distribution function of an envelope of the prior: a mixture of
fewer Beta distributions, each with the smallest shape parameters of a
group of components. The prior density of \\z_c\\ is that of the
standard normal distribution times the ratio of the densities of the
prior and the envelope, which is bounded. When a component has \\shape1
\le 1\\, the posterior mean of the ratio \\\pi_d / \pi_c\\ is infinite
(`Inf`).

A prior derived with the RBesT package can be used: the rows `w`, `a`
and `b` of a `mixbeta` object are the columns `w`, `shape1` and `shape2`
of `prior_c`.

## References

Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
Penard, N. (2023). A Bayesian method for safety signal detection in
ongoing blinded randomised controlled trials. *Pharmaceutical
Statistics*, 22(2), 378–395.
[doi:10.1002/pst.2278](https://doi.org/10.1002/pst.2278)

## See also

[`map_beta`](https://gosukehommaex.github.io/bssd/reference/map_beta.md)

## Examples

``` r
# MAP prior from three historical control arms of 100 subjects each
prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
set.seed(2)
n <- 200
drug <- rbinom(n, 1, 0.75)
y <- rbinom(n, 1, ifelse(drug == 1, 0.18, 0.095))
mapbern_post(y, t = rep(1, n), k = 3, prior_c = prior)
#> Model: MAPBERN 
#> Tested quantity: ratio rr = pi_d / pi_c
#> P(rr > 1 | data) = 0.8602
#> Signal (posterior probability > 0.95): no
#> Posterior of rr: mean Inf, median 1.627, 95% interval (0.6569, 4.263)
#> Posterior of pi_d: mean 0.166, median 0.1648, 95% interval (0.09953, 0.2396)
```
