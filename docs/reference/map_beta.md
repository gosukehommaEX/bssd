# Meta-analytic predictive prior of an event probability

Derives the meta-analytic predictive (MAP) prior of the event
probability per unit time in a new control arm from historical control
arms with the hierarchical Beta model of Brock et al. (2023), and
returns it as a mixture of Beta distributions for
[`mapbern_post`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md).

## Usage

``` r
map_beta(
  events = NULL,
  n = NULL,
  y = NULL,
  t = NULL,
  study = NULL,
  a_prior = c(1/3, 1/3),
  b_prior = c(5, 5),
  tol = 1e-10,
  control = list()
)
```

## Arguments

- events:

  Numbers of subjects with the event in the historical arms (summary
  data with one unit of time at risk per subject).

- n:

  Numbers of subjects in the historical arms (summary data).

- y, t, study:

  Subject-level data instead of `events` and `n`: event indicators,
  times at risk, and arm identifiers.

- a_prior:

  Shape parameters of the Beta prior of \\a\\; the default is that of
  Brock et al. (2023).

- b_prior:

  Shape parameters of the Beta prior of \\b\\; the default is that of
  Brock et al. (2023).

- tol:

  Total weight of the components that may be dropped.

- control:

  List of settings of the numerical integration: `n_coarse` (default
  41), `n_ab` (Gauss-Legendre nodes per coordinate, default 32), `n_pi`
  (nodes per panel for subject-level data, default 16), `eps` (default
  1e-8), `max_iter` (default 12) and `check` (default `TRUE`).

## Value

A data frame with columns `w`, `shape1` and `shape2`, one row per
component of the mixture, ordered by decreasing weight. Its attributes
are `n_total` (\\N\\), `hyper` (posterior means of \\a\\ and \\b\\),
`moments` (mean and standard deviation of the MAP prior) and `error`
(with `control$check = TRUE`, the largest absolute difference of the MAP
prior distribution function at its 2.5, 50 and 97.5 percent points
between the default and a finer grid).

## Details

Historical arm \\j\\ of \\J\\ has \\N_j\\ subjects, and \\N = \sum_j
N_j\\. With the event probability \\\pi_j\\ per unit time and a constant
hazard, a subject with time at risk \\t\\ has the event with probability
\\1 - (1 - \pi_j)^t\\. The model is \$\$\pi_j \mid a, b \sim
\mathrm{Beta}(a b N, (1 - a) b N), \quad a \sim \mathrm{Beta}(a_1, a_2),
\quad b \sim \mathrm{Beta}(b_1, b_2),\$\$ so \\a\\ is the mean event
probability and \\b\\ discounts the historical information. The MAP
prior of the probability \\\pi^\*\\ in a new arm is \$\$p(\pi^\*) =
\int\int \mathrm{Beta}(\pi^\*; a b N, (1 - a) b N) \\ p(a, b \mid
\mathrm{data}) \\ da \\ db.\$\$ The posterior of \\(a, b)\\ is
integrated by Gauss-Legendre quadrature on \\(\mathrm{logit}(a),
\mathrm{logit}(b))\\ after a coarse grid has located the region of the
posterior mass; each quadrature node gives one Beta component of the
mixture, with the quadrature weight times the posterior density as its
weight. The components with the smallest weights, which together have a
weight of at most `tol`, are dropped. Brock et al. (2023) sample the
same model by MCMC and approximate the predictive distribution by a
mixture fitted to the draws.

The likelihood of an arm is the Beta-binomial probability when only the
numbers of events (`events`, `n`) are given, which assumes one unit of
time at risk for every subject. With subject-level data (`y`, `t`,
`study`), \\\pi_j\\ is integrated out by Gauss-Legendre quadrature on
\\\mathrm{logit}(\pi_j)\\.

## References

Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
Penard, N. (2023). A Bayesian method for safety signal detection in
ongoing blinded randomised controlled trials. *Pharmaceutical
Statistics*, 22(2), 378–395.
[doi:10.1002/pst.2278](https://doi.org/10.1002/pst.2278)

## See also

[`mapbern_post`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)

## Examples

``` r
prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
attr(prior, "moments")
#>       mean         sd 
#> 0.10334757 0.03403532 
nrow(prior)
#> [1] 1422
```
