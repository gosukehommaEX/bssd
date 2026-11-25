# MAPBERN: meta-analytic predictive prior and mixture likelihood

``` r

library(bssd)
```

## Model

Brock et al. (2023) estimate the posterior probability that the event
probability per unit time on study drug, \\\pi_d\\, exceeds that on
control, \\\pi_c\\, in an ongoing blinded study. Each subject with time
at risk \\t_i\\ has the event with probability \\1 - (1 - \pi)^{t_i}\\
on the arm received, and the likelihood of the blinded data is the
mixture over the two arms weighted by the randomization probabilities.
The prior of \\\pi_c\\ is the meta-analytic predictive (MAP) prior from
historical control arms, \\ \pi_j \mid a, b \sim \mathrm{Beta}(a b N,
(1 - a) b N), \quad a \sim \mathrm{Beta}(1/3, 1/3), \quad b \sim
\mathrm{Beta}(5, 5), \\ where \\N\\ is the total number of historical
subjects, \\a\\ is the mean event probability and \\b\\ discounts the
historical information. The prior of \\\pi_d\\ is Beta(1/3, 1/3).

[`map_beta()`](https://gosukehommaex.github.io/bssd/reference/map_beta.md)
integrates the posterior of \\(a, b)\\ by quadrature, which gives the
MAP prior as a mixture of Beta distributions, and
[`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md)
integrates the two-dimensional posterior of \\(\pi_c, \pi_d)\\. Brock et
al. (2023) use Stan for both steps and fit a mixture to the MCMC draws
of the MAP prior. The likelihood is that of the exponential Bernoulli
mixture model
([`expbern_post()`](https://gosukehommaex.github.io/bssd/reference/expbern_post.md))
with hazards \\-\log(1 - \pi)\\; only the priors differ.

## Example

``` r

prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
attr(prior, "moments")
#>       mean         sd 
#> 0.10334757 0.03403532
nrow(prior)
#> [1] 1422
head(prior, 3)
#>            w   shape1   shape2
#> 1 0.04203266 17.25062 148.2004
#> 2 0.04030330 14.91639 128.1470
#> 3 0.03702683 15.25908 150.1920
```

``` r

set.seed(2)
n <- 200
drug <- rbinom(n, 1, 0.75)
y <- rbinom(n, 1, ifelse(drug == 1, 0.18, 0.095))
fit <- mapbern_post(y, t = rep(1, n), k = 3, prior_c = prior)
fit
#> Model: MAPBERN 
#> Tested quantity: ratio rr = pi_d / pi_c
#> P(rr > 1 | data) = 0.8602
#> Signal (posterior probability > 0.95): no
#> Posterior of rr: mean 1.84, median 1.627, 95% interval (0.6569, 4.263)
#> Posterior of pi_d: mean 0.166, median 0.1648, 95% interval (0.09953, 0.2396)
```

A prior derived with the RBesT package can be used as well: the rows
`w`, `a` and `b` of a `mixbeta` object are the columns `w`, `shape1` and
`shape2` of `prior_c`.

``` r

# mix: a mixbeta object of the RBesT package
prior_c <- data.frame(w = mix["w", ], shape1 = mix["a", ],
                      shape2 = mix["b", ])
```

## Numerical integration and comparison with MCMC

Against independent Python computations (dense composite Gauss-Legendre
rules split exactly at \\\pi_d = rr_0 \pi_c\\ for
[`mapbern_post()`](https://gosukehommaex.github.io/bssd/reference/mapbern_post.md),
and a dense grid over \\(a, b)\\ for
[`map_beta()`](https://gosukehommaex.github.io/bssd/reference/map_beta.md)),
the largest absolute difference of the 43 compared quantities is
1.1e-06. A JAGS model with the mixture prior and latent treatment
assignments gives:

| quantity                                        |   MCMC |   bssd | judgement |
|:------------------------------------------------|-------:|-------:|:----------|
| P(pi_d \> pi_c \| data)                         | 0.9574 | 0.9570 | PASS      |
| posterior mean of pi_c                          | 0.1011 | 0.1014 | PASS      |
| MCMC proportion below the bssd posterior median | 0.5000 | 0.4976 | PASS      |

## Reproduction of Brock et al. (2023)

The simulation of Tables 1 and 3
(`inst/reproduce/reproduce_Brock2023_simulation.R`) has not been run for
this installation.

## References

Brock, K., Chen, C., Ho, S., Fuller, G., Woolfolk, J., McShea, C. and
Penard, N. (2023). A Bayesian method for safety signal detection in
ongoing blinded randomised controlled trials. *Pharmaceutical
Statistics*, 22(2), 378–395. <https://doi.org/10.1002/pst.2278>
