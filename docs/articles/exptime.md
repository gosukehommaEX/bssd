# EXPTIME: exposure-time Poisson model

``` r

library(bssd)
```

## Model

Schnell and Ball (2016) monitor the rate of an adverse event of special
interest with blinded data: for each enrolled subject, the follow-up
time \\t_i\\ and the number of events \\y_i\\. With the unknown
treatment assignment \\a_i\\ (1 for study drug), \\ Y_i \mid a_i \sim
\mathrm{Poisson}(\lambda\_{a_i} t_i), \quad a_i \sim
\mathrm{Bernoulli}(k / (k + 1)), \quad \lambda_c \sim
\mathrm{Gamma}(\alpha_c, \beta_c), \quad \lambda_d \sim
\mathrm{Gamma}(\alpha_d, \beta_d), \\ where \\k\\ is the randomization
ratio of study drug to control. The prior of the control rate is
informative (from historical control data), and the prior of the study
drug rate is weaker. A safety signal is raised when the posterior
probability that \\\lambda_d\\ exceeds a critical rate is larger than a
threshold.

Summing over the assignment of each subject gives the likelihood \\
\prod_i \left\[ \frac{1}{k + 1} \lambda_c^{y_i} e^{-\lambda_c t_i} +
\frac{k}{k + 1} \lambda_d^{y_i} e^{-\lambda_d t_i} \right\], \\ so the
posterior is two-dimensional. Schnell and Ball (2016) sample it by Gibbs
sampling;
[`exptime_post()`](https://gosukehommaex.github.io/bssd/reference/exptime_post.md)
integrates it numerically. With event indicators as `y` and the times at
risk (up to the first event) as `t`, the model is the exponential model
for the time to the first event, which is how
[`sim_oc()`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md)
uses it (model `"EXPTIME"`).

## Example

The priors of the case study of Schnell and Ball (2016): 339 events in
26,064 person-months of historical control data, and for the study drug
the same data discounted to one hundredth plus 4 events in 167
person-months from a phase 2 study. The critical rate is 20% above the
control rate 0.013, and the threshold is 0.62.

``` r

set.seed(1)
n <- 600
t <- runif(n, 2, 18)
drug <- rbinom(n, 1, 2 / 3)
y <- rpois(n, ifelse(drug == 1, 0.016, 0.011) * t)
c(subjects = n, events = sum(y), person_months = round(sum(t)))
#>      subjects        events person_months 
#>           600            92          6019
```

``` r

fit <- exptime_post(y, t, k = 2, prior_c = c(339, 26064),
                    prior_d = c(339 / 100 + 4, 26064 / 100 + 167),
                    rate0 = 1.2 * 0.013, cutoff = 0.62)
fit
#> Model: EXPTIME 
#> Tested quantity: event rate lambda_d of study drug
#> P(lambda_d > 0.0156 | data) = 0.6252
#> Signal (posterior probability > 0.62): yes
#> Posterior of r: mean 1.262, median 1.253, 95% interval (0.9105, 1.667)
#> Posterior of lambda_d: mean 0.01638, median 0.0163, 95% interval (0.01223, 0.02097)
```

The posterior probability of the rate ratio is also returned:

``` r

fit$probs
#>      drug     ratio 
#> 0.6251912 0.9207968
```

Priors can be elicited from two quantiles with
[`elicit_prior()`](https://gosukehommaex.github.io/bssd/reference/elicit_prior.md);
for example, 50% sure that a rate is below 0.5 and 75% sure that it is
below 2:

``` r

elicit_prior(c(0.5, 2), c(0.5, 0.75), "gamma")
#>     shape      rate 
#> 0.3460050 0.2089832
```

## Numerical integration

The posterior is integrated over \\\log \lambda_c\\ by Gauss-Legendre
quadrature. For each of its nodes, the interval of \\\lambda_d\\ that
holds the conditional posterior mass is located on the scale \\z =
\Phi^{-1}(F_d(\lambda_d))\\, where \\F_d\\ is the prior distribution
function of \\\lambda_d\\, and the density is interpolated at Chebyshev
points in \\z\\ or in \\F_d(\lambda_d)\\, whichever represents it more
accurately. The first scale keeps the far tails of an informative prior,
which matters when the data put \\\lambda_d\\ many prior standard
deviations above its prior mean; the second keeps a vague prior from
creating long tails. The posterior probability of a threshold that
depends on \\\lambda_c\\ (such as the rate ratio) is obtained from the
Chebyshev series of each row. Compared with independent Python
computations on dense grids split exactly at the thresholds
(`inst/validation/validate_methods_python.R`, three data sets including
one with a study drug rate far in the upper tail of the prior of Schnell
and Ball, 2016), the largest absolute difference of the 5 posterior
probabilities is 2.9e-09, and all 32 compared quantities are PASS.

## Comparison with MCMC

A JAGS model with the latent treatment assignments
(`inst/validation/validate_jags_methods.R`, 4 chains of 50,000 draws)
gives the following results for a data set of 60 subjects.

| quantity                                        |   MCMC |   bssd | judgement |
|:------------------------------------------------|-------:|-------:|:----------|
| P(lambda_d \> 0.03 \| data)                     | 0.3926 | 0.3912 | PASS      |
| P(r \> 1.2 \| data)                             | 0.6383 | 0.6385 | PASS      |
| posterior mean of lambda_c                      | 0.0204 | 0.0204 | PASS      |
| MCMC proportion below the bssd posterior median | 0.5000 | 0.5001 | PASS      |

## Reproduction of Schnell and Ball (2016)

`inst/reproduce/reproduce_Schnell2016.R` generates the simulated case
study with the data generation of the authors’ code and seed 1, and
compares Table 2 and the prior elicitation example: 27 PASS, 3
EXPLAINED, 1 FAIL and 16 INFO. The posterior summaries in Table 2 are
Monte Carlo estimates of the authors’ Gibbs sampler, whose Metropolis
step uses Poisson probabilities of the arm totals; bssd computes the
posterior of the subject-level model above. The reasons for the items
that are not PASS are in `inst/reproduce/results/SUMMARY.md`.

The authors’ code was run outside the package (the code has no stated
license): with seed 1 it reproduces the published probabilities, and
with 200,000 draws its sampler and a sampler with the subject-level
likelihood give the following probabilities that the study drug rate
exceeds 0.0156.

| month | published | authors_code | authors_sampler | subject_level_sampler |  bssd |
|------:|----------:|-------------:|----------------:|----------------------:|------:|
|     4 |     0.710 |        0.710 |           0.679 |                 0.720 | 0.721 |
|     8 |     0.518 |        0.518 |           0.498 |                 0.601 | 0.603 |
|    12 |     0.872 |        0.872 |           0.839 |                 0.931 | 0.932 |
|    16 |     0.968 |        0.968 |           0.973 |                 0.996 | 0.996 |

The simulation of Table 1
(`inst/reproduce/reproduce_Schnell2016_simulation.R`) has not been run
for this installation.

## References

Schnell, P. M. and Ball, G. (2016). A Bayesian exposure-time method for
clinical trial safety monitoring with blinded data. *Therapeutic
Innovation & Regulatory Science*, 50(6), 833–838.
<https://doi.org/10.1177/2168479016656702>
