# Prior parameters from two quantiles

Finds the parameters of a Gamma, Beta or lognormal distribution with two
given quantiles, or of a lognormal distribution with a given mean and a
given probability content of an interval.

## Usage

``` r
elicit_prior(q, p, family = c("gamma", "beta", "lnorm"), mean = NULL)
```

## Arguments

- q:

  Two increasing quantiles, or the interval when `mean` is given.

- p:

  Two increasing probabilities, or the probability content of the
  interval when `mean` is given.

- family:

  One of `"gamma"`, `"beta"` and `"lnorm"`.

- mean:

  Optional mean of the lognormal distribution.

## Value

A named numeric vector: `c(shape, rate)` for the Gamma distribution,
`c(shape1, shape2)` for the Beta distribution and `c(meanlog, sdlog)`
for the lognormal distribution.

## Details

With `q = c(q1, q2)` and `p = c(p1, p2)`, the parameters solve \\F(q_1)
= p_1\\ and \\F(q_2) = p_2\\. For the Gamma distribution, the rate is
\\G^{-1}(p_1; s) / q_1\\ for shape \\s\\, and the shape solves the
second equation; for the Beta distribution, the second shape parameter
solves the first equation for each value of the first one, which then
solves the second equation; for the lognormal distribution the solution
is explicit. The equations are solved by
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), so the quantiles are
matched to about 10 significant digits. Schnell and Ball (2016) and
Gould and Wang (2017, Table I) elicit their priors in this way; their
own software minimizes the squared differences of the probabilities
instead.

For `family = "lnorm"` with `mean` given, `q` is an interval
`c(lower, upper)` and `p` is its probability content: the lognormal
distribution has mean `mean` and probability `p` between `lower` and
`upper`. This is the prior of the risk ratio of Gould and Wang (2017).

## References

Gould, A. L. and Wang, W. B. (2017). Monitoring potential adverse event
rate differences using data from blinded trials: the canary in the coal
mine. *Statistics in Medicine*, 36(1), 92–104.
[doi:10.1002/sim.7129](https://doi.org/10.1002/sim.7129)

Schnell, P. M. and Ball, G. (2016). A Bayesian exposure-time method for
clinical trial safety monitoring with blinded data. *Therapeutic
Innovation & Regulatory Science*, 50(6), 833–838.
[doi:10.1177/2168479016656702](https://doi.org/10.1177/2168479016656702)

## Examples

``` r
# 50 percent sure that the rate is below 0.5 and 75 percent sure that it
# is below 2 (Schnell and Ball, 2016)
elicit_prior(c(0.5, 2), c(0.5, 0.75), "gamma")
#>     shape      rate 
#> 0.3460050 0.2089832 
# median 0.02 and 95th percentile 0.04 of a proportion
elicit_prior(c(0.02, 0.04), c(0.5, 0.95), "beta")
#>     shape1     shape2 
#>   4.549595 207.162323 
# lognormal prior of the risk ratio of Gould and Wang (2017)
elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)
#>    meanlog      sdlog 
#> -0.5475868  1.0465054 
```
