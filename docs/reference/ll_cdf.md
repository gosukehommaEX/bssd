# Cumulative incidence of the log-logistic distribution

Computes the cumulative incidence (distribution function) of the
log-logistic time-to-event distribution in the parameterization used by
Duncan (2026), which is the parameterization of the log-logistic model
of the SAS LIFEREG procedure: \$\$F(t) = 1 - \frac{1}{1 + \exp\[(\log
t - \mu) / \sigma\]}.\$\$

## Usage

``` r
ll_cdf(t, mu, sigma)
```

## Arguments

- t:

  Numeric vector of times (positive).

- mu:

  Location (intercept) parameter on the log-time scale.

- sigma:

  Scale parameter (positive). The hazard is decreasing when `sigma >= 1`
  and first increases and then decreases when `sigma < 1`.

## Value

Numeric vector of cumulative incidences. The arguments are recycled to a
common length.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`ll_mu`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md)

## Examples

``` r
# incidence at week 104 under the placebo assumption of scenario 6
ll_cdf(104, mu = 7.3, sigma = 0.8)
#> [1] 0.03490787
```
