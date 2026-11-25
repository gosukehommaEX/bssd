# Randomization ratio of blinded data pooled from several studies

Computes the effective randomization ratio of study drug to control for
blinded data pooled from studies with different randomization ratios, as
described by Waterhouse et al. (2022, Appendix A). The pooled numbers of
events and the pooled exposure can then be analyzed with
[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)
and this ratio as `k`.

## Usage

``` r
pooled_ratio(exposure, k)
```

## Arguments

- exposure:

  Exposures of the studies (positive).

- k:

  Randomization ratios of the studies (study drug to control).

## Value

The pooled randomization ratio.

## Details

Study \\j\\ with exposure \\E_j\\ and randomization ratio \\k_j\\
contributes the expected exposures \\E_j k_j / (k_j + 1)\\ on study drug
and \\E_j / (k_j + 1)\\ on control. The pooled ratio is the ratio of the
totals, \$\$k = \sum_j E_j k_j / (k_j + 1) \Big/ \sum_j E_j / (k_j +
1).\$\$

## References

Waterhouse, B., Hartford, A., Mukhopadhyay, S., Ferguson, R. and
Hendrickson, B. A. (2022). Using the Bayesian detection of potential
risk using inference on blinded safety data (BDRIBS) method to support
the decision to refer an event for unblinded evaluation. *Pharmaceutical
Statistics*, 21(2), 372–385.
[doi:10.1002/pst.2175](https://doi.org/10.1002/pst.2175)

## See also

[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md),
[`bdribs_grid`](https://gosukehommaex.github.io/bssd/reference/bdribs_grid.md)

## Examples

``` r
# Example 3 of Waterhouse et al. (2022): 250 patient-years at 1:1 and
# 80 patient-years at 2:1
k <- pooled_ratio(c(250, 80), c(1, 2))
k
#> [1] 1.175824
bdribs_post(y = 16, exposure = 330, k = k, lambda_hp = 0.031, pt_hp = Inf)
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.9101
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.905, median 1.846, 95% interval (0.6727, 3.469)
```
