# Monitoring with BDRIBS: contour plots, sensitivity and pooled trials

``` r

library(bssd)
```

Waterhouse et al. (2022) describe how a safety management team uses
BDRIBS
([`vignette("bdribs1")`](https://gosukehommaex.github.io/bssd/articles/bdribs1.md))
in ongoing blinded trials: a contour plot of the posterior probability
\\P(RR \> 1)\\ over the numbers of events and the patient-years of
exposure prepared before the first patient enters the trial, a
sensitivity analysis over the background rate, and pooling of blinded
data from several trials with different randomization ratios. The
examples use a fixed background rate and the Jeffreys prior Beta(0.5,
0.5) for \\p\\.

## Contour plot

Example 1 of Waterhouse et al. (2022): MACE with a background rate of
3.1 events per 100 patient-years and randomization 1:1.

``` r

g <- bdribs_grid(y = 4:10, exposure = seq(40, 240, by = 10),
                 lambda_hp = 0.031, p_prior = c(0.5, 0.5))
z <- matrix(g$prob, nrow = 7)
contour(4:10, seq(40, 240, by = 10), z, levels = c(0.6, 0.7, 0.8, 0.9),
        xlab = "Number of events", ylab = "Patient-years of exposure")
```

![](bdribs-monitoring_files/figure-html/grid-1.png)

The probabilities for six and seven events (Figures 3 and 4):

``` r

subset(g, y %in% 6:7 & exposure %in% c(80, 100, 140))
#>    y exposure lambda_hp      prob
#> 31 6       80     0.031 0.9037555
#> 32 7       80     0.031 0.9649264
#> 45 6      100     0.031 0.8076505
#> 46 7      100     0.031 0.9128135
#> 73 6      140     0.031 0.5678571
#> 74 7      140     0.031 0.7324799
```

## Sensitivity to the background rate

With 10 events in 200 patient-years, the probabilities at the limits of
the 90% credible interval of the background rate from 496 events in
15,730 patient-years (Figures 5 and 6) are

``` r

bg <- qgamma(c(0.05, 0.95), 496, 15730)
bdribs_grid(y = 10, exposure = 200, lambda_hp = bg, p_prior = c(0.5, 0.5))
#>    y exposure  lambda_hp      prob
#> 1 10      200 0.02923996 0.8654598
#> 2 10      200 0.03389652 0.7554039
```

The background rate can also be given a Gamma prior instead of a fixed
value (`pt_hp` in
[`bdribs_post()`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)).

## Pooled trials

Example 3: trial A with 250 patient-years at 1:1 and trial B with 80
patient-years at 2:1 are pooled; the effective randomization ratio is
the ratio of the expected exposures on study drug and on control.

``` r

k <- pooled_ratio(c(250, 80), c(1, 2))
k
#> [1] 1.175824
bdribs_post(y = 16, exposure = 330, k = k, lambda_hp = 0.031, pt_hp = Inf,
            p_prior = c(0.5, 0.5))
#> Model: BDRIBS 
#> Relative risk: rate ratio r of study drug to placebo
#> P(r > 1 | data) = 0.9101
#> Signal (posterior probability > 0.95): no
#> Posterior of r: mean 1.905, median 1.846, 95% interval (0.6727, 3.469)
```

## Reproduction of Waterhouse et al. (2022)

`inst/reproduce/reproduce_Waterhouse2022.R` compares the probabilities,
exposures and the pooled randomization ratio given in the text and
figures (11 PASS, 1 FAIL). The published values were read from an R
Shiny application that samples by MCMC; the tolerance is 0.01 for
probabilities given in percent and 0.02 for those given as “about” (the
other tolerances are in the table).

| quantity | published | bssd | tolerance | judgement |
|:---|---:|---:|---:|:---|
| P(RR \> 1), 6 events, 80 patient-years | 0.90 | 0.904 | 0.020 | PASS |
| patient-years at which P(RR \> 1) falls to 0.80 with 6 events | 100.00 | 101.390 | 5.000 | PASS |
| patient-years at which P(RR \> 1) falls to 0.80 with 7 events | 125.00 | 126.954 | 5.000 | PASS |
| P(RR \> 1), 7 events, 140 patient-years | 0.74 | 0.732 | 0.010 | PASS |
| indicator P(RR \> 1) \> 0.80, 10 events, 200 patient-years | 1.00 | 1.000 | 0.500 | PASS |
| P(RR \> 1), 10 events, 200 PY, lower limit of the background rate | 0.86 | 0.865 | 0.020 | PASS |
| P(RR \> 1), 10 events, 200 PY, upper limit of the background rate | 0.76 | 0.755 | 0.020 | PASS |
| P(RR \> 1), 22 events, 504 patient-years, k = 2 | 0.70 | 0.697 | 0.020 | PASS |
| pooled randomization ratio | 1.18 | 1.176 | 0.005 | PASS |
| P(RR \> 1), 16 events, 330 PY, pooled ratio | 0.90 | 0.910 | 0.010 | FAIL |
| P(RR \> 1), 16 events, 330 PY, lower limit of the background rate | 0.85 | 0.847 | 0.020 | PASS |
| max \|bdribs_grid - bdribs_post\| on the grid of Figure 3 | 0.00 | 0.000 | 0.000 | PASS |

## References

Waterhouse, B., Hartford, A., Mukhopadhyay, S., Ferguson, R. and
Hendrickson, B. A. (2022). Using the Bayesian detection of potential
risk using inference on blinded safety data (BDRIBS) method to support
the decision to refer an event for unblinded evaluation. *Pharmaceutical
Statistics*, 21(2), 372–385. <https://doi.org/10.1002/pst.2175>
