# Simulate blinded clinical trial data on an adverse event

Generates subject-level data of a randomized two-group trial at one or
more analysis times, as used in the simulations of Duncan (2026).

## Usage

``` r
sim_blinded(
  n,
  k = 1,
  accrual,
  follow_up,
  analysis_time = accrual + follow_up,
  dist = c("exponential", "loglogistic"),
  par_p,
  par_d,
  dropout = 0
)
```

## Arguments

- n:

  Total number of subjects.

- k:

  Randomization ratio of study drug to placebo (`k:1`).

- accrual:

  Length of the accrual period (0 if all subjects enter at time 0).

- follow_up:

  Planned follow-up time of each subject.

- analysis_time:

  Calendar times (from the start of accrual) of the blinded analyses.
  The default is the time at which all subjects have completed
  follow-up.

- dist:

  Distribution of the time to the event: `"exponential"` or
  `"loglogistic"`.

- par_p, par_d:

  Parameters of the distribution for placebo and study drug: the hazard
  rate for `"exponential"`, or `c(mu, sigma)` for `"loglogistic"`.

- dropout:

  Hazard rate of dropout (0 for no dropout).

## Value

A data frame with one row per subject and analysis, and the columns
`look` (index of the analysis), `analysis_time`, `id`, `group`
(`"placebo"` or `"drug"`; not used by the blinded analyses), `entry`,
`follow_up`, `event` (1 if the event occurred during follow-up),
`event_time` (time from entry to the event, `NA` without the event) and
`time_at_risk` (`event_time` for subjects with the event, `follow_up`
otherwise).

## Details

Subjects enter the trial uniformly over `[0, accrual]`, and
`round(n * k / (k + 1))` of them, chosen at random, are assigned to the
study drug. The time to the first occurrence of the event follows the
exponential distribution (`dist = "exponential"`; `par_p` and `par_d`
are the hazard rates) or the log-logistic distribution of
[`ll_cdf`](https://gosukehommaex.github.io/bssd/reference/ll_cdf.md)
(`dist = "loglogistic"`; `par_p` and `par_d` are `c(mu, sigma)`). When
`dropout > 0`, the time to dropout follows the exponential distribution
with rate `dropout`.

At the analysis time \\\tau\\, the follow-up time of a subject who
entered at time \\e\\ is the smallest of `follow_up`, \\\tau - e\\ and
the time to dropout, and the subject has the event if the time to the
event is not larger than the follow-up time. Subjects who have not
entered by \\\tau\\ are not included. The same subjects and event times
are used at all analysis times.

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`sim_oc`](https://gosukehommaex.github.io/bssd/reference/sim_oc.md),
[`ll_mu`](https://gosukehommaex.github.io/bssd/reference/ll_mu.md)

## Examples

``` r
set.seed(1)
# 52-week incidence of 10% (placebo) and 20% (study drug)
d <- sim_blinded(n = 200, k = 1, accrual = 26, follow_up = 52,
                 analysis_time = c(39, 78), dist = "exponential",
                 par_p = -log(0.9) / 52, par_d = -log(0.8) / 52)
table(d$look, d$event)
#>    
#>       0   1
#>   1 183  17
#>   2 169  31
```
