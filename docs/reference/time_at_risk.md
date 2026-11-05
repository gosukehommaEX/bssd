# Total exposure time and total time at risk

Computes the two exposure measures used by the BDRIBS models of Duncan
(2026) from subject-level data. BDRIBS(1) uses the total exposure time
(the sum of the follow-up times), and BDRIBS(2) uses the total time at
risk, in which the time after the first occurrence of the event is
removed.

## Usage

``` r
time_at_risk(follow_up, event, event_time)
```

## Arguments

- follow_up:

  Numeric vector of follow-up times (positive).

- event:

  Vector of event indicators (1 if the subject had at least one
  occurrence of the event, 0 otherwise).

- event_time:

  Numeric vector of times to the first occurrence of the event. Values
  for subjects without the event are not used and may be `NA`.

## Value

A named numeric vector with elements `exposure` (total exposure time),
`time_at_risk` (total time at risk), and `n_event` (number of subjects
with the event).

## References

Duncan, B. (2026). Robust safety monitoring and signal detection using
alternatives to the standard Poisson distribution. *Journal of
Biopharmaceutical Statistics*, 36(2), 227–244.
[doi:10.1080/10543406.2024.2395532](https://doi.org/10.1080/10543406.2024.2395532)

## See also

[`bdribs_post`](https://gosukehommaex.github.io/bssd/reference/bdribs_post.md)

## Examples

``` r
time_at_risk(follow_up = c(52, 52, 30, 12),
             event = c(1, 0, 1, 0),
             event_time = c(10, NA, 25, NA))
#>     exposure time_at_risk      n_event 
#>          146           99            2 
```
