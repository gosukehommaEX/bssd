# Data generation of the simulated case study and the simulation study of
#   Schnell PM, Ball G (2016). A Bayesian exposure-time method for clinical
#   trial safety monitoring with blinded data. Therapeutic Innovation &
#   Regulatory Science 50(6):833-838. doi:10.1177/2168479016656702
# written from the functions generate.data() and data.at() of the authors'
# code (supplementary file 0-functions.R), with the same calls to the random
# number generator in the same order, so that a seed gives the same data.
# Sourced by reproduce_Schnell2016.R and reproduce_Schnell2016_simulation.R.

# All subjects, their arms (1 = treatment) and the trial times of their
# events.
schnell_generate <- function(seed, n, r, accrual, follow_up, rate_c,
                             rate_d) {
  set.seed(seed)
  enrolled <- runif(n, 0, accrual)
  arm <- rbinom(n, 1, r)
  events <- rpois(n, ((1 - arm) * rate_c + arm * rate_d) * follow_up)
  ev_id <- integer(0)
  ev_time <- numeric(0)
  for (i in seq_len(n)) {
    if (events[i] > 0) {
      ts <- runif(events[i], 0, follow_up)
      ev_id <- c(ev_id, rep(i, events[i]))
      ev_time <- c(ev_time, enrolled[i] + ts)
    }
  }
  list(enrolled = enrolled, arm = arm, ev_id = ev_id, ev_time = ev_time,
       follow_up = follow_up)
}

# Blinded data at a trial time: follow-up time and number of events of each
# enrolled subject (as data.at() of the authors; the follow-up is not capped
# at follow_up, which does not matter before the end of accrual plus
# follow-up).
schnell_data_at <- function(dat, time) {
  idx <- which(dat$enrolled <= time)
  past <- dat$ev_time <= time
  y <- tabulate(dat$ev_id[past], nbins = length(dat$enrolled))[idx]
  data.frame(t = time - dat$enrolled[idx], y = y)
}
