# Reproduction of Table 1 and Figure 1 of
#   Duncan B (2026). Robust safety monitoring and signal detection using
#   alternatives to the standard Poisson distribution. Journal of
#   Biopharmaceutical Statistics 36(2):227-244.
#   doi:10.1080/10543406.2024.2395532
#
# Run from the package root after devtools::load_all() or library(bssd):
#   source("inst/reproduce/reproduce_Duncan2026.R")
# Output (inst/reproduce/results/):
#   reproduce_Duncan2026.csv   judgement of each item
#   Duncan2026_figure1.csv     values of Figure 1
#   Duncan2026_figure1.pdf     Figure 1 redrawn
# Run time: a few seconds.
#
# Table 1 gives the log-logistic parameters (mu, sigma) of each scenario and
# the cumulative incidence (IP) at week 52 or 104. The incidences are
# recomputed with ll_cdf() and compared with the published values, which are
# rounded to 0.1% (tolerance 0.0005). The descriptions of the scenarios in the
# text are rounded to 1% (tolerance 0.005). Time is in weeks.
#
# Figures 2 to 10 are simulation results; they are reproduced approximately
# by inst/reproduce/reproduce_Duncan2026_simulation.R.

source("inst/reproduce/judge.R")
explained <- read_explained()

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}
tol_t <- 0.0005 + 1e-9
tol_x <- 0.005 + 1e-9

# Table 1: decreasing hazard (scenarios 4 and 5)
add("D_T1_04P", "Table 1", "scenario 4 placebo IP at week 104 (9.85, 2.0)",
    0.050, ll_cdf(104, 9.85, 2.0), tol_t)
add("D_T1_04P52", "Table 1",
    "scenario 4 placebo IP at week 52 (9.85, 2.0)", 0.050,
    ll_cdf(52, 9.85, 2.0), tol_t,
    "published 5.0% read as the week 52 value")
add("D_T1_04D", "Table 1",
    "scenario 4 largest study drug IP at week 104 (7.4, 2.0)", 0.151,
    ll_cdf(104, 7.4, 2.0), tol_t)
add("D_T1_04D52", "Table 1",
    "scenario 4 largest study drug IP at week 52 (7.4, 2.0)", 0.151,
    ll_cdf(52, 7.4, 2.0), tol_t,
    "published 15.1% read as the week 52 value")
add("D_T1_04Dm", "Table 1",
    "scenario 4 middle study drug IP at weeks 52 and 104 (8.35, 2.0)", NA,
    ll_cdf(52, 8.35, 2.0), NA,
    sprintf("week 52: %.4f; week 104: %.4f", ll_cdf(52, 8.35, 2.0),
            ll_cdf(104, 8.35, 2.0)))
add("D_T1_051", "Table 1", "scenario 5.1 IP at week 52 (6.17, 1.6)", 0.200,
    ll_cdf(52, 6.17, 1.6), tol_t)
add("D_T1_052", "Table 1",
    "scenario 5.2 constant hazard rate with the same week 52 IP", NA,
    -log(1 - ll_cdf(52, 6.17, 1.6)) / 52, NA, "hazard rate per week")

# Table 1: increasing hazard (scenarios 6 to 8)
add("D_T1_06P", "Table 1", "scenario 6 placebo IP at week 104 (7.3, 0.8)",
    0.035, ll_cdf(104, 7.3, 0.8), tol_t)
add("D_T1_06D", "Table 1",
    "scenario 6 largest study drug IP at week 104 (5.85, 0.8)", 0.181,
    ll_cdf(104, 5.85, 0.8), tol_t)
add("D_T1_06Dm", "Table 1",
    "scenario 6 other study drug IPs at week 104 (6.5, 6.1; 0.8)", NA,
    ll_cdf(104, 6.5, 0.8), NA,
    sprintf("mu = 6.5: %.4f; mu = 6.1: %.4f", ll_cdf(104, 6.5, 0.8),
            ll_cdf(104, 6.1, 0.8)))
add("D_T1_07P", "Table 1", "scenario 7 placebo IP at week 104 (5.85, 0.7)",
    0.152, ll_cdf(104, 5.85, 0.7), tol_t)
add("D_T1_07D", "Table 1",
    "scenario 7 largest study drug IP at week 104 (5.07, 0.7)", 0.353,
    ll_cdf(104, 5.07, 0.7), tol_t)
add("D_T1_07Dm", "Table 1",
    "scenario 7 other study drug IPs at week 104 (5.5, 0.7; 5.15, 0.6)", NA,
    ll_cdf(104, 5.5, 0.7), NA,
    sprintf("(5.5, 0.7): %.4f; (5.15, 0.6): %.4f", ll_cdf(104, 5.5, 0.7),
            ll_cdf(104, 5.15, 0.6)))
add("D_T1_081P", "Table 1",
    "scenario 8.1 true placebo IP at week 104 (6.5, 0.8)", 0.090,
    ll_cdf(104, 6.5, 0.8), tol_t)
add("D_T1_081A", "Table 1 footnote",
    "scenario 8.1 assumed placebo IP at week 104 (6.85, 0.8)", 0.060,
    ll_cdf(104, 6.85, 0.8), tol_t)
add("D_T1_082P", "Table 1",
    "scenario 8.2 true placebo IP at week 104 (6.845, 0.8)", 0.060,
    ll_cdf(104, 6.845, 0.8), tol_t)
add("D_T1_082A", "Table 1 footnote",
    "scenario 8.2 assumed placebo IP at week 104 (6.6, 0.8)", 0.080,
    ll_cdf(104, 6.6, 0.8), tol_t)
add("D_T1_082D", "Table 1",
    "scenario 8.2 study drug IP at week 104 (5.85, 0.8)", 0.181,
    ll_cdf(104, 5.85, 0.8), tol_t)

# descriptions of the scenarios in the text
add("D_TX_04P", "Section 4.2.1", "scenario 4 placebo IP at week 104", 0.04,
    ll_cdf(104, 9.85, 2.0), tol_x)
add("D_TX_04D", "Section 4.2.1",
    "scenario 4 largest study drug IP at week 104", 0.12,
    ll_cdf(104, 7.4, 2.0), tol_x)
add("D_TX_05", "Section 4.2.1", "scenario 5 IP at week 52", 0.20,
    ll_cdf(52, 6.17, 1.6), tol_x)
add("D_TX_06P", "Section 4.2.2", "scenario 6 placebo IP at week 104", 0.05,
    ll_cdf(104, 7.3, 0.8), tol_x)
add("D_TX_06D", "Section 4.2.2",
    "scenario 6 largest study drug IP at week 104", 0.15,
    ll_cdf(104, 5.85, 0.8), tol_x)
add("D_TX_07P", "Section 4.2.2", "scenario 7 placebo IP at week 104", 0.15,
    ll_cdf(104, 5.85, 0.7), tol_x)
add("D_TX_07D", "Section 4.2.2",
    "scenario 7 largest study drug IP at week 104", 0.35,
    ll_cdf(104, 5.07, 0.7), tol_x)

# Figure 1: expected mean and variance of event counts and incidence counts
# for N = 100 subjects followed for 52 weeks, with rates from 0.001 to 0.018
# per patient-week under a constant hazard (Section 2.2); the rate 0 is added
# because the curves in Figure 1 start at the origin
lambda <- c(0, seq(0.001, 0.018, by = 0.001))
ip <- 1 - exp(-lambda * 52)
fig1 <- data.frame(lambda = lambda, incidence_proportion = ip,
                   event_mean = 100 * 52 * lambda,
                   event_variance = 100 * 52 * lambda,
                   incidence_mean = 100 * ip,
                   incidence_variance = 100 * ip * (1 - ip))
add("D_F1_1", "Figure 1", "largest incidence proportion on the x-axis", NA,
    max(ip), NA, "the x-axis of Figure 1 runs to 70%")
add("D_F1_2", "Figure 1", "largest expected event count", NA,
    max(fig1$event_mean), NA, "the y-axis of Figure 1 runs to 100")
add("D_F1_3", "Section 2.2",
    "variance-to-mean ratio of the incidence count at 15% incidence", NA,
    1 - 0.15, NA, "equals 1 - incidence proportion")

tab <- judge(do.call(rbind, rows), explained)
dir.create("inst/reproduce/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/reproduce/results/reproduce_Duncan2026.csv",
                 row.names = FALSE)
utils::write.csv(fig1, "inst/reproduce/results/Duncan2026_figure1.csv",
                 row.names = FALSE)

grDevices::pdf("inst/reproduce/results/Duncan2026_figure1.pdf", width = 7,
               height = 5)
graphics::plot(100 * fig1$incidence_proportion, fig1$event_mean, type = "l",
               lwd = 2, col = "black", ylim = c(0, 100), xlim = c(0, 70),
               xlab = "Incidence proportion (%)", ylab = "Counts",
               main = paste("Expected mean and variance of event counts and",
                            "incidence counts\n(N = 100 subjects, t = 52",
                            "weeks)"))
graphics::lines(100 * fig1$incidence_proportion, fig1$incidence_mean,
                lwd = 2, col = "steelblue")
graphics::lines(100 * fig1$incidence_proportion, fig1$incidence_variance,
                lwd = 2, lty = 2, col = "steelblue")
graphics::legend("topleft", bty = "n", lwd = 2,
                 col = c("black", "steelblue", "steelblue"),
                 lty = c(1, 1, 2),
                 legend = c("Event count (mean = variance)",
                            "Incidence count (mean)",
                            "Incidence count (variance)"))
grDevices::dev.off()

print_counts(tab, "reproduce_Duncan2026")
self_check(tab)
