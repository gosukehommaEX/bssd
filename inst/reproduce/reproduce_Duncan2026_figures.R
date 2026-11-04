# Figures and judgements for the approximate reproduction of Figures 2 to 10
# of Duncan (2026) from the results of reproduce_Duncan2026_simulation.R.
#
# Run from the package root after the simulation script:
#   source("inst/reproduce/reproduce_Duncan2026_figures.R")
# Output (inst/reproduce/results/):
#   reproduce_Duncan2026_simulation.csv  judgement of each item
#   Duncan2026_figures.pdf               reproduced figures (one per page)
#
# The settings of the simulations in the article are not fully reported (see
# the header of reproduce_Duncan2026_simulation.R), and the published values
# were read from the figures (inst/reproduce/Duncan2026_digitized.csv, about
# 0.015 reading error). The reproduced values are therefore listed next to
# the published values as INFO items, and only the statements of the article
# about the comparison of the models are judged (PASS or FAIL). The criteria
# of these statements are written in the items below.

source("inst/reproduce/judge.R")
explained <- read_explained()
sim <- readRDS("inst/reproduce/results/Duncan2026_simulation.rds")
dig <- utils::read.csv("inst/reproduce/Duncan2026_digitized.csv",
                       stringsAsFactors = FALSE)
out_by_id <- stats::setNames(sim$out, vapply(sim$out, function(o) o$id, ""))

task_id <- function(figure, point) {
  switch(as.character(figure),
         "2" = sprintf("F2_S1_P%d", point),
         "3" = sprintf("F3_S2_P%d", point),
         "4" = sprintf("F4_S3_P%d", point),
         "5" = sprintf("F5_S4_P%d", point),
         "6" = sprintf("F5_S4_P%d", point),
         "7" = "F7_S5.1_P1",
         "8" = "F8_S5.2_P1",
         "9" = sprintf("F9_S6_P%d", point),
         "10" = sprintf("F10_S7_P%d", point),
         "text_4" = "F5_S4_P1",
         "text_8.1" = "Ftext_8.1_S8.1_P1",
         "text_8.2" = "Ftext_8.2_S8.2_P1")
}
res_of <- function(figure, point, model, look) {
  r <- out_by_id[[task_id(figure, point)]]$results
  r[r$model == model & r$look == look, ]
}
power <- function(figure, point, model, cutoff, look = 1) {
  x <- res_of(figure, point, model, look)$prob
  x <- x[!is.na(x)]
  p <- mean(x > cutoff)
  c(value = p, se = sqrt(p * (1 - p) / length(x)))
}
est_rd <- function(figure, point, model, look = 1) {
  x <- res_of(figure, point, model, look)$ip_difference
  c(value = mean(x, na.rm = TRUE), se = stats::sd(x, na.rm = TRUE) /
      sqrt(sum(!is.na(x))))
}
obs_rd <- function(figure, point, look = 1) {
  r <- out_by_id[[task_id(figure, point)]]$results
  x <- r$observed_difference[r$look == look & r$model == r$model[1]]
  c(value = mean(x), se = stats::sd(x) / sqrt(length(x)))
}

rows <- list()
add <- function(id, source, quantity, published, reproduced, tolerance,
                note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    id = id, source = source, quantity = quantity, published = published,
    reproduced = reproduced, tolerance = tolerance, note = note,
    stringsAsFactors = FALSE)
}

# published values read from the figures and the text, next to the
# reproduced values (INFO)
for (i in seq_len(nrow(dig))) {
  d <- dig[i, ]
  v <- switch(d$quantity,
              power = power(d$figure, d$point, d$model, d$cutoff, d$look),
              rd = est_rd(d$figure, d$point, d$model, d$look),
              expected_rd = obs_rd(d$figure, d$point, d$look))
  src <- if (grepl("^text", d$figure)) "Section 4.2" else
    paste("Figure", d$figure)
  add(sprintf("D_SIM_%s_P%d_L%d_%s_%s", d$figure, d$point, d$look, d$model,
              d$quantity),
      src, sprintf("%s of %s (point %d, look %d%s)", d$quantity, d$model,
                   d$point, d$look,
                   if (is.na(d$cutoff)) "" else
                     sprintf(", threshold %g", d$cutoff)),
      d$value, unname(v["value"]), NA,
      sprintf("Monte Carlo SE %.4f%s", v["se"],
              if (is.na(d$note) || d$note == "") "" else
                paste0("; ", d$note)))
}

# statements of the article about the comparison of the models (PASS/FAIL)
statement <- function(id, source, text, ok) {
  add(id, source, text, 1, as.numeric(isTRUE(ok)), 0.5)
}
pw <- function(f, p, m, c, l = 1) power(f, p, m, c, l)[["value"]]
rd <- function(f, p, m, l = 1) est_rd(f, p, m, l)[["value"]]
ob <- function(f, p, l = 1) obs_rd(f, p, l)[["value"]]

statement("D_Q_F2", "Section 4.1, Figure 2",
          paste("low incidence: the three models give nearly the same power",
                "(all pairwise differences <= 0.03)"),
          all(vapply(1:4, function(p) {
            v <- c(pw(2, p, "EXP_BERN", 0.925), pw(2, p, "BDRIBS1", 0.925),
                   pw(2, p, "BDRIBS2", 0.925))
            max(v) - min(v) <= 0.03
          }, logical(1))))
statement("D_Q_F3a", "Section 4.1, Figure 3",
          paste("medium incidence: BDRIBS(1) has lower power than EXP_BERN",
                "and BDRIBS(2) at the two largest effects (by >= 0.05)"),
          all(vapply(3:4, function(p) {
            pw(3, p, "BDRIBS1", 0.925) <=
              min(pw(3, p, "EXP_BERN", 0.925), pw(3, p, "BDRIBS2", 0.925)) -
              0.05
          }, logical(1))))
statement("D_Q_F3b", "Section 4.1, Figure 3",
          paste("medium incidence: EXP_BERN and BDRIBS(2) give similar",
                "power (difference <= 0.05)"),
          all(vapply(1:4, function(p) {
            abs(pw(3, p, "EXP_BERN", 0.925) - pw(3, p, "BDRIBS2", 0.925)) <=
              0.05
          }, logical(1))))
statement("D_Q_F4a", "Section 4.1, Figure 4",
          "high incidence: BDRIBS(1) has almost no power (<= 0.02)",
          all(vapply(1:4, function(p) pw(4, p, "BDRIBS1", 0.95) <= 0.02,
                     logical(1))))
statement("D_Q_F4b", "Section 4.1, Figure 4",
          paste("high incidence: EXP_BERN has at least the power of",
                "BDRIBS(2) (slack 0.01)"),
          all(vapply(2:4, function(p) {
            pw(4, p, "EXP_BERN", 0.95) >= pw(4, p, "BDRIBS2", 0.95) - 0.01
          }, logical(1))))
statement("D_Q_F5a", "Section 4.2.1, Figure 5",
          paste("decreasing hazard: the Type I errors of EXP_BERN and",
                "BDRIBS(2) are inflated (> 0.10)"),
          pw(5, 1, "EXP_BERN", 0.95) > 0.10 && pw(5, 1, "BDRIBS2", 0.95) > 0.10)
statement("D_Q_F5b", "Section 4.2.1, Figure 5",
          paste("decreasing hazard: the Type I error of LL_BERN is",
                "acceptable (<= 0.06)"),
          pw(5, 1, "LL_BERN", 0.95) <= 0.06)
statement("D_Q_F6a", "Section 4.2.1, Figure 6",
          paste("decreasing hazard: EXP_BERN and BDRIBS(2) overestimate the",
                "risk difference (by >= 0.01)"),
          all(vapply(1:3, function(p) {
            min(rd(6, p, "EXP_BERN"), rd(6, p, "BDRIBS2")) - ob(6, p) >= 0.01
          }, logical(1))))
statement("D_Q_F6b", "Section 4.2.1, Figure 6",
          paste("decreasing hazard: LL_BERN estimates the risk difference",
                "within 0.015"),
          all(vapply(1:3, function(p) abs(rd(6, p, "LL_BERN") - ob(6, p)) <=
                       0.015, logical(1))))
statement("D_Q_F7a", "Section 4.2.1, Figure 7",
          paste("decreasing hazard, no effect: EXP_BERN and BDRIBS(2)",
                "overestimate the risk difference at look 1 (>= 0.03) and",
                "less at look 3"),
          all(vapply(c("EXP_BERN", "BDRIBS2"), function(m) {
            rd(7, 1, m, 1) - ob(7, 1, 1) >= 0.03 &&
              rd(7, 1, m, 3) - ob(7, 1, 3) < rd(7, 1, m, 1) - ob(7, 1, 1)
          }, logical(1))))
statement("D_Q_F7b", "Section 4.2.1, Figure 7",
          paste("decreasing hazard, no effect: LL_BERN estimates the risk",
                "difference within 0.01 at every look"),
          all(vapply(1:3, function(l) {
            abs(rd(7, 1, "LL_BERN", l) - ob(7, 1, l)) <= 0.01
          }, logical(1))))
statement("D_Q_F7c", "Section 4.2.1",
          paste("decreasing hazard, no effect: false positive rate of",
                "EXP_BERN or BDRIBS(2) > 0.40 at some look (threshold",
                "0.925)"),
          max(vapply(1:3, function(l) {
            max(pw(7, 1, "EXP_BERN", 0.925, l), pw(7, 1, "BDRIBS2", 0.925, l))
          }, numeric(1))) > 0.40)
statement("D_Q_F8a", "Section 4.2.1, Figure 8",
          paste("constant hazard, no effect: LL_BERN underestimates the",
                "risk difference at look 1 (<= -0.02) and less at look 3"),
          rd(8, 1, "LL_BERN", 1) - ob(8, 1, 1) <= -0.02 &&
            abs(rd(8, 1, "LL_BERN", 3) - ob(8, 1, 3)) <
            abs(rd(8, 1, "LL_BERN", 1) - ob(8, 1, 1)))
statement("D_Q_F8b", "Section 4.2.1, Figure 8",
          paste("constant hazard, no effect: EXP_BERN and BDRIBS(2)",
                "estimate the risk difference within 0.015"),
          all(vapply(1:3, function(l) {
            max(abs(rd(8, 1, "EXP_BERN", l) - ob(8, 1, l)),
                abs(rd(8, 1, "BDRIBS2", l) - ob(8, 1, l))) <= 0.015
          }, logical(1))))
statement("D_Q_F9", "Section 4.2.2, Figure 9",
          paste("increasing hazard: power of LL_BERN > EXP_BERN >=",
                "BDRIBS(2) (slack 0.01) at the two largest effects"),
          all(vapply(3:4, function(p) {
            pw(9, p, "LL_BERN", 0.95) > pw(9, p, "EXP_BERN", 0.95) &&
              pw(9, p, "EXP_BERN", 0.95) >= pw(9, p, "BDRIBS2", 0.95) - 0.01
          }, logical(1))))
statement("D_Q_F10", "Section 4.2.2, Figure 10",
          paste("increasing hazard, high incidence: power of LL_BERN >",
                "EXP_BERN >= BDRIBS(2) (slack 0.01) at the two largest",
                "effects"),
          all(vapply(3:4, function(p) {
            pw(10, p, "LL_BERN", 0.925) > pw(10, p, "EXP_BERN", 0.925) &&
              pw(10, p, "EXP_BERN", 0.925) >= pw(10, p, "BDRIBS2", 0.925) - 0.01
          }, logical(1))))
statement("D_Q_S81", "Section 4.2.2",
          paste("placebo rate underestimated: false positive rates of the",
                "three models are small (<= 0.06, threshold 0.95)"),
          max(pw("text_8.1", 1, "BDRIBS2", 0.95),
              pw("text_8.1", 1, "EXP_BERN", 0.95),
              pw("text_8.1", 1, "LL_BERN", 0.95)) <= 0.06)
statement("D_Q_S82", "Section 4.2.2",
          paste("placebo rate overestimated: LL_BERN detects the signal",
                "more often than BDRIBS(2) (threshold 0.925)"),
          pw("text_8.2", 1, "LL_BERN", 0.925) >
            pw("text_8.2", 1, "BDRIBS2", 0.925))

tab <- judge(do.call(rbind, rows), explained)
utils::write.csv(tab,
                 "inst/reproduce/results/reproduce_Duncan2026_simulation.csv",
                 row.names = FALSE)

# figures: reproduced values (filled symbols and lines) and values read from
# the article (open symbols)
cols <- c(EXP_BERN = "firebrick", BDRIBS1 = "darkorange", BDRIBS2 = "steelblue",
          LL_BERN = "forestgreen", SIM = "black")
draw_power <- function(figure, points, models, cutoff, xlab, xs, main) {
  ys <- sapply(models, function(m) vapply(points, function(p) {
    pw(figure, p, m, cutoff)
  }, numeric(1)))
  graphics::matplot(xs, ys, type = "b", pch = 19, lty = 1,
                    col = cols[models], ylim = c(0, 1), xlab = xlab,
                    ylab = "Probability of exceeding the threshold",
                    main = main)
  dd <- dig[dig$figure == as.character(figure) & dig$quantity == "power", ]
  for (m in models) {
    v <- dd$value[dd$model == m][points]
    graphics::points(xs, v, pch = 1, cex = 1.3, col = cols[m])
  }
  graphics::legend("topleft", bty = "n", col = cols[models], pch = 19,
                   lty = 1, legend = models, cex = 0.8)
  graphics::legend("left", bty = "n", pch = c(19, 1),
                   legend = c("reproduced", "read from the article"),
                   cex = 0.8)
}
draw_rd <- function(figure, models, xs, looks, points, xlab, main) {
  ys <- sapply(models, function(m) {
    vapply(seq_along(xs), function(i) rd(figure, points[i], m, looks[i]),
           numeric(1))
  })
  obs <- vapply(seq_along(xs), function(i) ob(figure, points[i], looks[i]),
                numeric(1))
  rng <- range(c(ys, obs, dig$value[dig$figure == as.character(figure)]))
  graphics::matplot(xs, cbind(obs, ys), type = "b", pch = 19, lty = 1,
                    col = cols[c("SIM", models)], ylim = rng + c(-0.01, 0.01),
                    xlab = xlab, ylab = "Estimated risk difference",
                    main = main)
  dd <- dig[dig$figure == as.character(figure), ]
  for (m in c("SIM", models)) {
    v <- dd$value[dd$model == m]
    graphics::points(xs, v, pch = 1, cex = 1.3, col = cols[m])
  }
  graphics::legend("topright", bty = "n", col = cols[c("SIM", models)],
                   pch = 19, lty = 1, cex = 0.8,
                   legend = c("risk difference in the simulated data",
                              models))
}
xs_rd <- function(figure, points) {
  vapply(points, function(p) ob(figure, p), numeric(1))
}
grDevices::pdf("inst/reproduce/results/Duncan2026_figures.pdf", width = 7,
               height = 5.5)
draw_power(2, 1:4, c("EXP_BERN", "BDRIBS1", "BDRIBS2"), 0.925,
           "Total sample size", c(150, 225, 300, 375),
           "Figure 2: low incidence (threshold 0.925, interim 2)")
draw_power(3, 1:4, c("EXP_BERN", "BDRIBS1", "BDRIBS2"), 0.925,
           "Risk difference in the simulated data", xs_rd(3, 1:4),
           "Figure 3: medium incidence (threshold 0.925, interim 4)")
draw_power(4, 1:4, c("EXP_BERN", "BDRIBS1", "BDRIBS2"), 0.95,
           "Risk difference in the simulated data", xs_rd(4, 1:4),
           "Figure 4: high incidence (threshold 0.95, final)")
draw_power(5, 1:3, c("EXP_BERN", "BDRIBS2", "LL_BERN"), 0.95,
           "Risk difference in the simulated data", xs_rd(5, 1:3),
           "Figure 5: decreasing hazard (threshold 0.95, interim 3)")
draw_rd(6, c("EXP_BERN", "BDRIBS2", "LL_BERN"), xs_rd(5, 1:3), c(1, 1, 1),
        1:3, "Risk difference in the simulated data",
        "Figure 6: estimation, decreasing hazard (interim 3)")
draw_rd(7, c("EXP_BERN", "BDRIBS2", "LL_BERN"), 1:3, 1:3, c(1, 1, 1),
        "Interim look", "Figure 7: estimation, decreasing hazard, no effect")
draw_rd(8, c("EXP_BERN", "BDRIBS2", "LL_BERN"), 1:3, 1:3, c(1, 1, 1),
        "Interim look", "Figure 8: estimation, constant hazard, no effect")
draw_power(9, 1:4, c("EXP_BERN", "BDRIBS2", "LL_BERN"), 0.95,
           "Risk difference in the simulated data", xs_rd(9, 1:4),
           "Figure 9: increasing hazard (threshold 0.95, interim 3)")
draw_power(10, 1:4, c("EXP_BERN", "BDRIBS2", "LL_BERN"), 0.925,
           "Risk difference in the simulated data", xs_rd(10, 1:4),
           "Figure 10: increasing hazard (threshold 0.925, final)")
grDevices::dev.off()

print_counts(tab, "reproduce_Duncan2026_simulation")
self_check(tab)
