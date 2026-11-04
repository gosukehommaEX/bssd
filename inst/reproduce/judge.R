# Shared helpers for the reproduction scripts in inst/reproduce.
# Source this file from the package root after devtools::load_all() or
# library(bssd).

# Judge reproduced values against published values.
#   tab: data frame with columns id, source, quantity, published, reproduced,
#        tolerance and note
#   explained: data frame with columns id, script and reason
# PASS if |reproduced - published| <= tolerance; EXPLAINED if the item is
# registered in explained.csv; FAIL otherwise; INFO if nothing is published.
judge <- function(tab, explained) {
  tab$abs_diff <- abs(tab$reproduced - tab$published)
  tab$judgement <- ifelse(is.na(tab$published), "INFO",
                          ifelse(tab$abs_diff <= tab$tolerance, "PASS",
                                 "FAIL"))
  hit <- tab$judgement == "FAIL" & tab$id %in% explained$id
  tab$judgement[hit] <- "EXPLAINED"
  tab$note[hit] <- explained$reason[match(tab$id[hit], explained$id)]
  tab[, c("id", "source", "quantity", "published", "reproduced", "abs_diff",
          "tolerance", "judgement", "note")]
}

# Check that the judgement is not vacuous: shifting each published value of
# a PASS item by five times its tolerance must turn the item into FAIL.
self_check <- function(tab) {
  pass <- tab[tab$judgement == "PASS", ]
  shifted <- pass
  shifted$published <- pass$published + 5 * pass$tolerance
  re <- judge(shifted, data.frame(id = character(0), script = character(0),
                                  reason = character(0)))
  n_fail <- sum(re$judgement == "FAIL")
  cat(sprintf("Self-check: %d of %d PASS items become FAIL after shifting ",
              n_fail, nrow(pass)), "the published values.\n", sep = "")
  if (n_fail != nrow(pass)) {
    stop("Self-check failed: the tolerances do not detect a shift.")
  }
  invisible(n_fail)
}

# Read the registry of explained differences.
read_explained <- function(path = "inst/reproduce/explained.csv") {
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# Print the counts of the judgements.
print_counts <- function(tab, script) {
  cnt <- table(factor(tab$judgement,
                      levels = c("PASS", "EXPLAINED", "FAIL", "INFO")))
  cat(script, ":", paste(names(cnt), cnt, sep = " = ", collapse = ", "),
      "\n")
  invisible(cnt)
}
