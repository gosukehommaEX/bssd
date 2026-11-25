# Check that a transcription of the authors' code of Gould and Wang (2017)
# reproduces the published critical counts of Tables III and V. The
# transcription (Rpackage/tools/external/gould2017_author_grid.py, outside
# the package because the code has no stated license) wrote its counts to
# python/ref_gould2017_author_grid.csv. It evaluates the posterior on
# the fixed grids of the authors' code; the differences between these counts
# and the exact counts of bssd (inst/reproduce/reproduce_Gould2017.R) are
# therefore explained by the grids.
#
# Run from the package root:
#   source("inst/validation/validate_gould2017_author_grid.R")
# Output: inst/validation/results/validate_gould2017_author_grid.csv
# Run time: seconds (reads python/ref_gould2017_author_grid.csv).

source("inst/reproduce/judge.R")
explained <- read_explained("inst/validation/explained.csv")
g <- utils::read.csv("inst/validation/python/ref_gould2017_author_grid.csv",
                     stringsAsFactors = FALSE)
# counts are integers: the tolerance 0.5 requires equal counts
tab <- data.frame(id = paste0("AG_", g$id), source = "Gould and Wang (2017)",
                  quantity = "critical count of the authors' code (Python transcription)",
                  published = g$published, reproduced = g$author_grid,
                  tolerance = 0.5, note = "", stringsAsFactors = FALSE)
tab <- judge(tab, explained)
print_counts(tab, "validate_gould2017_author_grid")
out <- tab
names(out)[names(out) == "published"] <- "reference"
names(out)[names(out) == "reproduced"] <- "bssd"
dir.create("inst/validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out,
                 "inst/validation/results/validate_gould2017_author_grid.csv",
                 row.names = FALSE)
# the results are saved first so that a failed self-check can be inspected
self_check(tab)
