#!/usr/bin/env Rscript

## =====================================================================
##  compare-results.R: two runs of run-benchmarks.R side by side.
## =====================================================================

##      Rscript benchmarks/compare-results.R <base>/results.csv <new>/results.csv

##  Per arm the geometric mean of new / base over the problems both ran, then
##  one row per problem and arm. Ratios survive a change of machine load better
##  than absolute times; still compare runs from one machine.

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2L) stop("usage: compare-results.R <base results.csv> <new results.csv>")
base <- utils::read.csv(a[1]); new <- utils::read.csv(a[2])
m <- merge(base, new, by = c("model", "arm"), suffixes = c(".base", ".new"))
if (!nrow(m)) stop("the two runs share no problem")
m$time <- m$sec.new / m$sec.base
m$err  <- m$err_rel.new / m$err_rel.base
gm <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))

cat("\ngeometric mean new / base over", length(unique(m$model)), "problems\n")
for (arm in unique(m$arm)) {
  z <- m[m$arm == arm, ]
  cat(sprintf("  %-8s time %6.3f   error %s\n", arm, gm(z$time),
              if (arm == "value") "-" else sprintf("%6.3f", gm(z$err))))
}
cat("\nper problem\n")
out <- m[order(m$n_theta.base, m$arm), c("model", "arm", "sec.base", "sec.new", "time",
                                          "err_rel.base", "err_rel.new")]
out[-(1:2)] <- lapply(out[-(1:2)], signif, 3)
print(out, row.names = FALSE)
gone <- setdiff(unique(base$model), unique(new$model))
if (length(gone)) cat("\nonly in base:", paste(gone, collapse = ", "), "\n")
