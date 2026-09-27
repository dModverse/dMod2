# Which route closed each general direction: the default run next to one with the jets
# tried first (DMOD_SYM_JETFIRST). Rscript routes.R <default dir> <jets-first dir>
args <- commandArgs(trailingOnly = TRUE)
rd <- function(dir) {
  fs <- list.files(dir, pattern = "rds$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    o <- readRDS(f)
    g <- Filter(function(d) !identical(d$type, "scaling"), o$dirs)
    if (!length(g)) return(NULL)
    data.frame(motif = o$motif, support = vapply(g, function(d) length(d$support), 1L),
               key = vapply(g, function(d) paste(sort(d$support), collapse = ","), ""),
               route = vapply(g, function(d) if (is.null(d$route)) "?" else d$route, ""),
               stringsAsFactors = FALSE)
  }))
}
a <- rd(args[1]); b <- rd(args[2])
m <- merge(a, b, by = c("motif", "key", "support"), suffixes = c(".default", ".jetsFirst"),
           all = TRUE)
m$key <- NULL
print(m, row.names = FALSE, right = FALSE)
kind <- function(r) ifelse(grepl("^jets", r), "jets", ifelse(grepl("fit|read-off", r), "fit", r))
cat("\ndefault:     ", paste(names(table(kind(m$route.default))), table(kind(m$route.default)),
                              collapse = ", "), "\n")
cat("jets first:  ", paste(names(table(kind(m$route.jetsFirst))), table(kind(m$route.jetsFirst)),
                              collapse = ", "), "\n")
