#!/usr/bin/env Rscript

## =====================================================================
##  run-benchmarks.R: the objective of PEtab benchmark problems in dMod2.
## =====================================================================

##      Rscript benchmarks/run-benchmarks.R --tier tiny
##      Rscript benchmarks/run-benchmarks.R --tier full --asa
##      Rscript benchmarks/run-benchmarks.R --models Boehm,Raia --nrep 9

##  Run `--help` for the options. Results land in
##  benchmarks/results/<tag>/ as results.csv, README.md and two figures.

suppressPackageStartupMessages(library(dMod2))

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()
})

source(file.path(ROOT, "benchmarks", "R", "harness.R"))
bench_source(file.path(ROOT, "benchmarks", "R"))


## ---------------------------------------------------------------------
##  Options
## ---------------------------------------------------------------------

OPTS <- list(
  tier         = "tiny",     # tiny | medium | full
  models       = "",         # comma-separated short names; "" = the tier
  skip         = "",         # comma-separated short names to leave out
  tol          = "1e-6",     # atol = rtol of the timed solves
  nrep         = "5",        # bursts per timing, the minimum is reported
  target       = "0.25",     # seconds per burst
  asa          = "FALSE",    # also CVODES adjoint sensitivities (Sundials backend)
  cores        = "4",        # compile jobs per import
  `petab-root` = "",         # default benchmarks/cache/petab
  outdir       = "",         # default benchmarks/results
  builddir     = "",         # default a fresh temporary directory
  plots        = "TRUE",
  help         = "FALSE"
)

parse_args <- function(args, opts) {
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) { warning("ignoring argument: ", a); i <- i + 1L; next }
    key <- sub("^--", "", a)
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1L]]
      key <- kv[1L]; val <- paste(kv[-1L], collapse = "=")
    } else if (key %in% c("asa", "help", "no-plots")) {
      val <- "TRUE"
    } else {
      i <- i + 1L; val <- if (i <= length(args)) args[i] else ""
    }
    if (key == "no-plots") opts$plots <- "FALSE"
    else if (!key %in% names(opts)) stop("unknown option --", key)
    else opts[[key]] <- val
    i <- i + 1L
  }
  opts
}

OPTS <- parse_args(commandArgs(trailingOnly = TRUE), OPTS)

if (isTRUE(as.logical(OPTS$help))) {
  cat("
dMod2 benchmark suite: value, forward and reverse gradients of PEtab problems

  --tier <tiny|medium|full>   problem set                         [tiny]
        tiny    Boehm, Raia, Lucarelli
        medium  + Bachmann, Isensee
        full    + Chen, Lang
  --models <a,b>              these problems instead of the tier
  --skip <a,b>                leave these out
  --tol <x>                   atol = rtol of the timed solves     [1e-6]
  --nrep <n>                  bursts per timing (minimum reported) [5]
  --target <s>                seconds per burst                    [0.25]
  --asa                       also CVODES adjoints on the Sundials backend
  --cores <n>                 compile jobs per import              [4]
  --petab-root <dir>          the collection  [benchmarks/cache/petab]
  --outdir <dir>              results root    [benchmarks/results]
  --builddir <dir>            where models are compiled  [a temp dir]
  --no-plots                  skip the figures

Every timing runs on one core. The gradient error is scored against the
forward gradient at atol 1e-12, rtol 1e-10.
")
  quit(save = "no")
}

options(dMod.cores = 1)
tol <- as.numeric(OPTS$tol)
root <- if (nzchar(OPTS$`petab-root`)) OPTS$`petab-root` else
  file.path(ROOT, "benchmarks", "cache", "petab")
if (!dir.exists(root))
  stop("no PEtab collection at ", root, "; run benchmarks/fetch-models.R or pass --petab-root")
names_run <- if (nzchar(OPTS$models)) strsplit(OPTS$models, ",")[[1L]] else BENCH_TIERS[[OPTS$tier]]
if (is.null(names_run)) stop("unknown tier ", OPTS$tier)
names_run <- setdiff(names_run, strsplit(OPTS$skip, ",")[[1L]])
problems <- petab_select(petab_index(root), names_run)

outroot <- if (nzchar(OPTS$outdir)) OPTS$outdir else file.path(ROOT, "benchmarks", "results")
outdir <- file.path(outroot, bench_tag(OPTS))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
builddir <- if (nzchar(OPTS$builddir)) OPTS$builddir else tempfile("dmod2-bench-")
dir.create(builddir, recursive = TRUE, showWarnings = FALSE)
writeLines(bench_info(OPTS), file.path(outdir, "run-info.txt"))
cat(bench_info(OPTS), sep = "\n")


## ---------------------------------------------------------------------
##  Run
## ---------------------------------------------------------------------

rows <- list(); traits <- list(); skipped <- character(0)
for (i in seq_len(nrow(problems))) {
  e <- problems[i, ]
  cat(sprintf("\n[%d/%d] %s\n", i, nrow(problems), e$name)); utils::flush.console()
  res <- tryCatch({
    pet <- petab_import(e, builddir, cores = as.integer(OPTS$cores), tol = tol)
    traits[[e$short]] <- bench_traits(pet, e)
    asa <- if (isTRUE(as.logical(OPTS$asa)))
      tryCatch(petab_import(e, builddir, backend = "Sundials",
                            cores = as.integer(OPTS$cores), tol = tol),
               error = function(err) { cat("  no ASA:", conditionMessage(err), "\n"); NULL })
    bench_problem(pet, e, tol, as.integer(OPTS$nrep), as.numeric(OPTS$target), asa)
  }, error = function(err) err)
  if (inherits(res, "error")) {
    skipped[[e$short]] <- conditionMessage(res)
    cat("  skipped:", conditionMessage(res), "\n")
    next
  }
  print(res[, c("arm", "sec", "cost", "err_rel", "err_tau")], digits = 3, row.names = FALSE)
  rows[[e$short]] <- res
  utils::write.csv(do.call(rbind, rows), file.path(outdir, "results.csv"), row.names = FALSE)
}

df <- if (length(rows)) do.call(rbind, rows) else NULL
bench_coverage(traits)
bench_print(df)
figures <- if (isTRUE(as.logical(OPTS$plots))) bench_plots(df, outdir) else character(0)
bench_readme(outdir, OPTS, df, traits, skipped, figures)
cat("\nresults in", outdir, "\n")
