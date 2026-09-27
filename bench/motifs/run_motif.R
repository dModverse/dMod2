# Run one motif: detection with reconstruction, then reduction. Writes
# <outdir>/<motif>.rds with timings and outcome, prints a short report.
#   Rscript run_motif.R <motif> [outdir] [extra symmetryDetection args as R code]
args <- commandArgs(trailingOnly = TRUE)
# DMOD_LOADER: a script that attaches dMod2 (e.g. a dev overlay); else library()
if (nzchar(Sys.getenv("DMOD_LOADER"))) source(Sys.getenv("DMOD_LOADER")) else
  suppressMessages(library(dMod2))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE))))
source(file.path(here, "motifs.R"))
motif  <- args[1]
outdir <- if (length(args) >= 2) args[2] else file.path(here, "results")
extra  <- if (length(args) >= 3) eval(parse(text = paste0("list(", args[3], ")"))) else list()
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

m <- .motifs[[motif]]()
nDir <- m$nDir; m$nDir <- NULL
dargs <- c(m, list(method = "observability", reconstruct = TRUE, cores = 4L,
                   verbose = FALSE), extra)
dargs <- dargs[!duplicated(names(dargs), fromLast = TRUE)]

out <- list(motif = motif, extra = args[3], version = as.character(packageVersion("dMod2")))
t0 <- proc.time()[["elapsed"]]
res <- tryCatch(do.call(symmetryDetection, dargs), error = function(e) e)
out$tDetect <- proc.time()[["elapsed"]] - t0
if (inherits(res, "error")) {
  out$error <- conditionMessage(res)
} else {
  out$rank <- res$rank; out$dim <- res$dim
  out$nDir <- length(res$symmetries); out$nDirExpected <- nDir
  out$dirs <- lapply(res$symmetries, function(d)
    list(type = d$type, support = d$support, explicit = isTRUE(d$explicit),
         reason = d$reason, generator = d$generator, route = d$route))
  out$nOpen <- sum(!vapply(out$dirs, `[[`, logical(1), "explicit"))
  cat(sprintf("[detect done %.1fs: rank %d/%d, %d dirs, %d open]\n", out$tDetect,
              out$rank, out$dim, out$nDir, out$nOpen))
  saveRDS(out, file.path(outdir, paste0(motif, ".rds")))
  if (nzchar(Sys.getenv("SKIP_REDUCE"))) quit(save = "no")
  t1 <- proc.time()[["elapsed"]]
  red <- tryCatch(symmetryReduction(res), error = function(e) e)
  out$tReduce <- proc.time()[["elapsed"]] - t1
  if (inherits(red, "error")) out$reduceError <- conditionMessage(red)
  else {
    out$removed <- red$removed; out$remaining <- red$remaining
    out$trafo <- red$trafo
    # the reduced model must be identifiable: compose the reduction into the motif's
    # own trafo (initial values) and detect again
    if (!length(red$remaining) && !is.null(red$trafo)) {
      rt <- red$trafo
      tr <- if (is.null(m$trafo)) rt else {
        mt <- m$trafo
        inner <- setNames(cOde::replaceSymbols(names(rt), rt, as.character(mt)), names(mt))
        do.call(eqnvec, as.list(c(inner, rt[setdiff(names(rt), names(mt))])))
      }
      d2 <- dargs; d2$trafo <- tr; d2$reconstruct <- FALSE
      t2 <- proc.time()[["elapsed"]]
      chk <- tryCatch(do.call(symmetryDetection, d2), error = function(e) e)
      out$tCheck <- proc.time()[["elapsed"]] - t2
      out$reducedIdentifiable <- if (inherits(chk, "error")) NA else isTRUE(chk$identifiable)
    }
  }
}
saveRDS(out, file.path(outdir, paste0(motif, ".rds")))
cat(sprintf("%-16s detect %7.1fs  rank %s/%s  dirs %s (open %s)  reduce %s  remaining: %s  reduced identifiable: %s %s\n",
            motif, out$tDetect, out$rank %||% NA, out$dim %||% NA, out$nDir %||% NA,
            out$nOpen %||% NA,
            if (is.null(out$tReduce)) "-" else sprintf("%.1fs", out$tReduce),
            paste(out$remaining, collapse = ","), format(out$reducedIdentifiable %||% NA),
            if (!is.null(out$error)) paste("ERROR:", out$error)
            else if (!is.null(out$reduceError)) paste("RED-ERROR:", out$reduceError) else ""))
