# -------------------------------------------------------------------------#
# Hessian source on Bachmann: does the Boehm result extend to 113 parameters?
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The Optimisation vignette measures the Hessian sources on Boehm, nine
# parameters, and states as a limitation that the quasi-Newton advantage is
# expected to narrow as the parameter count grows. This script runs the same
# comparison on Bachmann, 113 parameters, over one shared set of starting
# points, and prints the same two tables.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Fri 05 Sep 2026
#
# [Info]
# The model is built by the shipped example, which builds it by hand and needs
# no libsbml. Fitting runs on helix through distributedComputing(). Run from
# the package root; the job folders go to a scratch directory.
# -------------------------------------------------------------------------#

library(dMod2)

# BACHMANN_SUBMIT=1 submits batch .batch, BACHMANN_COLLECT=1 fetches finished
# results without a model, neither builds only the model. Both are opt in
# because submitting deletes the job folder a running array writes into.
# BACHMANN_LOCAL=1 runs the first block on this machine instead, as a check.
.submit  <- nzchar(Sys.getenv("BACHMANN_SUBMIT"))
.local   <- nzchar(Sys.getenv("BACHMANN_LOCAL"))
.recover <- FALSE    # TRUE re-attaches to a job already in the queue
# One job per Hessian source, each an array of .blocks tasks. Small tasks over
# many array slots rather than few wide jobs: the scheduler backfills them one
# by one, so the fits start without waiting for a block of free cores.
.taskCores <- 4L     # fits per array task, one core each; the demo sets .cores of its own
.blocks  <- 50L      # array tasks per source
.nstart  <- .taskCores * .blocks   # 200 shared starting points per Hessian source
.jobname <- "bachmann_hessianSource"

# Starts are drawn in batches. A batch has its own seed but the same seed for
# every source, so the sources stay comparable; batches accumulate, so two of
# them give 2 * .blocks starts per source. Batch 1 keeps the bare job name.
.batch   <- 4L       # the batch to submit
.batches <- 1:4      # the batches to collect
.jobFor  <- function(source, batch)
  paste0(.jobname, if (batch > 1L) paste0("_b", batch) else "", "_", source)
.collect <- nzchar(Sys.getenv("BACHMANN_COLLECT"))

# Set by the slowest arm, not the typical one: a quasi-Newton method from an
# identity seed needs of the order of n_theta accepted steps, and a truncated
# arm shrinks the shared block set for every other arm too.
.cluster <- if (.local) localNode(cores = .taskCores) else
  slurmCluster("helix", cores = .taskCores, walltime = "03:00:00", memPerCore = 2, partition = "cpu-single")
.jobdir  <- file.path(tempdir(), "bachmann_jobs")
dir.create(.jobdir, showWarnings = FALSE)
inJobdir <- function(expr) {
  old <- setwd(.jobdir)
  on.exit(setwd(old))
  expr
}

# From the shared catalogue; the source is substituted into the task as a literal.
source(file.path("bench", "hessian", "hessianSourceSettings.R"))
settings <- hessianSourceSettings[c("gn", "gn_sr1", "sr1_id", "sr1_gn", "sr1_id_gn")]
# BACHMANN_ONLY names a single source, which is how a source whose remote build
# failed is sent again without touching the ones already in the queue.
.only    <- if (nzchar(Sys.getenv("BACHMANN_ONLY"))) Sys.getenv("BACHMANN_ONLY") else NULL


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Collect only
#
# The manifest of a job lies in its folder on helix, so the results come back
# without the model, from any machine, also while the arrays still run.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
fetchResults <- function(source) {
  unlist(lapply(.batches, function(batch) {
    job <- inJobdir(distributedComputing(NULL, jobname = .jobFor(source, batch), on = .cluster,
                                         recover = TRUE))
    lapply(Filter(Negate(is.null), inJobdir(job$get(combine = FALSE))), function(out) {
      if (is.list(out)) out$batch <- batch
      out
    })
  }), recursive = FALSE)
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Tables, the same two the vignette prints
#
# Blocks of one source are glued back together before ranking, so `index` runs
# over all starts of that source.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The walltime cuts the array off, so sources can come back with different
# blocks. Only blocks that every source returned are kept, which leaves all
# sources on exactly the same starting points.
collect <- function(raw) {
  ok <- Filter(function(x) is.list(x) && !is.null(x$fits) &&
                 inherits(x$fits, "parframe"), raw)
  if (!length(ok)) stop("no block came back")
  bySource <- split(ok, vapply(ok, `[[`, "", "source"))
  # A source that returned nothing at all is dropped rather than compared
  # against, which is what happens while its array is still queued.
  present <- intersect(names(settings), names(bySource))
  if (length(present) < length(settings))
    message("no result yet from: ",
            paste(setdiff(names(settings), present), collapse = ", "))
  bySource <- bySource[present]
  # Block numbers restart in every batch, so a slice is identified by the pair.
  key <- function(g) paste(vapply(g, `[[`, 0, "batch"),
                           vapply(g, `[[`, 0, "block"), sep = ":")
  keys <- Reduce(intersect, lapply(bySource, key))
  message(length(keys), " of ", .blocks * length(.batches),
          " slices complete in every source, ",
          length(keys) * .taskCores, " starts each")
  # mstrust() numbers its fits within a task, so `index` restarts at 1 in every
  # block. The draw order of a start is recovered from the block it sat in,
  # which is what any budget truncation has to follow.
  lapply(bySource, function(g) {
    keep <- g[key(g) %in% keys]
    f <- do.call(rbind, lapply(keep, function(k) {
      fits <- k$fits
      fits$index <- (k$batch - 1L) * .blocks * .taskCores +
        k$block + (seq_len(nrow(fits)) - 1L) * .blocks
      fits
    }))
    f[order(f$value), ]
  })
}

report <- function(runs) {
  best <- min(vapply(runs, function(f) min(f$value), 0.0))

  # evalPerHit is the price of one start that reaches the optimum, and it
  # factorises into the cost of a start and the rate at which starts succeed.
  perStart <- do.call(rbind, lapply(runs, function(f) {
    hits <- sum(f$value <= best + 0.1)
    data.frame(conv       = sum(f$converged),
               atCap      = sum(f$stopReason == "iterlim"),
               starts     = nrow(f),
               hits       = hits,
               neval      = sum(f$neval),
               medEval    = stats::median(f$neval),
               medValue   = stats::median(f$value),
               evalPerHit = if (hits > 0) round(sum(f$neval) / hits) else NA_real_)
  }))
  print(cbind(source = names(runs), perStart, row.names = NULL))

  invisible(perStart)
}

# runs <- collect(raw); report(runs)


if (.collect) {
  report(collect(unlist(lapply(names(settings), fetchResults), recursive = FALSE)))
  quit(save = "no")
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model and objective
#
# The example is sourced up to its PEtab cross-check, which is the part that
# would pull in reticulate and libsbml. distributedComputing() ships the model
# sources the objective names, from the example's output directory.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.exampleFile <- system.file("demo", "BachmannMSB2011.R",
                            package = "dMod2")
.exampleCode <- readLines(.exampleFile)
.upto <- grep("^petab <- importPEtab", .exampleCode)[1] - 1L
# Invisible: the excerpt ends on a plot, which Rscript would print to Rplots.pdf.
invisible(eval(parse(text = .exampleCode[seq_len(.upto)]), envir = globalenv()))

parlower <- setNames(log10(.pars$lowerBound), .pars$parameterId)
parupper <- setNames(log10(.pars$upperBound), .pars$parameterId)
# Prior centre and sampling centre in one. A flat -1 is wrong for a model whose
# published parameters span twelve decades, so anything further than a decade
# away is moved onto its own, which is the order of magnitude a modeller brings.
.published <- setNames(log10(.pars$nominalValue), .pars$parameterId)[outerpars]
stopifnot(!anyNA(.published))
pouter   <- structure(rep(-1, length(outerpars)), names = outerpars)
.moved   <- abs(.published + 1) > 1
pouter[.moved] <- round(.published[.moved])

# The published problem has one informative prior, on the receptor pool.
# A multi-start needs the weak prior of the Boehm setup as well, over everything
# but the error model, or starts drift into the bounds instead of converging.
# The prior must not be able to reorder the optima the comparison is about: with
# the centre above it costs about 1 at the published parameters, a quarter of
# the gap between the local optima seen here. Out of it: the error parameters,
# whose variance it would bias, and init_EpoRJAK2, which already has the
# published informative prior.
.weak <- setdiff(outerpars,
                 c(grep("^sd_", outerpars, value = TRUE), "init_EpoRJAK2"))
obj <- obj + constraintL2(pouter[.weak], sigma = 4, attrName = "prior")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# One shared set of starts, so a comparison differs only in the Hessian source.
# The spread is set on its own rather than from the prior: it decides how hard
# the multi-start problem is, the prior decides where the optimum sits.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
set.seed(20260905 + .batch - 1L)
starts <- msParframe(pouter, n = .nstart, sd = 4)
# At this width a normal draw leaves the box on nearly half of all coordinates,
# and clipping would pile those onto the bounds and turn a random start into a
# corner. Redraw instead, so the spread is kept and every start is interior.
for (.nm in outerpars) {
  .x <- starts[[.nm]]
  .out <- .x <= parlower[[.nm]] | .x >= parupper[[.nm]]
  while (any(.out)) {
    .x[.out] <- rnorm(sum(.out), pouter[[.nm]], 4)
    .out <- .x <= parlower[[.nm]] | .x >= parupper[[.nm]]
  }
  starts[[.nm]] <- .x
}

# Task i of an array takes one block of the shared starts; the source is fixed
# per job.
.slice <- function(block)
  starts[seq(block, nrow(starts), by = .blocks), , drop = FALSE]

fitBlock <- function(source, block)
  list(source = source, block = block,
       fits = as.parframe(do.call(mstrust, c(
         list(objfun = obj, center = .slice(block), fits = nrow(.slice(block)),
              cores = .taskCores, rinit = 0.1, rmax = 10, iterlim = 5000,
              parlower = parlower, parupper = parupper),
         settings[[source]]))))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Run
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# What a task needs; the dot names are not in the default `input`.
.input <- c("fitBlock", ".slice", "starts", ".blocks", ".taskCores", "obj", "settings",
            "parlower", "parupper")

if (.submit || .local) {

  # The source is substituted into the expression, so each job gets its own
  # as a literal rather than a variable the next pass of the loop overwrites.
  # A task is the number of its block.
  submitSource <- function(source) inJobdir(eval(substitute(
    distributedComputing(
      { fitBlock(SRC, task) },
      tasks   = as.list(if (.local) 1L else seq_len(.blocks)),
      on      = .cluster,
      jobname = .jobFor(SRC, .batch),
      input   = .input,
      build   = "compile",
      recover = .recover),
    list(SRC = source)), envir = globalenv()))

  .which <- if (is.null(.only)) names(settings) else .only
  jobs <- lapply(setNames(nm = .which), submitSource)

  if (.recover || .local) {
    if (.local) for (j in jobs) inJobdir(j$wait())
    report(collect(unlist(lapply(jobs, function(j)
      lapply(inJobdir(j$get(combine = FALSE)), function(out) {
        if (is.list(out)) out$batch <- .batch
        out
      })), recursive = FALSE)))
  }

}


