# numSteady bench: Benchmark-Models-PEtab problems.
#
#   Rscript bench/numSteady/petab.R <benchmark dir> <out.rds> [problems] [nevals]
#
# Per problem: import, objective at the published parameters, predictions
# against the collection's simulatedData table (observable scale, relative to
# max(|sim|, 1e-8)), and nevals objective evaluations with gradient at
# parameters perturbed by N(0, 0.1) on the estimation scale, counting failures.
# The dMod2 on the library path decides which steady-state method runs.

args  <- commandArgs(TRUE)
bdir  <- args[1]
out   <- args[2]
probs <- if (length(args) >= 3L && nzchar(args[3])) strsplit(args[3], ",")[[1]] else list.dirs(bdir, FALSE, FALSE)
nev   <- if (length(args) >= 4L) as.integer(args[4]) else 20L
suppressMessages(library(dMod2))

simCheck <- function(pp, sim, conds) {
  meta <- attr(pp, "petab_meta")
  scm  <- meta$sub_cond_map
  # condition columns: ids or names, whichever the importer used
  if ("conditionName" %in% names(conds)) {
    nm <- setNames(as.character(conds$conditionName), conds$conditionId)
    for (col in intersect(c("simulationConditionId", "preequilibrationConditionId"),
                          intersect(names(sim), names(scm)))) {
      v <- as.character(sim[[col]])
      swap <- !(v %in% scm[[col]]) & v %in% names(nm)
      v[swap] <- nm[v[swap]]
      sim[[col]] <- v
    }
  }
  # observable ids with a per-condition suffix
  obs <- names(meta$obs_meta$obs_trafo)
  bad <- !(sim$observableId %in% obs)
  sim$observableId[bad] <- sub("__.*$", "", sim$observableId[bad])
  trf  <- meta$obs_meta$obs_trafo
  cn   <- names(pp$dataList)
  times <- sort(unique(unlist(lapply(pp$dataList, function(d) d$time))))
  pred <- pp$prd(times, pp$bestfit, fixed = meta$fixed, deriv = FALSE)
  pre <- if ("preequilibrationConditionId" %in% names(sim)) sim$preequilibrationConditionId else NA
  pre[is.na(pre) | pre == "nan"] <- ""
  extra <- intersect(c("observableParameters", "noiseParameters"), intersect(names(scm), names(sim)))
  keyOf <- function(d, pre) do.call(paste, c(list(d$simulationConditionId, pre), lapply(extra, function(e) d[[e]])))
  key_sim <- keyOf(sim, pre)
  prep <- if ("preequilibrationConditionId" %in% names(scm)) scm$preequilibrationConditionId else ""
  prep[is.na(prep) | prep == "nan"] <- ""
  key_sub <- keyOf(scm, prep)
  rel <- rep(NA_real_, nrow(sim))
  for (i in seq_len(nrow(sim))) {
    k <- which(key_sub == key_sim[i])
    if (length(k) != 1L) next
    pc <- pred[[scm$sub_condition[k]]]
    if (is.null(pc) || !sim$observableId[i] %in% colnames(pc)) next
    tt <- if (is.infinite(sim$time[i])) max(pc[, "time"]) else sim$time[i]
    j  <- which(abs(pc[, "time"] - tt) < 1e-9)
    if (!length(j)) next
    v  <- pc[j[1], sim$observableId[i]]
    tr <- trf[[sim$observableId[i]]] %||% "lin"
    v  <- switch(tr, log = exp(v), log10 = 10^v, v)
    rel[i] <- abs(v - sim$simulation[i]) / max(abs(sim$simulation[i]), 1e-8)
  }
  rel
}

res <- list()
for (pr in probs) {
  cat("==", pr, "\n")
  yml <- file.path(bdir, pr, paste0(pr, ".yaml"))
  wd  <- file.path(tempdir(), pr); dir.create(wd, showWarnings = FALSE)
  r <- list(problem = pr, version = as.character(packageVersion("dMod2")))
  r$preeq <- {
    m <- utils::read.delim(list.files(file.path(bdir, pr), "^measurementData", full.names = TRUE)[1])
    "preequilibrationConditionId" %in% names(m) && any(!is.na(m$preequilibrationConditionId) &
                                                        nzchar(as.character(m$preequilibrationConditionId)))
  }
  r$t_import <- system.time(pp <- tryCatch(
    importPEtab(yml, backend = "cppDE", cores = 8L, outdir = wd,
                options = list(atol = 1e-12, rtol = 1e-10),
                optionsSens = list(atol = 1e-8, rtol = 1e-6)),
    error = function(e) conditionMessage(e)))[3]
  if (is.character(pp)) { r$error <- pp; res[[pr]] <- r; cat("  import failed:", pp, "\n"); next }
  r$npar <- length(pp$bestfit); r$ncond <- length(pp$dataList)
  r$t_obj <- system.time(v <- tryCatch(pp$obj(pp$bestfit), error = function(e) conditionMessage(e)))[3]
  r$value <- if (is.character(v)) NA_real_ else v$value
  if (is.character(v)) r$error <- v
  simf <- list.files(file.path(bdir, pr), "^simulatedData", full.names = TRUE)
  if (length(simf)) {
    cf <- list.files(file.path(bdir, pr), "^experimentalCondition", full.names = TRUE)[1]
    rel <- tryCatch(simCheck(pp, utils::read.delim(simf[1], stringsAsFactors = FALSE),
                             utils::read.delim(cf, stringsAsFactors = FALSE)),
                    error = function(e) { message("  sim check: ", conditionMessage(e)); NA })
    r$sim_n <- sum(!is.na(rel)); r$sim_med <- stats::median(rel, na.rm = TRUE)
    r$sim_max <- suppressWarnings(max(rel, na.rm = TRUE))
  }
  set.seed(1); fails <- 0L; msgs <- character(0)
  r$t_evals <- system.time(for (k in seq_len(nev)) {
    pk <- pp$bestfit + stats::rnorm(length(pp$bestfit), 0, 0.1)
    vk <- tryCatch(pp$obj(pk), error = function(e) conditionMessage(e))
    if (is.character(vk) || !is.finite(vk$value)) { fails <- fails + 1L; msgs <- c(msgs, substr(as.character(vk), 1, 150)) }
  })[3] / nev
  r$fails <- fails; r$fail_msg <- unique(msgs)
  cat(sprintf("  value %.4f  sim max rel %.2e (n=%d)  %.3f s/eval  fails %d/%d\n",
              r$value, r$sim_max %||% NA, r$sim_n %||% 0L, r$t_evals, fails, nev))
  res[[pr]] <- r
  saveRDS(res, out)
}
