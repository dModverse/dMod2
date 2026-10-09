# numSteady bench: Pimpl against the symbolic steady state of steadyStates().
#
#   Rscript bench/numSteady/exact.R [n per arm] [outdir]
#
# Per arm the perturbation parameters are substituted and steadyStates()
# parametrises the steady states of the model without ligand states (they are 0
# at steady state). Random values of the free symbols (rate constants, states,
# flux ratios r_*) give exact steady states with their parameters. Pimpl, built
# on the full model, gets the same parameters and the totals of the exact state
# and has to return it, cold (all states start at 1) and started at the exact
# state. A cold result that differs is checked to be a stable steady state of
# its own (a second basin). Exact states that are not stable are counted apart.

args   <- commandArgs(TRUE)
nper   <- if (length(args) >= 1L) as.integer(args[1]) else 40L
outdir <- if (length(args) >= 2L) args[2] else file.path(tempdir(), "numSteady_exact")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(dMod2); library(data.table)})
here <- if (file.exists("bench/numSteady/model.R")) "bench/numSteady" else "."
source(file.path(here, "model.R"))
full <- numSteadyModel(here)
bare <- numSteadyModel(here, ligand = FALSE)
setwd(outdir)

pss <- Pimpl(full$reactions, forcings = full$forcings, compile = TRUE, modelname = "nse_ss",
             keepRoot = FALSE)
E   <- environment(attr(pss, "mappings")[[1]])
dep <- E$dependent

arms <- lapply(list(Ctrl = c(0, 0, 0, 0), R1Knd = c(0, 0, 1, 0), R2Knd = c(0, 0, 0, 1),
                    OE_R1 = c(1, 0, 0, 0), OE_R2 = c(0, 1, 0, 0)), setNames, full$perturbation)

symbolic <- function(pert) {
  r <- bare$reactions
  r$rates <- cOde::replaceSymbols(names(pert), as.character(pert), r$rates)
  for (v in c("1.4", "1.3")) {
    ss <- tryCatch(suppressMessages({
      utils::capture.output(y <- steadyStates(r, forcings = bare$forcings, version = v, verbose = FALSE))
      y }), error = function(e) NULL)
    if (is.character(ss) && length(ss) > 1L) return(list(ss = ss, version = v))
  }
  NULL
}

evalSS <- function(ss, free) {
  env <- as.list(free); out <- setNames(rep(NA_real_, length(ss)), names(ss))
  for (sweep in seq_along(ss)) {
    todo <- names(out)[is.na(out)]
    if (!length(todo)) break
    for (k in todo) {
      v <- tryCatch(eval(parse(text = ss[[k]]), env), error = function(e) NA_real_)
      if (length(v) == 1L && !is.na(v)) { out[k] <- v; env[[k]] <- v }
    }
  }
  out
}

set.seed(1); res <- list()
for (a in names(arms)) {
  sym <- symbolic(arms[[a]])
  if (is.null(sym)) { cat(a, ": steadyStates() found no parametrisation\n"); next }
  ss <- sym$ss
  # free: symbols that are not solved for; a free state maps to itself
  solved   <- names(ss)[unlist(ss) != names(ss)]
  freeSyms <- setdiff(cOde::getSymbols(unlist(ss)), c(solved, bare$forcings))
  for (i in seq_len(nper)) {
    free <- setNames(10^runif(length(freeSyms), -2, 1), freeSyms)
    free[intersect(freeSyms, c("volumeC", "volumeEC"))] <- c(volumeC = 1, volumeEC = 10)[intersect(freeSyms, c("volumeC", "volumeEC"))]
    val <- evalSS(ss, free)
    vals <- c(free, val)
    if (anyNA(val) || any(!is.finite(val)) || any(val < 0)) {
      res[[length(res) + 1L]] <- data.table(arm = a, version = sym$version, valid = FALSE); next
    }
    x <- setNames(rep(0, length(dep)), dep)
    x[intersect(dep, names(vals))] <- vals[intersect(dep, names(vals))]
    p <- setNames(rep(1, length(E$parms_all)), E$parms_all)
    p[intersect(E$parms_all, names(vals))] <- vals[intersect(E$parms_all, names(vals))]
    p[names(arms[[a]])] <- arms[[a]]
    p[c("volumeC", "volumeEC")[c("volumeC", "volumeEC") %in% names(p)]] <- c(1, 10)[c("volumeC", "volumeEC") %in% names(p)]
    p[rownames(E$C_dep)] <- drop(E$C_dep %*% x)
    act <- dep[x > 0]
    st  <- E$stabilityOf(x, p, act)
    stable <- st$maxRe <= 1e-8 * st$scale
    big <- x > 1e-12 * max(x)
    relErr <- function(o) {
      if (is.character(o)) return(NA_real_)
      xo <- unclass(o)[dep]
      max(abs(xo[big] - x[big]) / x[big], abs(xo[!big] - x[!big]) / max(x))
    }
    t <- system.time(o <- tryCatch(pss(p, deriv = TRUE)[[1]], error = function(e) conditionMessage(e)))[3]
    ok <- !is.character(o)
    err <- relErr(o)
    # a differing cold root: a steady state of its own and stable?
    other <- NA
    if (ok && err > 1e-6) {
      xo <- unclass(o)[dep]; pv <- p
      f  <- E$eval_G(xo, pv)[1, dep]
      J  <- E$eval_J(xo, pv)[dep, dep]
      rel <- max(abs(f) / pmax(abs(J) %*% abs(xo), 1e-300)[, 1][xo > 0])
      sto <- E$stabilityOf(xo, pv, dep[xo > 0])
      other <- rel < 1e-6 && sto$maxRe <= 1e-8 * sto$scale
    }
    ow <- tryCatch(pss(c(p, x), deriv = FALSE)[[1]], error = function(e) conditionMessage(e))
    res[[length(res) + 1L]] <- data.table(arm = a, version = sym$version, valid = TRUE, stable = stable,
                                          ok = ok, time = t, relerr = err, other_stable = other,
                                          relerr_warm = relErr(ow),
                                          msg = if (ok) "" else substr(o, 1, 100))
  }
}
res <- rbindlist(res, fill = TRUE)
saveRDS(res, file.path(outdir, "exact.rds"))
st <- res[stable %in% TRUE]
print(st[, .(steadyStates = version[1], stable = .N, cold_equal = sum(relerr <= 1e-6, na.rm = TRUE),
             cold_other_basin = sum(other_stable %in% TRUE), cold_wrong = sum(other_stable %in% FALSE),
             cold_failed = sum(!ok), cold_relerr_med = median(relerr, na.rm = TRUE),
             warm_relerr_max = max(relerr_warm, na.rm = TRUE), t_med = median(time)), by = arm])
cat("exact states that are not stable:", sum(res$stable %in% FALSE), "\n")
