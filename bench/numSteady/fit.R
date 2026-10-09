# numSteady bench: fit of the minimal model with steady states from Pimpl.
#
#   Rscript bench/numSteady/fit.R [nfits] [cores] [outdir]
#
# prd = g * x * pSS * pCond. pCond maps log10 outer parameters to conditions,
# pSS computes the steady state per condition, x integrates with the events.
# Reports the objective, one trust fit, a multistart and the Pimpl statistics
# (calls, solves, memo hits, iterations per solve, successful starts).

args   <- commandArgs(TRUE)
nfits  <- if (length(args) >= 1L) as.integer(args[1]) else 16L
cores  <- if (length(args) >= 2L) as.integer(args[2]) else 8L
outdir <- if (length(args) >= 3L) args[3] else file.path(tempdir(), "numSteady")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(dMod2))
here <- if (file.exists("bench/numSteady/model.R")) "bench/numSteady" else "."
source(file.path(here, "model.R"))
m <- numSteadyModel(here)
setwd(outdir)

# ---- model ------------------------------------------------------------------
x <- odemodel(m$reactions, events = m$events, modelname = "ns_x", compile = FALSE) |>
  Xs(options = list(atol = 1e-10, rtol = 1e-8), optionsSens = list(atol = 1e-6, rtol = 1e-4))
g <- Y(m$observables, f = x, attachInput = FALSE, compile = FALSE, modelname = "ns_g")
pSS <- Pimpl(m$reactions, forcings = m$forcings, compile = FALSE, modelname = "ns_ss")

# condition layer: inputs of pSS and x, log10 outer parameters
inner <- union(getParameters(pSS), setdiff(getParameters(x, g), m$reactions$states))
inner <- union(inner, m$forcings)
trafo <- eqnvec() |>
  define("x~x", x = inner) |>
  branch(table = m$grid, conditions = rownames(m$grid), apply = "insert") |>
  insert("x~0", x = m$forcings) |>
  insert("volumeC ~ 1", "volumeEC ~ 10") |>
  insert("x~exp10(x)", x = .currentSymbols[!grepl("^offset_", .currentSymbols)])
pCond <- P(trafo, compile = FALSE, modelname = "ns_p")

compile(g, x, pSS, pCond, output = "ns_all", cores = cores, verbose = FALSE)
prd <- g * x * pSS * pCond

outer <- getParameters(prd)
cat(length(outer), "outer parameters,", length(m$data), "conditions,",
    sum(vapply(m$data, nrow, 1L)), "data points\n")
center <- setNames(rep(-1, length(outer)), outer)
center[grepl("^offset_", outer)] <- 0
center[grepl("^tSmad", outer)] <- log10(0.4)

obj <- normL2(m$data, prd, t0 = m$t0) +
  constraintL2(center, sigma = 3, attrName = "prior")

ss <- environment(attr(pSS, "mappings")[[1]])$stats
showStats <- function(label) {
  how <- if (length(ss$how)) paste(names(ss$how), ss$how, sep = "=", collapse = " ") else "-"
  cat(sprintf("%-10s calls %6d  solves %5d  memo %6d  iter/solve %5.1f  failed %d  starts: %s\n",
              label, ss$calls, ss$solves, ss$memo, ss$iter / max(ss$solves, 1L), ss$failed, how))
}

t1 <- system.time(v1 <- obj(center))[3]
t2 <- system.time(v2 <- obj(center))[3]
cat(sprintf("objective at the centre %.2f; first call %.2f s, repeated %.2f s\n", v1$value, t1, t2))
showStats("centre")

# ---- single trust fit, statistics of this process -----------------------------
statsReset <- environment(attr(pSS, "mappings")[[1]])$statsReset
statsReset()
set.seed(1)
t3 <- system.time(f1 <- trust(obj, center + rnorm(length(center), 0, 1), rinit = 0.1, rmax = 10,
                              iterlim = 400, hessianFallback = "sr1", fallbackLimit = 1L))[3]
cat(sprintf("\none trust fit: value %.2f, converged %s, %d iterations, %.0f s\n",
            f1$value, f1$converged, f1$iterations, t3))
showStats("trust")

# ---- multistart ----------------------------------------------------------------
set.seed(2)
tm <- system.time(fits <- mstrust(obj, center, studyname = "numSteady", fits = nfits, cores = cores,
                                  sd = 1, iterlim = 400, output = FALSE,
                                  hessianFallback = "sr1", fallbackLimit = 1L))[3]
pf <- as.parframe(fits)
ok <- vapply(fits, function(f) is.list(f) && is.null(f$error), TRUE)
cat(sprintf("\n%d of %d fits returned, %d converged, %.0f s\n",
            sum(ok), nfits, sum(pf$converged), tm))
print(head(pf[order(pf$value), c("value", "converged", "iterations")], 10))
saveRDS(list(fits = fits, center = center), file.path(outdir, "fits.rds"))
