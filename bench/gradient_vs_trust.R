# Counts what a gradient-only optimiser asks for against what the trust region
# asks for, from identical starts, and projects both onto solver runs.
#
# Value and gradient must come from one call: `deriv = FALSE` integrates the
# plain system and `deriv = TRUE` the extended one, and their values differ at
# the integrator tolerance (~1e-9 relative on Boehm). A line search fed that
# inconsistent pair needs more evaluations, 308 against 198 here.
#
#   trust today          n_full * (1 + n_theta)
#   L-BFGS with adjoint  n_call * 2, one forward run shared with one reverse
#
# Measured with the existing forward sensitivities. No derivative code changes.

library(dMod2)

# OMP_NUM_THREADS above 1 fails internally. Parallelism goes through `cores`.
Sys.setenv(OMP_NUM_THREADS = 1)
cores <- 4
setwd(tempdir())

solverOpts <- list(atol = 1e-9, rtol = 1e-9)
startSd <- 0.1

## --- counting wrappers ------------------------------------------------------

# Three counters, because the three request kinds have different projected cost.
counter <- new.env()
resetCounter <- function() {
  counter$n <- 0L; counter$lastPar <- NULL; counter$last <- NULL
}

# One evaluation per distinct point, serving optim()'s separate fn and gr.
fusedFn <- function(obj, ...) function(p) {
  if (is.null(counter$lastPar) || !identical(unname(p), counter$lastPar)) {
    counter$last <- obj(p, deriv = TRUE, ...)
    counter$lastPar <- unname(p)
    counter$n <- counter$n + 1L
  }
  counter$last
}

projectRuns <- function(npar, what)
  switch(what, trust = counter$n * (1 + npar), lbfgs = counter$n * 2)

# Central differences on a subset of coordinates. Both optimisers stopping for
# bad reasons points at the gradient before it points at either of them.
# The gradient's name order need not match the parameter vector's, so it is
# indexed by name. Positional indexing silently pairs the wrong coordinates.
gradCheck <- function(obj, p, idx, h = 1e-5, ...) {
  g <- obj(p, deriv = TRUE, ...)$gradient
  nms <- names(p)[idx]
  fd <- vapply(idx, function(k) {
    pu <- pd <- p; pu[k] <- pu[k] + h; pd[k] <- pd[k] - h
    (obj(pu, deriv = FALSE, ...)$value -
     obj(pd, deriv = FALSE, ...)$value) / (2 * h)
  }, numeric(1))
  data.frame(par = nms, analytic = unname(g[nms]), fd = fd,
             relerr = abs(unname(g[nms]) - fd) / pmax(abs(fd), 1e-8),
             row.names = NULL)
}

## --- Boehm ------------------------------------------------------------------

petabBoehm <- importPEtab(
  system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2"),
  backend = "cppDE", cores = cores, outdir = tempdir(),
  optionsOde = solverOpts, optionsSens = solverOpts)

set.seed(1)
startBoehm <- pmin(pmax(petabBoehm$bestfit +
                          rnorm(length(petabBoehm$bestfit), sd = startSd),
                        petabBoehm$parlower + 1e-3),
                   petabBoehm$parupper - 1e-3)

gradCheck(petabBoehm$obj, startBoehm, seq_along(startBoehm))

resetCounter()
fnBoehm <- fusedFn(petabBoehm$obj)
trustBoehm <- trust(function(p) fnBoehm(p), startBoehm, rinit = 1, rmax = 10,
                    iterlim = 500, gtol = 1e-4,
                    parlower = petabBoehm$parlower,
                    parupper = petabBoehm$parupper)
nTrustBoehm <- counter$n
runsBoehmTrust <- projectRuns(length(startBoehm), "trust")

resetCounter()
fnBoehm <- fusedFn(petabBoehm$obj)
lbfgsBoehm <- optim(startBoehm,
                    fn = function(p) fnBoehm(p)$value,
                    gr = function(p) fnBoehm(p)$gradient,
                    method = "L-BFGS-B",
                    lower = petabBoehm$parlower, upper = petabBoehm$parupper,
                    control = list(maxit = 5000, factr = 1e1))
nLbfgsBoehm <- counter$n
runsBoehmLbfgs <- projectRuns(length(startBoehm), "lbfgs")

## --- Bachmann ---------------------------------------------------------------

yamlBachmann <- path.expand(file.path(
  "~/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models",
  "Bachmann_MSB2011/Bachmann_MSB2011.yaml"))

petabBachmann <- importPEtab(yamlBachmann, backend = "cppDE",
                             cores = cores, outdir = tempdir(),
                             optionsOde = solverOpts, optionsSens = solverOpts)

set.seed(1)
startBachmann <- pmin(pmax(petabBachmann$bestfit +
                             rnorm(length(petabBachmann$bestfit), sd = startSd),
                           petabBachmann$parlower + 1e-3),
                      petabBachmann$parupper - 1e-3)

set.seed(2)
gradCheck(petabBachmann$obj, startBachmann,
          sort(sample(seq_along(startBachmann), 12)), cores = cores)

resetCounter()
fnBachmann <- fusedFn(petabBachmann$obj, cores = cores)
trustBachmann <- trust(function(p) fnBachmann(p), startBachmann,
                       rinit = 1, rmax = 10, iterlim = 500, gtol = 1e-4,
                       parlower = petabBachmann$parlower,
                       parupper = petabBachmann$parupper)
nTrustBachmann <- counter$n
runsBachmannTrust <- projectRuns(length(startBachmann), "trust")

resetCounter()
fnBachmann <- fusedFn(petabBachmann$obj, cores = cores)
lbfgsBachmann <- optim(startBachmann,
                       fn = function(p) fnBachmann(p)$value,
                       gr = function(p) fnBachmann(p)$gradient,
                       method = "L-BFGS-B",
                       lower = petabBachmann$parlower,
                       upper = petabBachmann$parupper,
                       control = list(maxit = 5000, factr = 1e1))
nLbfgsBachmann <- counter$n
runsBachmannLbfgs <- projectRuns(length(startBachmann), "lbfgs")


## --- Lang, 294 parameters ---------------------------------------------------
## Needs the observableFormula rename fix; without it importPEtab fails on the
## leading-underscore species ids.

yamlLang <- path.expand(file.path(
  "~/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models",
  "Lang_PLOSComputBiol2024/Lang_PLOSComputBiol2024.yaml"))

petabLang <- importPEtab(yamlLang, backend = "cppDE",
                         cores = cores, outdir = tempdir(),
                         optionsOde = solverOpts, optionsSens = solverOpts)

set.seed(1)
startLang <- pmin(pmax(petabLang$bestfit +
                         rnorm(length(petabLang$bestfit), sd = startSd),
                       petabLang$parlower + 1e-3),
                  petabLang$parupper - 1e-3)

resetCounter()
fnLang <- fusedFn(petabLang$obj, cores = cores)
trustLang <- trust(function(p) fnLang(p), startLang, rinit = 1, rmax = 10,
                   iterlim = 500, gtol = 1e-4,
                   parlower = petabLang$parlower,
                   parupper = petabLang$parupper)
nTrustLang <- counter$n
runsLangTrust <- projectRuns(length(startLang), "trust")

resetCounter()
fnLang <- fusedFn(petabLang$obj, cores = cores)
lbfgsLang <- optim(startLang,
                   fn = function(p) fnLang(p)$value,
                   gr = function(p) fnLang(p)$gradient,
                   method = "L-BFGS-B",
                   lower = petabLang$parlower, upper = petabLang$parupper,
                   control = list(maxit = 5000, factr = 1e1))
nLbfgsLang <- counter$n
runsLangLbfgs <- projectRuns(length(startLang), "lbfgs")

## --- Verdict ----------------------------------------------------------------
## Projected solver runs, and whether L-BFGS reaches a comparable objective.

## `valOpt` is the objective at the published estimate: a run that stops far
## short of it makes the evaluation counts a comparison of two failures.

data.frame(
  model     = c("Boehm", "Bachmann", "Lang"),
  npar      = c(length(startBoehm), length(startBachmann), length(startLang)),
  nTrust    = c(nTrustBoehm, nTrustBachmann, nTrustLang),
  nLbfgs    = c(nLbfgsBoehm, nLbfgsBachmann, nLbfgsLang),
  runsTrust = c(runsBoehmTrust, runsBachmannTrust, runsLangTrust),
  runsLbfgs = c(runsBoehmLbfgs, runsBachmannLbfgs, runsLangLbfgs),
  speedup   = c(runsBoehmTrust / runsBoehmLbfgs,
                runsBachmannTrust / runsBachmannLbfgs,
                runsLangTrust / runsLangLbfgs),
  valTrust  = c(trustBoehm$value, trustBachmann$value, trustLang$value),
  valLbfgs  = c(lbfgsBoehm$value, lbfgsBachmann$value, lbfgsLang$value),
  valOpt    = c(petabBoehm$obj(petabBoehm$bestfit)$value,
                petabBachmann$obj(petabBachmann$bestfit, cores = cores)$value,
                petabLang$obj(petabLang$bestfit, cores = cores)$value))
