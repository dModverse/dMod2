# The model of model.R at the best fit written as a PEtab problem, the state at rest as a
# pre-equilibration, and read back. Writes results/petab/.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- petabexport
# the estimated parameters at the best fit, the gauge as fixed values on the linear scale of the
# model
best <- as.parvec(readRDS("results/fits.rds"))
yaml <- exportPEtab(data, reactions, observables = g, p = p, pouter = best, errors = e,
                    fixed = exp(gauge), steadyState = pSS, events = events,
                    modelID = "oppelt", dir = "results/petab", overwrite = TRUE)
list.files("results/petab")

## ---- petabimport
petab <- importPEtab(yaml, backend = "cppDE", options = list(atol = 1e-10, rtol = 1e-8),
                     modelname = "oppelt_petab", compile = FALSE)

## ---- petabcompile
compile(petab$prd, petab$e, output = "oppelt_petab", cores = 6)

## ---- petabcompare
# the data term of model.R and the objective of the problem read back
c(model.R = attr(obj(c(best, gauge), deriv = FALSE), "data"),
  PEtab   = petab$obj(petab$bestfit, deriv = FALSE)$value)
