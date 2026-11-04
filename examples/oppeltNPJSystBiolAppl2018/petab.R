# The model of model.R at its best fit written as a PEtab problem, the state at rest as a
# pre-equilibration, and read back. Writes results/petab/.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- petabexport
# the estimated parameters at the best fit, the gauge as fixed values on the scale of the model:
# K_R is the one parameter not on the log scale
best <- as.parvec(readRDS("results/fits.rds"))
linear <- replace(exp(gauge), "K_R", gauge[["K_R"]])
yaml <- exportPEtab(data, reactions, observables = g, p = p, pouter = best, errors = e,
                    fixed = linear, steadyState = pSS, events = events,
                    modelID = "oppelt", dir = "results/petab", overwrite = TRUE)
list.files("results/petab")

## ---- petabimport
petab <- importPEtab(yaml, backend = "cppDE", options = list(atol = 1e-10, rtol = 1e-8),
                     modelname = "oppelt_petab")
# the data term of model.R and the objective of the problem read back
c(model.R = attr(obj(c(best, gauge), deriv = FALSE), "data"),
  PEtab   = petab$obj(petab$bestfit, deriv = FALSE)$value)
