# The same problem from its PEtab files, compared with model.R at the published values, and
# model.R written as a PEtab problem.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- import
petabDir <- system.file("extdata", "petab_bachmann", package = "dMod2")
petab <- importPEtab(file.path(petabDir, "Bachmann_MSB2011.yaml"), backend = "cppDE",
                     options = solverOptions, modelname = "bachmann_petab")

## ---- compare
# the values of the parameter table are the published optimum
published <- petab$bestfit[outer]
c(script = attr(obj(published, deriv = FALSE), "chi2"),
  PEtab  = attr(petab$obj(published, deriv = FALSE), "chi2"))

## ---- export
dir.create("results", showWarnings = FALSE)
# PEtab version 2 has no parameter scales: the values are written on the linear scale
exportPEtab(data, reactions, observables = g, p = p, pouter = published, errors = e,
            dir = "results/petab", overwrite = TRUE)
list.files("results/petab")
