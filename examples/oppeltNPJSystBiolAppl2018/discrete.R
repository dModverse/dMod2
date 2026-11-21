# Discrete symmetries of the reduced model of reduced.R, by symident: isolated parameter values with
# the output of a given one, and the chart of one side of the fixed set. Saves results/discrete.rds,
# which the guide reads.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- discrete
# the settings of symmetries.R; the resting state and the fixed values of reduced.R enter as
# initial values and parameters of the detection
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
unbleachedReduced <- setNames(eqnvec(sub(" - tau_.*$", "", observablesReduced)),
                              names(observablesReduced))
known <- c(fixedReduced[names(fixedReduced) != "tau_A20GFP"], resting[names(resting) != "TNF"])
discrete <- symmetryDetection(reactionsReduced, unbleachedReduced,
                              trafo = setNames(as.character(known), names(known)), events = events,
                              forcings = "TNF", conditions = doses, mode = "global",
                              cores = getOption("dMod2.cores"))

## ---- fundamental
# the side degrad_A20 < degrad_RnaA20 of the fixed set, in a chart of positive coordinates
fundamental <- symmetryReduction(discrete, branch = "degrad_A20 < degrad_RnaA20")

## ---- discretesave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "reduced.R", "discrete.R")
saved <- list(discrete = discrete, fundamental = fundamental)
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/discrete.rds")
