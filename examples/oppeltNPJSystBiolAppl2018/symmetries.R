# Structural identifiability of the full model: the Lie symmetries of the parameters, by
# symident, and their reduction. Saves results/symmetries.rds and results/reduction.rds, which the
# guide reads.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- symmetries
# at rest without TNFα, which an event adds at time 0 in the stimulated condition; the moieties
# of NFκB and IKK by their totals, the parameters of Pimpl()
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
# without the bleaching terms: in the resting control the states are constant, and the slope of
# an imaging readout is the bleaching rate alone
unbleached <- setNames(eqnvec(sub(" - tau_.*$", "", observables)), names(observables))
symmetries <- symmetryDetection(reactions, unbleached, events = events, forcings = "TNF",
                                conditions = doses, equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
                                cores = getOption("dMod2.cores"))

## ---- reduction
# a section of every orbit of the scalings
reduction <- symmetryReduction(symmetries, timeout = 120)

## ---- symmetriessave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "symmetries.R", "data/oppelt2018.csv")
attr(symmetries, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(symmetries, "results/symmetries.rds")
attr(reduction, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(reduction, "results/reduction.rds")
