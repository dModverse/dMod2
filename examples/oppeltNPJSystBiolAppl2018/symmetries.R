# Structural identifiability of the full model: the Lie symmetries of its parameters, by
# symident, and their reduction. Saves results/symmetries.rds and results/reduction.rds, which the
# guide reads.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- symmetries
# at rest without TNFα, which an event adds at time 0 in the stimulated condition; the moieties
# of NFκB and IKK by their totals, the parameters of Pimpl()
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
# without the bleaching terms: in the resting control the states are constant, and the slope of
# an imaging readout is its bleaching rate alone
unbleached <- eqnvec(sub(" - tau_.*$", "", observables), names = names(observables))
symmetries <- symmetryDetection(reactions, unbleached, events = events, forcings = "TNF",
                                conditions = doses, equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
                                cores = getOption("dMod2.cores"))

## ---- reduction
# the invariants of the general direction and a section of every orbit
reduction <- symmetryReduction(symmetries, timeout = 120)

dir.create("results", showWarnings = FALSE)
saveRDS(symmetries, "results/symmetries.rds")
saveRDS(reduction, "results/reduction.rds")
