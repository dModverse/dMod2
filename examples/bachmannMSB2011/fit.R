# Multistart: 200 fits from starts drawn uniformly within the bounds. Saves results/fits.rds,
# which plots.R reads. Takes about 100 core hours.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- starts
starts <- msParframe(lower, n = 200, seed = 1, samplefun = runif, keepfirst = FALSE,
                     min = lower, max = upper)

## ---- multistart
fits <- mstrust(obj, starts, cores = getOption("dMod2.cores"), iterlim = 5000,
                parlower = lower, parupper = upper) |>
  as.parframe()

dir.create("results", showWarnings = FALSE)
saveRDS(fits, "results/fits.rds")
