# Multistart of the full model with its symmetries removed by the gauge of model.R. Saves
# results/fits.rds, which plots.R and profiles.R read.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- multistart
# 500 fits from starts around the means of the prior, standard deviation 2 on the log scale; a
# start at which the model cannot be evaluated is replaced by the next draw
free <- setdiff(names(mu), names(gauge))
set.seed(1)
fits <- mstrust(obj, mu[free], fixed = gauge, sd = 2, fits = 500, iterlim = 1000,
                tolControl = tolerances, cores = getOption("dMod2.cores")) |>
  as.parframe()

dir.create("results", showWarnings = FALSE)
saveRDS(fits, "results/fits.rds")
