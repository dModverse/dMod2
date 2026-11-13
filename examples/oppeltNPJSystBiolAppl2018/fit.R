# Multistart of the full model with the symmetries removed by the gauge of model.R. Saves
# results/fits.rds, which plots.R and profiles.R read.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- multistart
# 200 fits from starts drawn from the prior; a start at which the model cannot be evaluated is
# replaced by the next draw. helix.R runs the same fits as 100 array tasks of 2 fits each.
free <- setdiff(names(mu), names(gauge))
set.seed(1)
fits <- mstrust(obj, mu[free], fixed = gauge, sd = sigma[free], fits = 200, iterlim = 1000,
                tolControl = tolerances, cores = getOption("dMod2.cores")) |>
  as.parframe()

## ---- fitsave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "fit.R", "data/oppelt2018.csv")
attr(fits, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(fits, "results/fits.rds")
