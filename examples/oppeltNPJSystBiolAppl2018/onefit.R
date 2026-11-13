# One fit from a start near the best fit of fit.R, with the gauge held fixed. Saves
# results/onefit.rds, which the guide reads.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- onefit
# the start drawn with standard deviation 0.5 around the best fit, and drawn again where the
# model cannot be evaluated
best <- as.parvec(readRDS("results/fits.rds"))
set.seed(1)
fit <- mstrust(obj, best, fixed = gauge, sd = 0.5, fits = 1, rinit = 0.1, rmax = 10,
               iterlim = 1000, tolControl = tolerances)[[1]]

## ---- onefitsave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "onefit.R", "data/oppelt2018.csv", "results/fits.rds")
attr(fit, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(fit, "results/onefit.rds")
