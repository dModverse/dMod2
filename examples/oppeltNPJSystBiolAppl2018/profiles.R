# Profile likelihood of every parameter of the full model, at the best fit of fit.R. Saves
# results/profiles.rds, which plots.R, refit.R and variants.R read.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- refit
# the best fit of the multistart, refitted to the tolerances of the profiles
refit <- trust(obj, as.parvec(readRDS("results/fits.rds")), fixed = gauge, rinit = 0.1, rmax = 10,
               iterlim = 1000, tolControl = tolerances)

## ---- profileprior
# the profiles replace the prior of model.R by one of standard deviation 20 at the refit, so that
# the data term shapes them; they start at the optimum of that objective
objProfile <- normL2(data, prd, e) + constraintL2(refit$argument, sigma = 20) + stability
# the best of 20 fits from `center` and from starts about it, sd 0.1, which reaches the floor of
# a flat valley where a single fit stops on it
bestFit <- function(center) {
  fits <- mstrust(objProfile, center, fixed = gauge, sd = 0.1, fits = 20, startFromCenter = TRUE,
                  rinit = 0.1, rmax = 10, iterlim = 1000, tolControl = tolerances,
                  cores = getOption("dMod2.cores"))
  values <- vapply(fits, function(f) if (is.numeric(f$value)) f$value else Inf, 0)
  fits[[which.min(values)]]
}
fitProfile <- bestFit(refit$argument)
best <- fitProfile$argument

## ---- profiles
# the profiles integrated, every step re-optimized, each side until the data term has risen by the
# threshold or the parameter has moved by 7 on the natural log scale, the settings of refit.R.
# Where a profile reaches a lower value than the fit, the fit restarts from there and all
# profiles run again. helix.R runs one profile per array task.
limits <- c(lower = -7, upper = 7)
repeat {
  profiles <- profile(objProfile, best, names(best), fixed = gauge, limits = limits,
                      stepControl = list(stop = "data"), algoControl = list(reoptimize = TRUE),
                      optControl = list(iterlim = 100, tolControl = tolerances),
                      cores = getOption("dMod2.cores"))
  if (min(profiles$value) > fitProfile$value - 1e-2) break
  fitProfile <- bestFit(as.parvec(profiles)[names(best)])
  best <- fitProfile$argument
}

## ---- profilessave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "profiles.R", "data/oppelt2018.csv", "results/fits.rds")
saved <- list(refit = refit, fitProfile = fitProfile, profiles = profiles)
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/profiles.rds")
