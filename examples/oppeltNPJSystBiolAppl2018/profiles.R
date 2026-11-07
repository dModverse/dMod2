# Profile likelihood of every parameter of the full model, at the best fit of fit.R refitted with
# a wide prior. Saves results/profiles.rds, which plots.R and reduced.R read.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- wideprior
# the prior widened to standard deviation 10, so that the profiles follow the data; the best fit
# of the multistart refitted with it
objWide <- normL2(data, prd, e) + constraintL2(mu, sigma = 10) + stability
best <- as.parvec(readRDS("results/fits.rds"))
refit <- trust(objWide, best, fixed = gauge, rinit = 0.1, rmax = 10, iterlim = 1000,
               tolControl = tolerances)
bestWide <- refit$argument

## ---- profiles
# every step re-optimised by up to 100 iterations of trust(); a profile stops where the data term
# has risen by the 95% threshold, or 7 log units from the optimum. Where a profile reaches a lower
# value than the fit, the fit restarts from there and all profiles run again.
repeat {
  profiles <- profile(objWide, bestWide, names(bestWide), fixed = gauge,
                      limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
                      algoControl = list(reoptimize = TRUE),
                      optControl = list(iterlim = 100, tolControl = tolerances),
                      cores = getOption("dMod2.cores"))
  if (min(profiles$value) > refit$value - 1e-2) break
  refit <- trust(objWide, as.parvec(profiles)[names(bestWide)], fixed = gauge, rinit = 0.1,
                 rmax = 10, iterlim = 1000, tolControl = tolerances)
  bestWide <- refit$argument
}

dir.create("results", showWarnings = FALSE)
saveRDS(list(refit = refit, profiles = profiles), "results/profiles.rds")
