# Fit and profiles of the reduced model of reduced.R, and its likelihood ratio against the full
# model. Saves results/reduced.rds, which plots.R and prediction.R read.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- reducedfit
# from the full model refitted in profiles.R, then a multistart around that fit
full <- readRDS("results/profiles.rds")$refit
start <- c(full$argument, gauge)[outerReduced]
freeReduced <- setdiff(outerReduced, names(gaugeReduced))
fit <- trust(objReduced, start[freeReduced], fixed = gaugeReduced, rinit = 0.1, rmax = 10,
             iterlim = 1000, tolControl = tolerances)
set.seed(1)
fitsReduced <- mstrust(objReduced, fit$argument, fixed = gaugeReduced, sd = 1, fits = 50,
                       startFromCenter = TRUE, iterlim = 1000, tolControl = tolerances,
                       cores = getOption("dMod2.cores")) |>
  as.parframe()
bestReduced <- as.parvec(fitsReduced)

## ---- likelihoodratio
# the data terms of the two models, the difference against the 95% quantile of chi-square with
# the number of parameters the reduction removed
dataFull <- normL2(data, prd, e)(full$argument, fixed = gauge)$value
dataReduced <- normL2(data, prdReduced, eReduced)(bestReduced, fixed = gaugeReduced)$value
removed <- length(full$argument) - length(bestReduced)
ratio <- c(full = dataFull, reduced = dataReduced, difference = dataReduced - dataFull,
           removed = removed, threshold = qchisq(0.95, removed))

## ---- reducedprofiles
profilesReduced <- profile(objReduced, bestReduced, names(bestReduced), fixed = gaugeReduced,
                           limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
                           algoControl = list(reoptimize = TRUE),
                           optControl = list(iterlim = 100, tolControl = tolerances),
                           cores = getOption("dMod2.cores"))

saveRDS(list(fits = fitsReduced, ratio = ratio, profiles = profilesReduced), "results/reduced.rds")
