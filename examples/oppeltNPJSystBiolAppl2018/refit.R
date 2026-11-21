# Fit and profiles of the reduced model of reduced.R, and the likelihood ratio against the full
# model. Saves results/reduced.rds, which plots.R and prediction.R read.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- reducedfit
# from the full model refitted in profiles.R, the basal phosphorylation in the pNFκB complex at
# the product of the two rates it replaces, then a multistart around that fit
full <- readRDS("results/profiles.rds")$refit
start <- c(full$argument, gauge, act_Ikb_complex = full$argument[["act_Ikb_Base"]] +
             full$argument[["eta_act_Ikb_Base_c_pNfkIkb"]])[outerReduced]
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

## ---- reducedprofileprior
# the prior of the profiles of profiles.R, standard deviation 20 at the best fit, and one fit from
# there
objReducedProfile <- normL2(data, prdReduced, eReduced) + constraintL2(bestReduced, sigma = 20)
fitReducedProfile <- trust(objReducedProfile, bestReduced, fixed = gaugeReduced, rinit = 0.1,
                           rmax = 10, iterlim = 1000, tolControl = tolerances)
bestReducedProfile <- fitReducedProfile$argument

## ---- reducedprofiles
# the settings of the full model
profilesReduced <- profile(objReducedProfile, bestReducedProfile, names(bestReducedProfile),
                           fixed = gaugeReduced, limits = c(lower = -7, upper = 7),
                           stepControl = list(stop = "data"), algoControl = list(reoptimize = TRUE),
                           optControl = list(iterlim = 100, tolControl = tolerances),
                           cores = getOption("dMod2.cores"))

## ---- reducedsave
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "reduced.R", "refit.R", "data/oppelt2018.csv", "results/profiles.rds")
saved <- list(fits = fitsReduced, ratio = ratio, fitProfile = fitReducedProfile,
              profiles = profilesReduced)
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/reduced.rds")
