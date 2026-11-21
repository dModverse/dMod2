# The reduced model of reduced.R in the chart of the fundamental domain degrad_A20 < degrad_RnaA20
# of the discrete symmetry of discrete.R: one fit from the best fit of refit.R and the profiles of
# all parameters. Saves results/fundamental.rds, which the guide reads.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- fundamentaltrafo
# the parameters of reduced.R with degrad_RnaA20 above degrad_A20 by exp(delta_RnaA20), the chart
# of symmetryReduction(branch = "degrad_A20 < degrad_RnaA20"); ODE, observables and errors of
# reduced.R
trafoFundamental <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(xReduced, gReduced, eReduced),
                               reactionsReduced$states)) |>
  define("x <- y", x = names(fixedReduced), y = fixedReduced) |>
  define("x <- y", x = names(resting), y = resting) |>
  define("degrad_RnaA20 <- degrad_A20 + delta_RnaA20") |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(data)) |>
  define("dose_TNF <- 0", conditionMatch = "control") |>
  define("dose_TNF <- 1", conditionMatch = "TNF")
pFundamental <- P(trafoFundamental, modelname = "oppelt_fund_par", compile = FALSE)
compile(pFundamental, output = "oppelt_fund", cores = getOption("dMod2.cores"))
prdFundamental <- gReduced * xReduced * pFundamental

## ---- fundamentalobjective
# the prior of reduced.R; the difference of the two rates with the prior of degrad_RnaA20
outerFundamental <- getParameters(prdFundamental)
muFundamental <- c(muReduced, delta_RnaA20 = muReduced[["degrad_RnaA20"]])[outerFundamental]
sigmaFundamental <- c(sigmaReduced,
                      delta_RnaA20 = sigmaReduced[["degrad_RnaA20"]])[outerFundamental]
objFundamental <- normL2(data, prdFundamental, eReduced) +
  constraintL2(muFundamental, sigma = sigmaFundamental)

## ---- fundamentalfit
# the best fit of refit.R, on the side degrad_A20 < degrad_RnaA20, in the chart
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
ordered <- sort(exp(unclass(bestReduced)[c("degrad_A20", "degrad_RnaA20")]))
start <- unclass(bestReduced)[setdiff(names(bestReduced), "degrad_RnaA20")]
start[c("degrad_A20", "delta_RnaA20")] <- log(c(ordered[[1]], ordered[[2]] - ordered[[1]]))
fitFundamental <- trust(objFundamental, start[setdiff(outerFundamental, names(gaugeReduced))],
                        fixed = gaugeReduced, rinit = 0.1, rmax = 10, iterlim = 1000,
                        tolControl = tolerances)
# the data terms of the best fit of refit.R, of the same point in the chart, and of the fit
dataFundamental <- normL2(data, prdFundamental, eReduced)
dataTermReduced <- normL2(data, prdReduced, eReduced)
dataTerms <- c(reduced = dataTermReduced(bestReduced, fixed = gaugeReduced)$value,
               chart = dataFundamental(start, fixed = gaugeReduced)$value,
               fit = dataFundamental(fitFundamental$argument, fixed = gaugeReduced)$value)

## ---- fundamentalprofileprior
# the prior of the profiles of refit.R, standard deviation 20 at the fit, and one fit from there
objFundamentalProfile <- normL2(data, prdFundamental, eReduced) +
  constraintL2(fitFundamental$argument, sigma = 20)
fitFundamentalProfile <- trust(objFundamentalProfile, fitFundamental$argument, fixed = gaugeReduced,
                               rinit = 0.1, rmax = 10, iterlim = 1000, tolControl = tolerances)
bestFundamentalProfile <- fitFundamentalProfile$argument

## ---- fundamentalprofiles
# the settings of refit.R
profilesFundamental <- profile(objFundamentalProfile, bestFundamentalProfile,
                               names(bestFundamentalProfile), fixed = gaugeReduced,
                               limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
                               algoControl = list(reoptimize = TRUE),
                               optControl = list(iterlim = 100, tolControl = tolerances),
                               cores = getOption("dMod2.cores"))

## ---- fundamentalsave
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "reduced.R", "fundamental.R", "data/oppelt2018.csv")
saved <- list(fit = fitFundamental, dataTerms = dataTerms, fitProfile = fitFundamentalProfile,
              profiles = profilesFundamental)
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/fundamental.rds")
