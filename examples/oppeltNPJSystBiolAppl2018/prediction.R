# Prediction band of nuclear IκBα after TNFα, which no experiment measured, from prediction
# profiles of the reduced model. Saves results/prediction.rds, which plots.R reads.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- predictionsetup
# the readouts of the reduced model and the logarithm of nuclear IκBα
gPrediction <- Y(c(observablesReduced, lognIkb = "log(nIkb)"), xReduced,
                 modelname = "oppelt_red_prediction", compile = FALSE)
ePrediction <- Y(errors, gPrediction, modelname = "oppelt_red_prediction_err", compile = FALSE)
compile(gPrediction, ePrediction, output = "oppelt_red_prediction")
prdPrediction <- gPrediction * xReduced * pSSReduced * pReduced
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
predicted <- function(t) {
  out <- prdPrediction(t, c(bestReduced, gaugeReduced), deriv = FALSE)$TNF
  unname(out[match(t, out[, "time"]), "lognIkb"])
}
# densest at the flanks of the pulse of nuclear IκBα
times <- c(0, 2, seq(5, 55, 5), 57:59, seq(60, 85, 5), 86:89, seq(90, 120, 5), seq(140, 360, 20))

## ---- predictionprofiles
# at time t the logarithm of the prediction is a parameter v, tied to the model by a data point
# of standard deviation 0.01 that enters the data term
predictionProfile <- function(t) {
  objV <- normL2(data, prdPrediction, ePrediction, times = t) +
    constraintL2(muReduced, sigma = 10) +
    datapointL2(name = "lognIkb", time = t, parameter = "v", sigma = 0.01, condition = "TNF",
                attrName = "data")
  start <- c(setNames(as.numeric(bestReduced), names(bestReduced)), v = predicted(t))
  profile(objV, start, "v", fixed = gaugeReduced,
          stepControl = list(stop = "data", stepsize = 0.01, min = 1e-3),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances))
}
predictionProfiles <- parallel::mclapply(times, predictionProfile, mc.cores = getOption("dMod2.cores"))

## ---- predictionband
# the interval of the logarithm, back on the scale of the concentration
band <- do.call(rbind, Map(function(t, pr) {
  ci <- confint(pr, level = 0.95, val.column = "data")
  data.frame(time = t, value = exp(predicted(t)), lower = exp(ci$lower), upper = exp(ci$upper))
}, times, predictionProfiles))

saveRDS(list(band = band, profiles = predictionProfiles), "results/prediction.rds")
