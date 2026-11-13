# Prediction band of nuclear IκBα after TNFα, which no experiment measured, from prediction
# profiles of the reduced model. Saves results/prediction.rds, which plots.R reads.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- predictionsetup
# the readouts of the reduced model and nuclear IκBα, which is 0 at rest
gPrediction <- Y(c(observablesReduced, nIkb = "nIkb"), xReduced,
                 modelname = "oppelt_red_prediction", compile = FALSE)
ePrediction <- Y(errors, gPrediction, modelname = "oppelt_red_prediction_err", compile = FALSE)
compile(gPrediction, ePrediction, output = "oppelt_red_prediction")
prdPrediction <- gPrediction * xReduced * pReduced
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
predicted <- function(t) {
  out <- prdPrediction(t, c(bestReduced, gaugeReduced), deriv = FALSE)$TNF
  unname(out[match(t, out[, "time"]), "nIkb"])
}
# densest at the flanks of the pulse of nuclear IκBα
times <- c(2, seq(5, 55, 5), 57:59, seq(60, 85, 5), 86:89, seq(90, 120, 5), seq(140, 360, 20))

## ---- predictionprofiles
# at time t the prediction is a parameter v, tied to the model by a data point of standard
# deviation 1 % of the prediction that enters the data term; steps relative to the prediction
predictionProfile <- function(t) {
  v <- predicted(t)
  objV <- normL2(data, prdPrediction, ePrediction, times = t) +
    constraintL2(muReduced, sigma = sigmaReduced) +
    datapointL2(name = "nIkb", time = t, parameter = "v", sigma = 0.01 * v, condition = "TNF",
                attrName = "data")
  start <- c(setNames(as.numeric(bestReduced), names(bestReduced)), v = v)
  profile(objV, start, "v", fixed = gaugeReduced,
          stepControl = list(stop = "data", stepsize = 0.01 * v, min = 1e-3 * v),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances))
}
predictionProfiles <- parallel::mclapply(times, predictionProfile, mc.cores = getOption("dMod2.cores"),
                                         mc.preschedule = FALSE)

## ---- predictionband
band <- do.call(rbind, Map(function(t, pr) {
  ci <- confint(pr, level = 0.95, val.column = "data")
  data.frame(time = t, value = predicted(t), lower = ci$lower, upper = ci$upper)
}, times, predictionProfiles))

## ---- predictionsave
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "reduced.R", "prediction.R", "data/oppelt2018.csv", "results/reduced.rds")
saved <- list(band = band, profiles = predictionProfiles)
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/prediction.rds")
