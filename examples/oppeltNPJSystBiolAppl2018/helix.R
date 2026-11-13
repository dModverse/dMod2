# The expensive runs of this project on the bwForCluster Helix, as SLURM array jobs of many tasks
# with one or two cores each, which start sooner than a task on many cores. The models compile on
# this machine; distributedComputing() copies the sources and builds them once on Helix against
# the cppDE installed there. Needs an account on Helix under the host alias helix of ~/.ssh/config,
# with dMod2 and cppDE installed there in the versions of the machine that submits. Run from the
# project directory, step by step; a job folder in the working directory means the job was
# submitted already, and the call recovers it. onefit.R, sensitivities.R and symmetries.R run on
# this machine.

library(dMod2)
# the job script dMod2 ships for Helix, set up once
slurmTemplate("helix", edit = FALSE)

## ---- helixhelpers
# the code of a script from the label `from` (or the top) up to the first line that matches `to`
# (or the end), evaluated in the session
run <- function(file, from = NULL, to = NULL) {
  l <- readLines(file)
  i <- if (is.null(from)) 1L else grep(paste0("^## ---- ", from, "$"), l) + 1L
  j <- c(if (!is.null(to)) grep(to, l), length(l) + 1L)
  j <- j[j > i][1]
  eval(parse(text = l[i:(j - 1L)]), globalenv())
}
# the results of the tasks of a finished job, joined by get() unless `combine = FALSE`: parts of
# one msem() multistart into one fit, data frames by rbind(); a task without a result is an error
results <- function(job, combine = TRUE)
  withCallingHandlers(job$get(combine = combine),
                      warning = function(w) stop(conditionMessage(w), call. = FALSE))


## ---- helixfits
# The multistart of fit.R as 100 array tasks of 2 fits on 2 cores each. distributedComputing()
# copies the model sources and compiles them once on Helix against the cppDE installed there;
# task i starts with the seed 1 + i. The fits come back as one parframe in results/fits.rds.
source("model.R")
free <- setdiff(names(mu), names(gauge))
fitsJob <- distributedComputing({
  mstrust(obj, mu[free], fixed = gauge, sd = sigma[free], fits = 2, iterlim = 1000,
          tolControl = tolerances, cores = cores) |> as.parframe()
}, tasks = 100, on = slurmCluster("helix", cores = 2, walltime = "08:00:00"),
jobname = "oppelt_fits", seed = 1L, build = "compile",
input = c("obj", "mu", "sigma", "free", "gauge", "tolerances"),
recover = dir.exists("oppelt_fits_folder"))

## ---- helixfitscollect
fitsJob$check()
fitsJob$wait()
fits <- as.parframe(results(fitsJob))
run("fit.R", "fitsave")

## ---- helixprofiles
# The profiles of profiles.R, one parameter per array task on 2 cores, from the fit under the prior
# of the profiles; with a profile below that fit, the fit restarts at the lowest point of the
# profiles (`lowest`, from helixprofilescollect) and the job runs again.
run("profiles.R", "refit", "^## ---- profiles$")
if (exists("lowest")) {
  fitProfile <- trust(objProfile, lowest, fixed = gauge, rinit = 0.1, rmax = 10, iterlim = 1000,
                      tolControl = tolerances)
  best <- fitProfile$argument
}
limits <- c(lower = -7, upper = 7)
profilesJob <- distributedComputing({
  profile(objProfile, best, task, fixed = gauge, limits = limits,
          stepControl = list(stop = "data"), algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances), cores = cores)
}, tasks = as.list(names(best)), on = slurmCluster("helix", cores = 2, walltime = "24:00:00"),
jobname = "oppelt_profiles", seed = 1L, build = "compile",
input = c("objProfile", "best", "gauge", "limits", "tolerances"),
recover = dir.exists("oppelt_profiles_folder"))

## ---- helixprofilescollect
profilesJob$check()
profilesJob$wait()
profiles <- results(profilesJob)
min(profiles$value) - fitProfile$value   # below -1e-2: helixprofiles again, from `lowest`
lowest <- unlist(as.data.frame(profiles)[which.min(profiles$value), names(best)])
run("profiles.R", "profilessave")


## ---- helixoptimizers
# optimizers.R: one fit per task, a start by a method; the exact Hessian on 16 GB
source("model.R")
run("sensitivities.R", "reversechain", "^## ---- gradients")
run("optimizers.R", "optimizers", "^optimizers <- do.call")
exact <- tasks$method == "exact"
optimizersJobs <- lapply(c(FALSE, TRUE), function(heavy) {
  distributedComputing({ fitOne(task$method, task$start) },
    tasks = split(tasks[exact == heavy, ], seq_len(sum(exact == heavy))),
    on = slurmCluster("helix", cores = 1, walltime = "06:00:00", memPerCore = if (heavy) 16 else 2),
    jobname = paste0("oppelt_optimizers", if (heavy) "_exact"), seed = 1L, build = "compile",
    input = c("fitOne", "objR", "starts", "free", "gauge", "tolerances", "methods"),
    recover = dir.exists(paste0("oppelt_optimizers", if (heavy) "_exact", "_folder")))
})

## ---- helixoptimizerscollect
for (job in optimizersJobs) job$wait()
# a task without result, a fit the queue ended for memory, is a fit without value
failedFit <- function(task) data.frame(method = task$method, start = task$start, value = NA_real_,
                                       converged = FALSE, iterations = NA_integer_, elapsed = NA_real_)
optimizers <- do.call(rbind, Map(function(job, heavy) {
  out <- suppressWarnings(job$get(combine = FALSE))
  part <- tasks[exact == heavy, ]
  do.call(rbind, lapply(seq_along(out), function(i)
    if (is.data.frame(out[[i]])) out[[i]] else failedFit(part[i, ])))
}, optimizersJobs, c(FALSE, TRUE)))
optimizers <- optimizers[order(optimizers$start, match(optimizers$method, names(methods))), ]
rownames(optimizers) <- NULL
run("optimizers.R", "optimizerssave")


## ---- helixreduced
# refit.R on the reduced model of reduced.R, once the limits are decided against the profiles:
# the fits on this machine, one profile per task; with a profile below the fit, as for the full
# model, again from `lowest`
source("reduced.R")
run("refit.R", "reducedfit", "^## ---- reducedprofiles$")
if (exists("lowestReduced")) {
  fitReducedProfile <- trust(objReducedProfile, lowestReduced, fixed = gaugeReduced, rinit = 0.1,
                             rmax = 10, iterlim = 1000, tolControl = tolerances)
  bestReducedProfile <- fitReducedProfile$argument
}
reducedJob <- distributedComputing({
  profile(objReducedProfile, bestReducedProfile, task, fixed = gaugeReduced,
          limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances), cores = cores)
}, tasks = as.list(names(bestReducedProfile)),
on = slurmCluster("helix", cores = 2, walltime = "12:00:00"),
jobname = "oppelt_reduced", seed = 1L, build = "compile",
input = c("objReducedProfile", "bestReducedProfile", "gaugeReduced", "tolerances"),
recover = dir.exists("oppelt_reduced_folder"))

## ---- helixreducedcollect
reducedJob$wait()
profilesReduced <- results(reducedJob)
min(profilesReduced$value) - fitReducedProfile$value   # below -1e-2: helixreduced again
lowestReduced <- unlist(as.data.frame(profilesReduced)[which.min(profilesReduced$value),
                                                       names(bestReducedProfile)])
run("refit.R", "reducedsave")


## ---- helixprediction
# prediction.R, after refit.R: the profile of the prediction at one time per task
run("prediction.R", "predictionsetup", "^predictionProfiles <- ")
predictionJob <- distributedComputing({ predictionProfile(task) }, tasks = as.list(times),
  on = slurmCluster("helix", cores = 1, walltime = "06:00:00"), jobname = "oppelt_prediction",
  seed = 1L, build = "compile",
  input = c("predictionProfile", "predicted", "data", "prdPrediction", "ePrediction", "muReduced", "sigmaReduced",
            "bestReduced", "gaugeReduced", "tolerances"),
  recover = dir.exists("oppelt_prediction_folder"))

## ---- helixpredictioncollect
predictionJob$wait()
predictionProfiles <- results(predictionJob, combine = FALSE)
run("prediction.R", "predictionband")


## ---- helixcompounds
# compounds.R, after refit.R: the 200 starts of msem(), 2 per task; every task draws all starts
# from the seed of the script and fits the two that `which` names
run("compounds.R", to = "^set.seed\\(3\\)")
compoundsJob <- distributedComputing({
  set.seed(3)
  msem(objCompounds + priorCompounds, random, init, fixed = gaugeReduced,
       control = list(laplace = list(maxOuter = 1000, trustControl = list(tolControl = tolerances))),
       fits = 200, sd = sdStart, keepFits = TRUE, which = task, cores = cores)
}, tasks = split(1:200, rep(1:100, each = 2)),
on = slurmCluster("helix", cores = 2, walltime = "24:00:00"), jobname = "oppelt_compounds",
seed = 1L, build = "compile",
input = c("objCompounds", "priorCompounds", "random", "init", "gaugeReduced", "sdStart",
          "tolerances"), recover = dir.exists("oppelt_compounds_folder"))

## ---- helixcompoundscollect
compoundsJob$wait()
fitCompounds <- results(compoundsJob)
run("compounds.R", "compoundsave")
