# The expensive runs of this project on the bwForCluster Helix, as SLURM array jobs of many tasks
# with one or two cores each, which start sooner than a task on many cores. The models compile on
# this machine; distributedComputing() copies the sources and builds them once on Helix against
# the cppDE installed there. Needs an account on Helix under the host alias helix of ~/.ssh/config,
# with dMod2 and cppDE installed there in ~/R/lib-devel4, in the versions of the machine that
# submits. Run from the project directory, step by step; a job folder in the working directory
# means the job was submitted already, and the call recovers it. onefit.R, sensitivities.R and
# symmetries.R run on this machine.

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
jobname = "oppelt_fits", seed = 1L, build = "compile", libs = "~/R/lib-devel4",
input = c("obj", "mu", "sigma", "free", "gauge", "tolerances"),
recover = dir.exists("oppelt_fits_folder"))

## ---- helixfitscollect
fitsJob$check()
fitsJob$wait()
fits <- as.parframe(results(fitsJob))
run("fit.R", "fitsave")

## ---- helixprofiles
# The profiles of profiles.R, one parameter per array task on 2 cores, from the fit under the prior
# of the profiles
source("model.R")
run("profiles.R", "refit", "^## ---- profiles$")
limits <- c(lower = -7, upper = 7)
profilesJob <- distributedComputing({
  profile(objProfile, best, task, fixed = gauge, limits = limits,
          stepControl = list(stop = "data"), algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances), cores = cores)
}, tasks = as.list(names(best)), on = slurmCluster("helix", cores = 2, walltime = "24:00:00"),
jobname = "oppelt_profiles", seed = 1L, build = "compile", libs = "~/R/lib-devel4",
input = c("objProfile", "best", "gauge", "limits", "tolerances"),
recover = dir.exists("oppelt_profiles_folder"))

## ---- helixprofilescollect
profilesJob$check()
profilesJob$wait()
profiles <- results(profilesJob)
min(profiles$value) - fitProfile$value   # how far a profile reaches below the fit
run("profiles.R", "profilessave")


## ---- helixoptimizers
# optimizers.R: one fit per task, a start by a method, at most 30 minutes; the exact Hessian on 16 GB
source("model.R")
run("sensitivities.R", "reversechain", "^## ---- gradients")
run("optimizers.R", "optimizers", "^optimizers <- do.call")
exact <- tasks$method == "exact"
optimizersJobs <- lapply(c(FALSE, TRUE), function(heavy) {
  distributedComputing({ fitOne(task$method, task$start) },
    tasks = split(tasks[exact == heavy, ], seq_len(sum(exact == heavy))),
    on = slurmCluster("helix", cores = 1, walltime = "00:30:00", memPerCore = if (heavy) 16 else 2),
    jobname = paste0("oppelt_optimizers", if (heavy) "_exact"), seed = 1L, build = "compile",
    libs = "~/R/lib-devel4",
    input = c("fitOne", "objR", "starts", "free", "gauge", "tolerances", "methods"),
    recover = dir.exists(paste0("oppelt_optimizers", if (heavy) "_exact", "_folder")))
})

## ---- helixoptimizerscollect
for (job in optimizersJobs) job$wait()
# a task without result, a fit the queue ended for time or memory, is a fit without value
failedFit <- function(task) data.frame(method = task$method, start = task$start, value = NA_real_,
                                       converged = FALSE, iterations = NA_integer_,
                                       stopReason = NA_character_, elapsed = NA_real_)
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
# the fits on this machine, one profile per task
source("reduced.R")
run("refit.R", "reducedfit", "^## ---- reducedprofiles$")
reducedJob <- distributedComputing({
  profile(objReducedProfile, bestReducedProfile, task, fixed = gaugeReduced,
          limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances), cores = cores)
}, tasks = as.list(names(bestReducedProfile)),
on = slurmCluster("helix", cores = 2, walltime = "12:00:00"),
jobname = "oppelt_reduced", seed = 1L, build = "compile", libs = "~/R/lib-devel4",
input = c("objReducedProfile", "bestReducedProfile", "gaugeReduced", "tolerances"),
recover = dir.exists("oppelt_reduced_folder"))

## ---- helixreducedcollect
reducedJob$wait()
profilesReduced <- results(reducedJob)
min(profilesReduced$value) - fitReducedProfile$value   # how far a profile reaches below the fit
run("refit.R", "reducedsave")


## ---- helixfundamental
# fundamental.R on the reduced model in the chart of the fundamental domain: the fits on this
# machine, one profile per task
source("reduced.R")
run("fundamental.R", "fundamentaltrafo", "^## ---- fundamentalprofiles$")
fundamentalJob <- distributedComputing({
  profile(objFundamentalProfile, bestFundamentalProfile, task, fixed = gaugeReduced,
          limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances), cores = cores)
}, tasks = as.list(names(bestFundamentalProfile)),
on = slurmCluster("helix", cores = 2, walltime = "12:00:00"),
jobname = "oppelt_fundamental", seed = 1L, build = "compile", libs = "~/R/lib-devel4",
input = c("objFundamentalProfile", "bestFundamentalProfile", "gaugeReduced", "tolerances"),
recover = dir.exists("oppelt_fundamental_folder"))

## ---- helixfundamentalcollect
fundamentalJob$wait()
profilesFundamental <- results(fundamentalJob)
# how far a profile reaches below the fit
min(profilesFundamental$value) - fitFundamentalProfile$value
run("fundamental.R", "fundamentalsave")

## ---- helixprediction
# prediction.R, after refit.R: the profile of the prediction at one time per task
source("reduced.R")
run("prediction.R", "predictionsetup", "^predictionProfiles <- ")
predictionJob <- distributedComputing({ predictionProfile(task) }, tasks = as.list(times),
  on = slurmCluster("helix", cores = 1, walltime = "06:00:00"), jobname = "oppelt_prediction",
  seed = 1L, build = "compile", libs = "~/R/lib-devel4",
  input = c("predictionProfile", "predicted", "data", "prdPrediction", "ePrediction", "muReduced", "sigmaReduced",
            "bestReduced", "gaugeReduced", "tolerances"),
  recover = dir.exists("oppelt_prediction_folder"))

## ---- helixpredictioncollect
predictionJob$wait()
predictionProfiles <- results(predictionJob, combine = FALSE)
run("prediction.R", "predictionband")


## ---- helixdcfspikeslab
# the spike-and-slab fit of dcf.R: the 1000 starts of msem(), one per task, at most 30 minutes each
run("dcf.R", to = "^## ---- dcfspikeslab$")
dcfJob <- distributedComputing({
  set.seed(3)
  msem(objSS, randomSS, c(initDCF, prob_compound = 0.2, tau_compound = 0.5), fixed = gaugeReduced,
       control = list(sparseFit = list(trustControl = list(tolControl = tolerances))),
       fits = 1000, sd = 1, which = task, cores = 1)
}, tasks = as.list(1:1000), on = slurmCluster("helix", cores = 1, walltime = "00:30:00", memPerCore = 4),
jobname = "oppelt_dcf_ss", seed = 1L, build = "compile", libs = "~/R/lib-016d",
input = c("objSS", "randomSS", "initDCF", "gaugeReduced", "tolerances"),
recover = dir.exists("oppelt_dcf_ss_folder"))

## ---- helixdcfspikeslabcollect
# a start at the time limit has no result; get() leaves it out with a warning
dcfJob$wait()
fitDCF <- dcfJob$get()

## ---- helixdcfspikeslabmore
# the second round of dcf.R: 200 starts around the best fit
bestSS <- fitDCF$argument
moreJob <- distributedComputing({
  set.seed(4)
  msem(objSS, randomSS, bestSS, fixed = gaugeReduced,
       control = list(sparseFit = list(trustControl = list(tolControl = tolerances))),
       fits = 200, sd = 1, which = task, cores = 1)
}, tasks = as.list(1:200), on = slurmCluster("helix", cores = 1, walltime = "00:30:00", memPerCore = 4),
jobname = "oppelt_dcf_ss_more", seed = 1L, build = "compile", libs = "~/R/lib-016d",
input = c("objSS", "randomSS", "bestSS", "gaugeReduced", "tolerances"),
recover = dir.exists("oppelt_dcf_ss_more_folder"))

## ---- helixdcfspikeslabmorecollect
moreJob$wait()
moreDCF <- moreJob$get()
starts <- rbind(cbind(fitDCF$multistart, round = 1L), cbind(moreDCF$multistart, round = 2L))
if (moreDCF$value < fitDCF$value) fitDCF <- moreDCF
fitDCF$multistart <- starts
run("dcf.R", "dcfspikeslabsave", "^## ---- dcfL1$")

## ---- helixdcfl1
# the L1 scan of dcf.R with fitTasks() on the cluster: the tasks of each call in array jobs of at
# most 1000 tasks, one fit per task, at most ten minutes; a job goes out when at most 1048 tasks of
# the user are queued, since the queue holds at most 2048 of one user, and comes back when it
# finished. The jobs of a call are numbered, so a second run of the section recovers them.
run("dcf.R", to = "^## ---- dcfspikeslab$")
run("dcf.R", "dcfL1", "^## ---- dcfL1fit$")
queued <- function() as.integer(system("ssh helix 'squeue -u $USER -h -r | wc -l'", intern = TRUE))
calls <- 0L
fitTasks <- function(tasks, fn) {
  calls <<- calls + 1L
  parts <- split(seq_along(tasks), ceiling(seq_along(tasks) / 1000))
  out <- vector("list", length(parts))
  jobs <- list()
  k <- 1L
  while (k <= length(parts) || length(jobs)) {
    while (k <= length(parts)) {
      name <- paste0("oppelt_dcf_l1_", calls, "_", k)
      if (!dir.exists(paste0(name, "_folder")) && queued() > 1048L) break
      jobs[[as.character(k)]] <- distributedComputing({ fn(task) }, tasks = tasks[parts[[k]]],
        on = slurmCluster("helix", cores = 1, walltime = "00:10:00", memPerCore = 2),
        jobname = name, seed = 1L, build = "compile", libs = "~/R/lib-016d",
        input = c("fn", "objDCF", "initDCF", "changesDCF", "changesL1", "offL1", "zeroDCF",
                  "gaugeReduced", "tolerances"), recover = dir.exists(paste0(name, "_folder")))
      k <- k + 1L
    }
    for (j in names(jobs)) {
      if (!isTRUE(suppressWarnings(jobs[[j]]$check()))) next
      out[[as.integer(j)]] <- suppressWarnings(jobs[[j]]$get(combine = FALSE))
      jobs[[j]] <- NULL
    }
    if (length(jobs)) Sys.sleep(60)
  }
  Filter(function(r) is.list(r) && !inherits(r, "try-error"), unlist(out, recursive = FALSE))
}
run("dcf.R", "dcfL1fit")
