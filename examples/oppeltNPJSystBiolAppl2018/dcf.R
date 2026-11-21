# The problem of the publication: TNFα-stimulated HepG2 cells of the first blot scale with and
# without DCF, every rate constant of the reduced model changed by DCF. The changes selected by
# spike-and-slab random effects and by an L1 scan as in the publication. Saves results/dcf.rds and
# results/dcfL1.rds.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- dcfdata
# the TNFα-stimulated cells of the first blot scale, without compound and with DCF
cellsDCF <- subset(raw, cell == "HepG2" & scale == 1 & treatment == "TNF" & compound %in% c("no", "DCF"))
cellsDCF$compound[cellsDCF$compound == "no"] <- NA
cellsDCF$condition <- ifelse(is.na(cellsDCF$compound), "TNF_none", "TNF_DCF")
dataDCF <- as.datalist(cellsDCF[c("name", "time", "value", "sigma", "condition", "treatment", "compound")],
                       splitBy = "condition", keepCovariates = c("treatment", "compound"))

## ---- dcfmodel
# every rate constant the publication does not fix has a change by DCF, 24 as in the publication,
# a spike-and-slab random effect over the compound; for the spike and slab without the four rate
# constants of the gauge, whose change equals that of a partner; second derivatives forward over
# reverse for em()
ratesDCF <- setdiff(getParameters(odemodel(reactionsReduced, compile = FALSE)),
                    c(reactionsReduced$states, names(fixedReduced), "Vnuc", "TNF", "dose_TNF"))
dm <- c("forward", "reverse", "forward-reverse")
modelDCF <- odemodel(reactionsReduced, events = events, modelname = "oppelt_dcf_ode", compile = FALSE,
                     derivMode = dm)
# solves without sensitivities at atol 1e-12: at atol 1e-10 they can leave the positive orthant, where
# the complex formation is unstable; those with sensitivities hold at 1e-10
xDCF <- Xs(modelDCF, options = list(atol = 1e-12, rtol = 1e-8), optionsSens = list(atol = 1e-10))
gDCF <- Y(observablesReduced, xDCF, modelname = "oppelt_dcf_obs", compile = FALSE, derivMode = dm)
eDCF <- Y(errors, gDCF, modelname = "oppelt_dcf_err", compile = FALSE, derivMode = dm)
trafoBase <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(xDCF, gDCF, eDCF), reactionsReduced$states)) |>
  define("x <- y", x = names(fixedReduced), y = fixedReduced) |>
  define("x <- y", x = names(resting), y = resting) |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(dataDCF)) |>
  define("dose_TNF <- 1")
ratesSS <- setdiff(ratesDCF, names(gauge))
pDCF <- P(insert(trafoBase, x ~ 1 + (1 | gr(compound, dist = spikeSlab())), x = ratesDCF),
          method = "explicit", modelname = "oppelt_dcf_par", compile = FALSE, derivMode = dm)
pSS <- P(insert(trafoBase, x ~ 1 + (1 | gr(compound, dist = spikeSlab())), x = ratesSS),
         method = "explicit", modelname = "oppelt_dcfss_par", compile = FALSE, derivMode = dm)
compile(gDCF, xDCF, pDCF, pSS, eDCF, output = "oppelt_dcf", cores = getOption("dMod2.cores"))
prdDCF <- gDCF * xDCF * pDCF

## ---- dcfobjective
# the prior of reduced.R on the rate constants and observation parameters; the changes start at 0
randomDCF <- randomEffects(pDCF)
changesDCF <- unique(na.omit(unlist(randomDCF$etas)))
populationDCF <- setdiff(getParameters(prdDCF), c(changesDCF, names(gaugeReduced)))
priorDCF <- constraintL2(muReduced[populationDCF], sigma = sigmaReduced[populationDCF])
objDCF <- normL2(dataDCF, prdDCF, eDCF) + priorDCF
randomSS <- randomEffects(pSS)
objSS <- normL2(dataDCF, gDCF * xDCF * pSS, eDCF) + priorDCF
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
initDCF <- bestReduced[populationDCF]

## ---- dcfspikeslab
# the changes as spike-and-slab random effects: 1000 starts of msem(), the population parameters
# drawn around the fit of the reduced model with sd 1 on the log scale, the probability of a change
# starting at 0.2 and the spread at 0.5
initSS <- c(initDCF, prob_compound = 0.2, tau_compound = 0.5)
set.seed(3)
fitDCF <- msem(objSS, randomSS, initSS, fixed = gaugeReduced,
               control = list(sparseFit = list(trustControl = list(tolControl = tolerances))),
               fits = 1000, sd = 1, cores = getOption("dMod2.cores"))

## ---- dcfspikeslabmore
# 200 more starts around the best fit, sd 1 as above: the waterfall of the optimum
set.seed(4)
moreDCF <- msem(objSS, randomSS, fitDCF$argument, fixed = gaugeReduced,
                control = list(sparseFit = list(trustControl = list(tolControl = tolerances))),
                fits = 200, sd = 1, cores = getOption("dMod2.cores"))
starts <- rbind(cbind(fitDCF$multistart, round = 1L), cbind(moreDCF$multistart, round = 2L))
if (moreDCF$value < fitDCF$value) fitDCF <- moreDCF
fitDCF$multistart <- starts

## ---- dcfspikeslabsave
# saved without the model, which an evaluation at the estimates restores
saved <- fitDCF
saved[c("design", "prdfn", "errfn", "data")] <- NULL
sources <- c("model.R", "reduced.R", "dcf.R", "data/oppelt2018.csv", "results/reduced.rds")
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/dcf.rds")

## ---- dcfL1
# The L1 scan over the 20 candidates of the spike and slab, the changes of the four partners of the
# gauge at 0: the penalty lambda |eta| on every candidate, lambda = 0, 31 values from 1e-2 to 1e4
# evenly on the log scale and four more between 10^0.8 and 10, per lambda 1000 starts, the
# population parameters drawn around the fit of the reduced model with sd 1, the changes at 0, the
# first start the fit itself
zeroDCF <- setNames(numeric(length(changesDCF)), changesDCF)
changesL1 <- paste0("eta_", ratesSS, "_DCF")
offL1 <- setdiff(changesDCF, changesL1)
stopifnot(all(changesL1 %in% changesDCF))
lambdasDCF <- sort(c(0, 10^seq(-2, 4, by = 0.2), 10^seq(0.84, 0.96, by = 0.04)))
startsL1 <- msParframe(c(initDCF, zeroDCF), n = 1000, seed = 12,
                       sd = c(setNames(rep(1, length(initDCF)), names(initDCF)), zeroDCF))
startL1 <- function(i) unlist(startsL1[i, c(names(initDCF), changesL1)])
# a start vector with the population parameters and candidates of `argument`, the others at the start
fullStart <- function(argument) {
  start <- c(initDCF, zeroDCF[changesL1])
  shared <- intersect(names(argument), names(start))
  start[shared] <- argument[shared]
  start
}
# the penalized fit of a task, with -2 log L without the penalty
penalizedDCF <- function(task) {
  fit <- try(trustL1(objDCF, task$start, mu = zeroDCF[changesL1], lambda = task$lambda, rinit = 1,
                     rmax = 10, iterlim = 1000, fixed = c(gaugeReduced, zeroDCF[offL1]),
                     tolControl = tolerances), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  list(lambda = task$lambda, argument = fit$argument,
       value = objDCF(c(fit$argument, zeroDCF[offL1]), fixed = gaugeReduced, deriv = FALSE)$value,
       converged = fit$converged)
}
# the candidates of a set fitted without penalty, the other changes at 0
refitDCF <- function(task) {
  off <- setdiff(changesDCF, task$set)
  start <- task$start[setdiff(names(task$start), off)]
  fit <- try(trust(objDCF, start, rinit = 1, rmax = 10, iterlim = 1000,
                   fixed = c(gaugeReduced, zeroDCF[off]), tolControl = tolerances), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  list(set = task$set, argument = fit$argument, value = fit$value, converged = fit$converged)
}
# the penalized objective of every fit in penalizedFits at `x`
penalizedAt <- function(x) vapply(penalizedFits, function(r) r$value + x * sum(abs(r$argument[changesL1])), 0)
# the fits of a list of tasks, one per task; helix.R runs them on the cluster instead
fitTasks <- function(tasks, fn)
  Filter(Negate(is.null), parallel::mclapply(tasks, fn, mc.cores = getOption("dMod2.cores"),
                                             mc.preschedule = FALSE))

## ---- dcfL1fit
tasksL1 <- expand.grid(start = seq_len(nrow(startsL1)), lambda = lambdasDCF)
penalizedFits <- fitTasks(lapply(seq_len(nrow(tasksL1)), function(i)
  list(start = startL1(tasksL1$start[i]), lambda = tasksL1$lambda[i])), penalizedDCF)

## ---- dcfL1continue
# every lambda also from the five best fits of the other lambdas, evaluated at it, and from the
# five best of the lambda itself, until no lambda improves by more than 0.01
bestAt <- function() vapply(lambdasDCF, function(x) {
  own <- vapply(penalizedFits, `[[`, 0, "lambda") == x
  min(penalizedAt(x)[own])
}, 0)
repeat {
  before <- bestAt()
  lambdaOf <- vapply(penalizedFits, `[[`, 0, "lambda")
  tasks <- unlist(lapply(lambdasDCF, function(x) {
    value <- penalizedAt(x)
    other <- which(lambdaOf != x); own <- which(lambdaOf == x)
    pick <- c(other[order(value[other])][1:5], own[order(value[own])][1:5])
    lapply(pick, function(i) list(start = penalizedFits[[i]]$argument[c(names(initDCF), changesL1)], lambda = x))
  }), recursive = FALSE)
  penalizedFits <- c(penalizedFits, fitTasks(tasks, penalizedDCF))
  if (all(before - bestAt() <= 0.01)) break
}

## ---- dcfL1sets
# per lambda the best penalized fit and the candidates it keeps; the sets refitted without penalty:
# all candidates, those of every lambda, the seven of the publication (Table S7 of the supplement,
# estimates and 95 % intervals) and those the spike and slab selects
lambdaOf <- vapply(penalizedFits, `[[`, 0, "lambda")
bestDCF <- lapply(lambdasDCF, function(x) {
  i <- which(lambdaOf == x)
  penalizedFits[[i[which.min(penalizedAt(x)[i])]]]
})
names(bestDCF) <- lambdasDCF
keptDCF <- lapply(bestDCF, function(r) changesL1[abs(r$argument[changesL1]) > 1e-8])
publicationDCF <- data.frame(
  rate = c("act_Ikk_by_TNF", "trigger_iIkk", "act_pIkk", "act_Ikb_by_Ikk", "prod_mIkb_by_nNfk",
           "degrad_mIkb", "shuttle_RnaA20"),
  value = c(-0.72, 1.56, 0.91, -0.94, -0.83, -1.77, -0.96),
  lower = c(-0.98, 1.10, 0.57, -1.66, -0.94, -2.06, -1.68),
  upper = c(-0.46, 2.05, 1.27, -0.45, -0.73, -1.50, -0.52))
fitSS <- readRDS("results/dcf.rds")
key <- function(set) paste(sort(set), collapse = " ")
setsDCF <- c(list(all = changesL1,
                  publication = paste0("eta_", publicationDCF$rate, "_DCF"),
                  spikeSlab = paste0("eta_", fitSS$random$eta[fitSS$random$p0 < 0.5], "_DCF")),
             keptDCF[lengths(keptDCF) > 0])
setsDCF <- setsDCF[!duplicated(vapply(setsDCF, key, ""))]

## ---- dcfL1refit
# per set the best penalized fit of every lambda and 300 starts around the fit of the reduced model
# (sd 1, the changes at 0); where fewer than five fits reach the best value within 0.1, 300 more,
# then 1000 more and the best refit of every set
drawn <- function(n, seed) {
  d <- msParframe(c(initDCF, zeroDCF[changesL1]), n = n, seed = seed,
                  sd = c(setNames(rep(1, length(initDCF)), names(initDCF)), zeroDCF[changesL1]))
  lapply(seq_len(nrow(d)), function(i) unlist(d[i, c(names(initDCF), changesL1)]))
}
refitsOf <- function(set) Filter(function(r) identical(key(r$set), key(set)), refitFits)
bestRefit <- function(set) { f <- refitsOf(set); f[[which.min(vapply(f, `[[`, 0, "value"))]] }
reached <- function(set) { v <- vapply(refitsOf(set), `[[`, 0, "value"); sum(v < min(v) + 0.1) }
refitFits <- list()
open <- seq_along(setsDCF)
for (round in 1:3) {
  starts <- switch(round,
    c(lapply(bestDCF, function(r) fullStart(r$argument)), drawn(300, 21)),
    drawn(300, 22),
    c(lapply(setsDCF, function(set) fullStart(bestRefit(set)$argument)), drawn(1000, 30)))
  tasks <- unlist(lapply(open, function(k) lapply(starts, function(s) list(set = setsDCF[[k]], start = s))),
                  recursive = FALSE)
  refitFits <- c(refitFits, fitTasks(tasks, refitDCF))
  open <- which(vapply(setsDCF, reached, 0L) < 5L)
  if (!length(open)) break
}

## ---- dcfL1select
# the likelihood ratio test of the set of each lambda against all candidates, the selection that
# of the largest lambda not worse at p < 0.001, the threshold of the publication
valueOf <- function(set) bestRefit(set)$value
scanDCF <- data.frame(lambda = lambdasDCF, changes = lengths(keptDCF),
                      data = vapply(bestDCF, `[[`, 0, "value"),
                      penalty = lambdasDCF * vapply(bestDCF, function(r) sum(abs(r$argument[changesL1])), 0),
                      refit = vapply(keptDCF, function(set) if (length(set)) valueOf(set) else NA_real_, 0),
                      row.names = NULL)
scanDCF$p <- pchisq(scanDCF$refit - valueOf(changesL1), df = length(changesL1) - scanDCF$changes,
                    lower.tail = FALSE)
scanDCF$p[scanDCF$changes == length(changesL1)] <- 1
lambdaDCF <- max(scanDCF$lambda[!is.na(scanDCF$p) & scanDCF$p >= 0.001])
selectedL1 <- keptDCF[[as.character(lambdaDCF)]]
pathDCF <- t(vapply(bestDCF, function(r) r$argument[changesL1], zeroDCF[changesL1]))

## ---- dcfL1save
# the penalized fits and refits as values for the waterfalls, the best ones in full
dcfL1 <- list(scan = scanDCF, path = pathDCF, lambda = lambdaDCF, selected = selectedL1,
              penalizedSelected = bestDCF[[as.character(lambdaDCF)]]$argument,
              refitSelected = bestRefit(selectedL1), publication = publicationDCF,
              refitPublication = bestRefit(setsDCF$publication),
              refitSpikeSlab = bestRefit(setsDCF$spikeSlab), all = bestRefit(changesL1),
              sets = setsDCF, refits = lapply(setsDCF, bestRefit),
              penalized = data.frame(lambda = lambdaOf, value = vapply(penalizedFits, `[[`, 0, "value"),
                                     L1 = vapply(penalizedFits, function(r) sum(abs(r$argument[changesL1])), 0),
                                     converged = vapply(penalizedFits, `[[`, TRUE, "converged")),
              refitValues = data.frame(set = vapply(refitFits, function(r) key(r$set), ""),
                                       value = vapply(refitFits, `[[`, 0, "value")))
sources <- c("model.R", "reduced.R", "dcf.R", "data/oppelt2018.csv", "results/reduced.rds",
             "results/dcf.rds")
attr(dcfL1, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(dcfL1, "results/dcfL1.rds")
