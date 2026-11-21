# The DILI compounds as a grouping of the reduced model: which rates each compound changes, by
# spike-and-slab random effects, and the L1 scan of the same changes. Saves results/compounds.rds,
# results/compoundsCheck.rds and results/compoundsL1.rds.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- compounddata
# HepG2 cells: the first blot scale with and without TNFα and DCF, the third with TNFα and the
# other four compounds; a condition without a compound belongs to no group
cells <- subset(raw, cell == "HepG2" & (scale == 1 | (scale == 3 & treatment == "TNF")) &
                  (scale == 3 | compound %in% c("no", "DCF")))
cells$compound[cells$compound == "no"] <- NA
cells$blot <- paste0("blot", cells$scale)
cells$condition <- with(cells, paste(ifelse(treatment == "TNF", "TNF", "control"),
                                     ifelse(is.na(compound), "none", compound), blot, sep = "_"))
dataCompounds <- as.datalist(cells[c("name", "time", "value", "sigma", "condition", "treatment",
                                     "compound", "blot")],
                             splitBy = "condition", keepCovariates = c("treatment", "compound", "blot"))


## ---- compoundmodel
# the reduced model with second derivatives forward over reverse, which the mixed-effects fit uses;
# the resting state of reduced.R holds for every value of the rate constants, so a compound
# changes the dynamics from the time of TNFα on
dm <- c("forward", "reverse", "forward-reverse")
candidates <- setdiff(getParameters(odemodel(reactionsReduced, compile = FALSE)),
                      c(reactionsReduced$states, names(gauge), names(fixedReduced), "dose_TNF"))
modelCompounds <- odemodel(reactionsReduced, events = events, modelname = "oppelt_cmp_ode",
                           compile = FALSE, derivMode = dm)
xCompounds <- Xs(modelCompounds, options = list(atol = 1e-10, rtol = 1e-8))
gCompounds <- Y(observablesReduced, xCompounds, modelname = "oppelt_cmp_obs", compile = FALSE,
                derivMode = dm)
eCompounds <- Y(errors, gCompounds, modelname = "oppelt_cmp_err", compile = FALSE, derivMode = dm)


## ---- compoundeffects
# every rate constant the reduced model estimates, with a spike-and-slab random effect over the
# compounds; the blots of the third scale with scales and offsets of their own
blot3 <- setdiff(c("s_IkBa", "off_IkBa", "s_pIkBa", "off_pIkBa", "s_pIKK", "off_pIKK", "s_pNFkB",
                   "off_pNFkB"), names(fixedReduced))
trafoCompounds <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(xCompounds, gCompounds, eCompounds),
                               reactionsReduced$states)) |>
  define("x <- y", x = names(fixedReduced), y = fixedReduced) |>
  define("x <- y", x = names(resting), y = resting) |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(dataCompounds)) |>
  define("dose_TNF <- 0", conditionMatch = "^control") |>
  define("dose_TNF <- 1", conditionMatch = "^TNF") |>
  insert("x <- y", x = blot3, y = paste0(blot3, "_blot3"), conditionMatch = "blot3$") |>
  insert(x ~ 1 + (1 | gr(compound, dist = spikeSlab())), x = candidates)
pCompounds <- P(trafoCompounds, method = "explicit", modelname = "oppelt_cmp_par", compile = FALSE,
                derivMode = dm)

compile(gCompounds, xCompounds, pCompounds, eCompounds, output = "oppelt_cmp",
        cores = getOption("dMod2.cores"))
prdCompounds <- gCompounds * xCompounds * pCompounds
objCompounds <- normL2(dataCompounds, prdCompounds, eCompounds)


## ---- compoundfit
# the fixed effects from the reduced model at the optimum, those of the third blot from the first;
# the effects of the compounds start at 0, the hyperparameters at a few changes of moderate size
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
random <- randomEffects(pCompounds)
population <- setdiff(getParameters(prdCompounds), c(unlist(random$etas), names(gaugeReduced)))
init <- bestReduced[intersect(population, names(bestReduced))]
init[paste0(blot3, "_blot3")] <- bestReduced[blot3]
init[c("prob_compound", "tau_compound")] <- c(0.2, 0.5)
# the prior of model.R, that of the first blot also for the third
third <- function(v) c(v, setNames(v[blot3], paste0(blot3, "_blot3")))[population]
priorCompounds <- constraintL2(third(muReduced), sigma = third(sigmaReduced))
# the spread of the starts: per fixed effect half the 68 % profile interval of the reduced model,
# between 0.1 and 0.5 on the log scale, 0.5 where the profile is open and for the hyperparameters
ci <- confint(readRDS("results/reduced.rds")$profiles, level = 0.68)
width <- setNames(pmin(pmax((ci$upper - ci$lower) / 2, 0.1), 0.5), ci$name)
width <- c(width, setNames(width[intersect(blot3, names(width))], paste0(intersect(blot3, names(width)), "_blot3")))
sdStart <- setNames(rep(0.5, length(init)), names(init))
known <- intersect(names(sdStart), names(width)[!is.na(width)])
sdStart[known] <- width[known]
# 200 starts: the one above and 199 moved by sdStart, the probability of a change on the logit
# scale and the spread by a factor; a start whose fit fails is drawn again; every start runs to the
# end, with up to 1000 iterations of the trust region, one core each
set.seed(3)
fitCompounds <- msem(objCompounds + priorCompounds, random, init, fixed = gaugeReduced,
                     control = list(sparseFit = list(maxOuter = 1000, trustControl = list(tolControl = tolerances))),
                     fits = 200, sd = sdStart, keepFits = TRUE,
                     cores = getOption("dMod2.cores"))

## ---- compoundsave
# the best fit with the table and the fits of all starts, saved without the model, which an
# evaluation at the estimates restores
saved <- fitCompounds
saved[c("design", "prdfn", "errfn", "data")] <- NULL
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "reduced.R", "compounds.R", "data/oppelt2018.csv", "results/reduced.rds")
attr(saved, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/compounds.rds")

## ---- compoundcheck
# -2 log L of the best fit by importance sampling on the model, 5000 draws per compound, against
# the approximation of the fit at the same estimates
set.seed(4)
checkCompounds <- logLik(fitCompounds, method = "importance", draws = 5000)
attr(checkCompounds, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(checkCompounds, "results/compoundsCheck.rds")

## ---- compoundsL1
# The L1 scan of the publication on the same data and model: every change of a rate constant by a
# compound is a parameter of the objective with the penalty lambda |eta|, one lambda for all
# compounds; changes that reach no data of the compound stay at 0
fitSS <- readRDS("results/compounds.rds")
reach <- fitSS$random$data %||% rep(TRUE, nrow(fitSS$random))
etasL1 <- fitSS$random$parameter[reach]
fixedL1 <- c(gaugeReduced, setNames(numeric(sum(!reach)), fitSS$random$parameter[!reach]))
objL1 <- objCompounds + priorCompounds
# two starts per penalty: the spike-and-slab estimates with the modes of the changes, and the start
# of compoundfit with all changes at 0
effects <- setdiff(names(init), c("prob_compound", "tau_compound"))
startsL1 <- list(spikeSlab = c(fitSS$argument[effects], setNames(fitSS$random$mode[reach], etasL1)),
                 reduced = c(init[effects], setNames(numeric(length(etasL1)), etasL1)))
lambdas <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 60, 100)
tasksL1 <- expand.grid(lambda = lambdas, start = names(startsL1), stringsAsFactors = FALSE)
# the penalized fit of a task, with -2 log L without the penalty
penalizedL1 <- function(task) {
  fit <- trustL1(objL1, startsL1[[task$start]], mu = setNames(numeric(length(etasL1)), etasL1),
                 lambda = task$lambda, rinit = 1, rmax = 10, iterlim = 1000, fixed = fixedL1,
                 tolControl = tolerances)
  list(lambda = task$lambda, start = task$start, argument = fit$argument,
       value = objL1(fit$argument, fixed = fixedL1, deriv = FALSE)$value,
       converged = fit$converged, iterations = fit$iterations)
}
# the changes a penalized fit keeps, fitted again without penalty, the others fixed at 0
refitL1 <- function(task) {
  off <- setdiff(etasL1, task$support)
  fixed <- c(fixedL1, setNames(numeric(length(off)), off))
  fit <- trust(objL1, task$argument[setdiff(names(task$argument), off)], rinit = 1, rmax = 10,
               iterlim = 1000, fixed = fixed, tolControl = tolerances)
  list(support = task$support, value = fit$value, converged = fit$converged,
       iterations = fit$iterations)
}

## ---- compoundsL1fit
# one fit per task, on one core each
penalized <- parallel::mclapply(split(tasksL1, seq_len(nrow(tasksL1))), penalizedL1,
                                mc.cores = getOption("dMod2.cores"))

## ---- compoundsL1refit
# per penalty the better start; each set of changes kept once, refitted from the penalized fit
values <- vapply(penalized, function(r) r$value + r$lambda * sum(abs(r$argument[etasL1])), 0)
bestL1 <- lapply(split(seq_along(penalized), tasksL1$lambda), function(i) penalized[[i[which.min(values[i])]]])
supports <- lapply(bestL1, function(r) etasL1[r$argument[etasL1] != 0])
first <- !duplicated(vapply(supports, paste, "", collapse = " "))
refits <- parallel::mclapply(Map(function(r, s) list(support = s, argument = r$argument),
                                 bestL1[first], supports[first]),
                             refitL1, mc.cores = getOption("dMod2.cores"))

## ---- compoundsL1save
# per penalty the number of changes kept and -2 log L of the refit; the likelihood ratio test of each
# refit against the refit with all changes, the selection that of the largest penalty not worse at
# p < 0.001, the threshold of the publication
key <- function(s) paste(s, collapse = " ")
refitValue <- setNames(vapply(refits, `[[`, 0, "value"), vapply(refits, function(r) key(r$support), ""))
scan <- data.frame(lambda = as.numeric(names(bestL1)), changes = lengths(supports),
                   penalized = vapply(bestL1, `[[`, 0, "value"),
                   converged = vapply(bestL1, `[[`, TRUE, "converged"),
                   refit = unname(refitValue[vapply(supports, key, "")]), row.names = NULL)
full <- scan$refit[scan$lambda == 0]
scan$p <- pchisq(scan$refit - full, df = max(scan$changes) - scan$changes, lower.tail = FALSE)
scan$p[scan$changes == max(scan$changes)] <- 1
chosen <- max(scan$lambda[scan$p >= 0.001])
path <- t(vapply(bestL1, function(r) r$argument[etasL1], numeric(length(etasL1))))
rownames(path) <- names(bestL1)
selectedL1 <- fitSS$random[match(supports[[as.character(chosen)]], fitSS$random$parameter),
                           c("group", "eta", "parameter")]
selectedL1$value <- bestL1[[as.character(chosen)]]$argument[selectedL1$parameter]
rownames(selectedL1) <- NULL
compoundsL1 <- list(scan = scan, path = path, lambda = chosen, selected = selectedL1,
                    penalized = penalized, refits = refits)
sourcesL1 <- c("model.R", "reduced.R", "compounds.R", "data/oppelt2018.csv", "results/reduced.rds",
               "results/compounds.rds")
attr(compoundsL1, "stamp") <- list(md5 = tools::md5sum(sourcesL1), dMod2 = packageVersion("dMod2"))
saveRDS(compoundsL1, "results/compoundsL1.rds")
