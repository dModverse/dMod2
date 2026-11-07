# The DILI compounds as a grouping of the full model of the publication, on which it selected the
# changes of DCF by an L1 penalty: which rates each compound changes, by spike-and-slab random
# effects. Saves results/compoundsFull.rds.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- fullcompounddata
# HepG2 cells: the first blot scale with and without TNFα and DCF, the third with TNFα and the
# other four compounds; a condition without a compound belongs to no group
cells <- subset(raw, cell == "HepG2" & (scale == 1 | (scale == 3 & treatment == "TNF")) &
                  (scale == 3 | compound %in% c("no", "DCF")))
cells$compound[cells$compound == "no"] <- NA
cells$blot <- paste0("blot", cells$scale)
cells$condition <- with(cells, paste(ifelse(treatment == "TNF", "TNF", "control"),
                                     ifelse(is.na(compound), "none", compound), blot, sep = "_"))
dataFull <- as.datalist(cells[c("name", "time", "value", "sigma", "condition", "treatment",
                                "compound", "blot")],
                        splitBy = "condition", keepCovariates = c("treatment", "compound", "blot"))


## ---- fullcompoundmodel
# the full model with second derivatives forward over reverse; the cells rest without a compound:
# the steady state reads the rate constants as copies `_rest` that the compounds do not change
dm <- c("forward", "reverse", "forward-reverse")
candidatesFull <- setdiff(getParameters(odemodel(reactions, compile = FALSE)),
                          c(reactions$states, names(gauge), "dose_TNF"))
reactionsRestFull <- reactions
reactionsRestFull$rates <- replaceSymbols(candidatesFull, paste0(candidatesFull, "_rest"),
                                          reactionsRestFull$rates)
pSSFull <- Pimpl(reactionsRestFull, modelname = "oppelt_full_ss", compile = FALSE, deriv2 = TRUE,
                 controlsPTC = list(stability = FALSE))
modelFull <- odemodel(reactions, events = events, modelname = "oppelt_full_ode", compile = FALSE,
                      derivMode = dm)
xFull <- Xs(modelFull, options = list(atol = 1e-10, rtol = 1e-8))
gFull <- Y(observables, xFull, modelname = "oppelt_full_obs", compile = FALSE, derivMode = dm)
eFull <- Y(errors, gFull, modelname = "oppelt_full_err", compile = FALSE, derivMode = dm)


## ---- fullcompoundeffects
# every rate constant of the full model but the gauge, with a spike-and-slab random effect over the
# compounds, from the time of TNFα on; the blots of the third scale with scales and offsets of
# their own
blot3Full <- c("s_IkBa", "off_IkBa", "s_pIkBa", "off_pIkBa", "s_pIKK", "off_pIKK", "s_pNFkB", "off_pNFkB")
trafoFull <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(pSSFull, xFull, gFull, eFull), reactions$states)) |>
  insert("x <- exp(x)", x = setdiff(.currentSymbols, "K_R")) |>
  define("TNF <- 0") |>
  branch(covariates(dataFull)) |>
  define("dose_TNF <- 0", conditionMatch = "^control") |>
  define("dose_TNF <- 1", conditionMatch = "^TNF") |>
  insert("x <- y", x = blot3Full, y = paste0(blot3Full, "_blot3"), conditionMatch = "blot3$") |>
  insert(x ~ 1 + (1 | gr(compound, dist = spikeSlab())), x = candidatesFull) |>
  define("x <- exp(y)", x = paste0(candidatesFull, "_rest"), y = candidatesFull)
pFull <- P(trafoFull, method = "explicit", modelname = "oppelt_full_par", compile = FALSE,
           derivMode = dm)

compile(gFull, xFull, pSSFull, pFull, eFull, output = "oppelt_full", cores = getOption("dMod2.cores"))
prdFull <- gFull * xFull * pSSFull * pFull
objFull <- normL2(dataFull, prdFull, eFull)
# the resting state of the cells without a compound stays stable, as in model.R
stabilityFull <- constraintStability(pSSFull, pFull, duration = max(as.data.frame(dataFull)$time))


## ---- fullcompoundfit
# the fixed effects from the full model refitted with the wide prior of profiles.R, those of the
# third blot from the first; the effects of the compounds start at 0, the hyperparameters at a few
# changes
bestFull <- readRDS("results/profiles.rds")$refit$argument
randomFull <- randomEffects(pFull)
populationFull <- setdiff(getParameters(prdFull), c(unlist(randomFull$etas), names(gauge)))
initFull <- bestFull[intersect(populationFull, names(bestFull))]
initFull[paste0(blot3Full, "_blot3")] <- bestFull[blot3Full]
initFull[c("prob_compound", "tau_compound")] <- c(0.2, 0.5)

## ---- fullprior
# the prior of model.R on the fixed effects, the third blot with that of the first
muFull <- c(mu, setNames(mu[blot3Full], paste0(blot3Full, "_blot3")))
sigmaFull <- c(sigma, setNames(sigma[blot3Full], paste0(blot3Full, "_blot3")))
priorFull <- constraintL2(muFull[populationFull], sigma = sigmaFull[populationFull])

## ---- fullcompoundmultistart
# the spread of the starts: per fixed effect half its 68 % profile interval, between 0.1 and 0.5
# on the log scale, 0.5 where the profile is open and for the hyperparameters
ciFull <- confint(readRDS("results/profiles.rds")$profiles, level = 0.68)
widthFull <- setNames(pmin(pmax((ciFull$upper - ciFull$lower) / 2, 0.1), 0.5), ciFull$name)
widthFull <- c(widthFull, setNames(widthFull[blot3Full], paste0(blot3Full, "_blot3")))
sdFull <- setNames(rep(0.5, length(initFull)), names(initFull))
knownFull <- intersect(names(sdFull), names(widthFull)[is.finite(widthFull)])
sdFull[knownFull] <- widthFull[knownFull]
# 200 starts: the one above and 199 moved by sdFull; a start whose fit fails is drawn again; every
# start runs to its end, with up to 300 iterations of its trust region, one core each
set.seed(3)
fitFull <- msem(objFull + priorFull + stabilityFull, randomFull, initFull, fixed = gauge,
                control = list(laplace = list(maxOuter = 300, trustControl = list(tolControl = tolerances))),
                fits = 200, sd = sdFull, keepFits = TRUE, cores = getOption("dMod2.cores"))

## ---- fullcompoundsave
# the best fit with the table and the fits of all starts, saved without the model
savedFull <- fitFull
savedFull[c("design", "prdfn", "errfn", "data")] <- NULL
saveRDS(savedFull, "results/compoundsFull.rds")
