# The DILI compounds as a grouping of the reduced model: which rates each compound changes, by
# spike-and-slab random effects. Saves results/compounds.rds.

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
