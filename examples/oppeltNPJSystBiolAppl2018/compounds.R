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
# the cells rest without a compound: the steady state reads the rate constants as copies `_rest`
# that the compounds do not change
dm <- c("forward", "reverse", "forward-reverse")
candidates <- setdiff(getParameters(odemodel(reactionsReduced, compile = FALSE)),
                      c(reactionsReduced$states, names(gauge), fast, "dose_TNF"))
reactionsRest <- reactionsReduced
reactionsRest$rates <- replaceSymbols(candidates, paste0(candidates, "_rest"), reactionsRest$rates)
pSSCompounds <- Pimpl(reactionsRest, modelname = "oppelt_cmp_ss", compile = FALSE, deriv2 = TRUE,
                      controlsPTC = list(stabilityTol = 1e-3))
modelCompounds <- odemodel(reactionsReduced, events = events, modelname = "oppelt_cmp_ode",
                           compile = FALSE, derivMode = dm)
xCompounds <- Xs(modelCompounds, options = list(atol = 1e-10, rtol = 1e-8))
gCompounds <- Y(observablesReduced, xCompounds, modelname = "oppelt_cmp_obs", compile = FALSE,
                derivMode = dm)
eCompounds <- Y(errors, gCompounds, modelname = "oppelt_cmp_err", compile = FALSE, derivMode = dm)


## ---- compoundeffects
# every rate constant the reduced model estimates, with a spike-and-slab random effect over the
# compounds, from the time of TNFα on; the blots of the third scale with scales and offsets of
# their own
blot3 <- setdiff(c("s_IkBa", "off_IkBa", "s_pIkBa", "off_pIkBa", "s_pIKK", "off_pIKK", "s_pNFkB",
                   "off_pNFkB"), vanishing)
trafoCompounds <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(pSSCompounds, xCompounds, gCompounds, eCompounds),
                               reactionsReduced$states)) |>
  define("x <- 1000", x = fast) |>
  define("x <- 0", x = vanishing) |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(dataCompounds)) |>
  define("dose_TNF <- 0", conditionMatch = "^control") |>
  define("dose_TNF <- 1", conditionMatch = "^TNF") |>
  insert("x <- y", x = blot3, y = paste0(blot3, "_blot3"), conditionMatch = "blot3$") |>
  insert(x ~ 1 + (1 | gr(compound, dist = spikeSlab())), x = candidates) |>
  define("x <- exp(y)", x = paste0(candidates, "_rest"), y = candidates)
pCompounds <- P(trafoCompounds, method = "explicit", modelname = "oppelt_cmp_par", compile = FALSE,
                derivMode = dm)

compile(gCompounds, xCompounds, pSSCompounds, pCompounds, eCompounds, output = "oppelt_cmp",
        cores = getOption("dMod2.cores"))
prdCompounds <- gCompounds * xCompounds * pSSCompounds * pCompounds
objCompounds <- normL2(dataCompounds, prdCompounds, eCompounds)


## ---- compoundfit
# the fixed effects from the reduced model at its optimum, those of the third blot from the first;
# the effects of the compounds start at 0, the hyperparameters at a few changes of moderate size
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
random <- randomEffects(pCompounds)
population <- setdiff(getParameters(prdCompounds), c(unlist(random$etas), names(gaugeReduced)))
init <- bestReduced[intersect(population, names(bestReduced))]
init[paste0(blot3, "_blot3")] <- bestReduced[blot3]
init[c("prob_compound", "tau_compound")] <- c(0.2, 0.5)
priorCompounds <- constraintL2(c(muReduced, setNames(muReduced[blot3], paste0(blot3, "_blot3")))[population],
                               sigma = 10)
# ten starts: the one above and nine with the fixed effects moved by a standard deviation of 0.5 on
# the log scale, the probability of a change on the logit scale and the spread by a factor
set.seed(3)
moved <- setdiff(names(init), c("prob_compound", "tau_compound"))
startsCompounds <- c(list(init), lapply(1:9, function(k) {
  s <- init
  s[moved] <- s[moved] + rnorm(length(moved), sd = 0.5)
  s[["prob_compound"]] <- plogis(qlogis(0.2) + rnorm(1, sd = 0.5))
  s[["tau_compound"]] <- 0.5 * exp(rnorm(1, sd = 0.5))
  s
}))
# the fit from start k, saved without the model, which an evaluation at the estimates restores
fitCompound <- function(k) {
  fit <- tryCatch(em(objCompounds + priorCompounds, random, startsCompounds[[k]], fixed = gaugeReduced,
                     control = list(laplace = list(trustControl = list(tolControl = tolerances))),
                     cores = getOption("dMod2.cores")), error = function(e) conditionMessage(e))
  if (is.character(fit)) return(fit)
  fit[c("design", "prdfn", "errfn", "data")] <- NULL
  fit
}
fitsCompounds <- lapply(seq_along(startsCompounds), fitCompound)

## ---- compoundsave
# the best fit and a table of all starts
ok <- !vapply(fitsCompounds, is.character, TRUE)
starts <- data.frame(start = seq_along(fitsCompounds),
                     value = vapply(fitsCompounds, function(f) if (is.character(f)) NA_real_ else f$value, 0),
                     converged = vapply(fitsCompounds, function(f) !is.character(f) && isTRUE(f$converged), TRUE),
                     error = vapply(fitsCompounds, function(f) if (is.character(f)) f else NA_character_, ""))
saved <- fitsCompounds[[which(ok)[which.min(starts$value[ok])]]]
saved$multistart <- starts
saveRDS(saved, "results/compounds.rds")
