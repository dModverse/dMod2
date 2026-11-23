#
# Three population models with random effects between subjects, fitted by FOCEI and by SAEM, each
# as a multistart of 100 fits, with the diagnostics of `em()` fits. Compartments and volumes are
# those of the reaction network; the rates are fluxes in amounts per time in a frame of volume 1.

library(dMod2)
dir.create("build", showWarnings = FALSE)
# generated sources and shared objects in build/, reused when nothing changed
options(dMod2.outdir = "build")
cores <- max(1L, min(20L, parallel::detectCores() - 2L))
nstart <- 100L
# second derivatives forward over reverse for FOCEI, forward over forward for the reported Hessian
dm <- c("forward", "reverse", "forward-forward", "forward-reverse")

# fits are saved in results/cache.rds under a key and read from there when present; delete the
# file to fit everything again
cached <- function(key, expr) {
  file <- "results/cache.rds"
  dir.create("results", showWarnings = FALSE)
  store <- if (file.exists(file)) readRDS(file) else list()
  if (is.null(store[[key]])) {
    store[[key]] <- expr
    saveRDS(store, file)
  }
  store[[key]]
}
# a fit of em() cached without the model, which an evaluation at the fit estimates restores
cachedEm <- function(key, obj, random, start, control = list(), ...) {
  model <- c("design", "prdfn", "errfn", "data")
  f <- cached(key, { f <- em(obj, random, start, control = control, ...); f[model] <- NULL; f })
  f[model] <- em(obj, random, f, control = modifyList(control, list(trustControl = list(iterlim = 0),
                                                                    boundary = FALSE)),
                 verbose = FALSE)[model]
  f
}

# Multistarts ----
#
# The starts: the initial values on the log scale and `nstart - 1` draws around them, N(0, 0.5)
# per parameter. Every fit runs on one core, `nstart` fits on `cores` cores at once.
drawStarts <- function(init, n = nstart) {
  set.seed(1)
  rbind(init, t(vapply(seq_len(n - 1L), function(i) init + rnorm(length(init), 0, 0.5), init)))
}

# dMod2 from every start; the fits as a table of -2 log L, convergence, seconds and estimates
runDMod2 <- function(obj, random, starts, method) {
  out <- parallel::mclapply(seq_len(nrow(starts)), function(i) {
    options(dMod2.cores = 1L)
    set.seed(i)
    t <- system.time(f <- tryCatch(em(obj, random, starts[i, ], method = method, verbose = FALSE),
                                   error = function(e) NULL))[["elapsed"]]
    if (is.null(f)) return(list(value = NA_real_, converged = FALSE, seconds = t, argument = NULL))
    list(value = f$value, converged = isTRUE(f$converged), seconds = t, argument = f$argument)
  }, mc.cores = cores, mc.preschedule = FALSE)
  data.frame(package = "dMod2", start = seq_along(out), value = vapply(out, `[[`, 0, "value"),
             converged = vapply(out, `[[`, TRUE, "converged"), seconds = vapply(out, `[[`, 0, "seconds"))
}

# per method the best value, the starts within 0.1 of it, the failed and unconverged fits, the
# median and total seconds per fit; `runs` a named list of the tables of runDMod2()
summarizeStarts <- function(runs) {
  do.call(rbind, lapply(names(runs), function(m) {
    r <- runs[[m]]
    best <- min(r$value, na.rm = TRUE)
    data.frame(method = m, best = best, hits = sum(r$value - best < 0.1, na.rm = TRUE),
               failed = sum(is.na(r$value)), unconverged = sum(!r$converged & !is.na(r$value)),
               medianSeconds = stats::median(r$seconds), totalSeconds = sum(r$seconds))
  }))
}

# the waterfall of every method, a panel each, and the seconds per fit
plotStarts <- function(runs, title) {
  w <- plotValues(runs, pseudoLog = TRUE) + ggplot2::labs(title = title)
  seconds <- do.call(rbind, lapply(names(runs), function(m) data.frame(method = m, seconds = runs[[m]]$seconds)))
  s <- ggplot2::ggplot(seconds, ggplot2::aes(method, seconds)) + ggplot2::geom_boxplot() +
    ggplot2::scale_y_log10(labels = scales::label_log()) +
    ggplot2::labs(x = NULL, y = "seconds per fit, one core") + theme_dMod2()
  cowplot::plot_grid(w, s, nrow = 1, rel_widths = c(2, 1))
}

# -2 log L of FOCEI at given estimates, the definition under which the SAEM fits are compared
focei <- function(obj, random, pars)
  em(obj, random, pars, control = list(trustControl = list(iterlim = 0)), verbose = FALSE)$value

# the multistarts of FOCEI and SAEM of one model, and the best SAEM fit with -2 log L of FOCEI at
# the SAEM estimates; `saem` the model SAEM fits
compareMethods <- function(runs, label, saem) {
  runs <- list(FOCEI = runs$focei, SAEM = runs$saem)
  print(summarizeStarts(runs))
  print(plotStarts(runs, label))
  sfit <- cachedEm(paste(label, "dMod2 saem best"), saem$obj, saem$random, saem$starts[which.min(runs$SAEM$value), ],
                   method = "saem", verbose = FALSE)
  print(round(sfit$argument, 4))
  print(c(SAEM = focei(saem$obj, saem$random, sfit$argument)))
  invisible(NULL)
}

# Theophylline ----
#
# Twelve subjects, one oral dose each, a one-compartment model with first-order absorption from a
# depot of volume 1 into the central compartment of volume V; absorption rate, volume and
# clearance vary between subjects.
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
theoData <- as.datalist(data.frame(name = "cp", time = Theoph$Time, value = Theoph$conc, sigma = NA_real_,
                                   condition = Theoph$Subject))
theoNet <- eqnlist() |>
  assignCompartment(Ag = "depot") |>
  assignCompartment(Cc = "central", volume = "V") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("Ag", "Cc", "Ka * Ag", "absorption") |>
  addReaction("Cc", "", "Cl * Cc", "elimination", rateCompartment = "unit")
theoNet
theoX <- Xs(odemodel(theoNet, modelname = "pm_theoph", compile = FALSE, derivMode = dm))
theoG <- Y(c(cp = "Cc"), theoX, modelname = "pm_theoph_g", compile = FALSE, derivMode = dm)
theoE <- Y(c(cp = "sigma"), theoG, attachInput = FALSE, modelname = "pm_theoph_e", compile = FALSE,
           derivMode = dm)
subjects <- sort(unique(Theoph$Subject))
theoTable <- data.frame(ID = subjects, dose = tapply(Theoph$Dose * Theoph$Wt, Theoph$Subject, `[`, 1)[subjects],
                        weight = tapply(Theoph$Wt, Theoph$Subject, `[`, 1)[subjects], row.names = subjects)
theoTrafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma")), table = theoTable)
theoTrafo <- insert(theoTrafo, "Ag <- dose", dose = dose)
theoTrafo <- insert(theoTrafo, "Cc <- 0")
theoTrafo <- insert(theoTrafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma"))
theoTrafo <- insert(theoTrafo, Ka + V + Cl ~ 1 + (1 | ID))
theoP <- P(theoTrafo, method = "explicit", modelname = "pm_theoph_p", compile = FALSE, derivMode = dm)
compileSeconds <- c(theophylline = system.time(compile(theoX, theoG, theoE, theoP,
                                                       output = "pm_theoph_all"))[["elapsed"]])
theoObj <- normL2(theoData, theoG * theoX * theoP, errmodel = theoE)
theoRandom <- randomEffects(theoP)
theoStarts <- drawStarts(c(Ka = 0, V = 3, Cl = 0.7, sigma = log(0.5), omega_Ka_Ka = log(0.3),
                           omega_V_V = log(0.3), omega_Cl_Cl = log(0.3)))

theoRuns <- cached("Theophylline dMod2", list(focei = runDMod2(theoObj, theoRandom, theoStarts, "focei"),
                                               saem = runDMod2(theoObj, theoRandom, theoStarts, "saem")))
theoFit <- cachedEm("Theophylline dMod2 focei best", theoObj, theoRandom,
                    theoStarts[which.min(theoRuns$focei$value), ], verbose = FALSE)
summary(theoFit)

# Diagnostics ----
plot(theoFit)
plotResiduals(theoFit, nsim = 500, seed = 1)
plotCombined(theoFit)
plotRandomEffects(theoFit)
plotRandomEffects(theoFit, covariates = theoTable["weight"])
plotVPC(theoFit, nsim = 500, bins = 8, seed = 1)

# FOCEI and SAEM ----
compareMethods(theoRuns, "Theophylline", list(obj = theoObj, random = theoRandom, starts = theoStarts))


# Warfarin ----
#
# 32 subjects, one oral dose; the plasma concentration `cp` and the prothrombin complex activity
# `pca`, whose synthesis the concentration inhibits (a turnover model with maximal inhibition 1).
# Volume, clearance, baseline and turnover rate of the activity vary between subjects.
data(warfarin)
warfData <- as.datalist(warfarin, splitBy = "condition", keepCovariates = c("dose", "wt", "age", "sex"))
warfNet <- eqnlist() |>
  assignCompartment(Ag = "depot") |>
  assignCompartment(Cc = "central", volume = "V") |>
  assignCompartment(PCA = "effect") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("Ag", "Cc", "Ka * Ag", "absorption") |>
  addReaction("Cc", "", "Cl * Cc", "elimination", rateCompartment = "unit") |>
  addReaction("", "PCA", "kout * pca0 * IC50 / (IC50 + Cc)", "synthesis") |>
  addReaction("PCA", "", "kout * PCA", "degradation")
warfNet
warfX <- Xs(odemodel(warfNet, modelname = "pm_warf", compile = FALSE, derivMode = dm))
warfG <- Y(c(cp = "Cc", pca = "PCA"), warfX, modelname = "pm_warf_g", compile = FALSE,
           derivMode = dm)
# the error of cp combined from an absolute and a relative part, positive also where the solver
# leaves cp at or below 0
warfE <- Y(c(cp = "sqrt(sigma_cp_abs^2 + (sigma_cp_rel * cp)^2)", pca = "sigma_pca"), warfG,
           attachInput = FALSE, modelname = "pm_warf_e", compile = FALSE, derivMode = dm)
warfTable <- attr(warfData, "condition.grid")
warfTable$ID <- rownames(warfTable)
warfTrafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "IC50", "kout", "pca0", "Ag", "Cc", "PCA",
                                                  "sigma_cp_abs", "sigma_cp_rel", "sigma_pca")), table = warfTable)
warfTrafo <- insert(warfTrafo, "Ag <- d", d = dose)
warfTrafo <- insert(warfTrafo, "Cc <- 0")
warfTrafo <- insert(warfTrafo, "PCA <- pca0")
warfTrafo <- insert(warfTrafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "IC50", "kout", "pca0", "sigma_cp_abs",
                                                    "sigma_cp_rel", "sigma_pca"))
warfTrafo <- insert(warfTrafo, V + Cl + pca0 + kout ~ 1 + (1 | ID))
warfP <- P(warfTrafo, method = "explicit", modelname = "pm_warf_p", compile = FALSE, derivMode = dm)
compileSeconds["warfarin"] <- system.time(compile(warfX, warfG, warfE, warfP, output = "pm_warf_all"))[["elapsed"]]
warfObj <- normL2(warfData, warfG * warfX * warfP, errmodel = warfE)
warfRandom <- randomEffects(warfP)
warfStarts <- drawStarts(c(Ka = 0, V = 2, Cl = -2, IC50 = 0, kout = -3, pca0 = 4.6, sigma_cp_abs = log(0.1),
                           sigma_cp_rel = log(0.1),
                           sigma_pca = log(5), omega_V_V = log(0.3), omega_Cl_Cl = log(0.3),
                           omega_pca0_pca0 = log(0.3), omega_kout_kout = log(0.3)))
warfRuns <- cached("Warfarin dMod2", list(focei = runDMod2(warfObj, warfRandom, warfStarts, "focei"),
                                           saem = runDMod2(warfObj, warfRandom, warfStarts, "saem")))
warfFit <- cachedEm("Warfarin dMod2 focei best", warfObj, warfRandom,
                    warfStarts[which.min(warfRuns$focei$value), ], verbose = FALSE)
summary(warfFit)

# Diagnostics ----
plot(warfFit)
plotResiduals(warfFit, nsim = 500, seed = 1)
plotCombined(warfFit, conditionsPerPage = 8)[[1]]
plotRandomEffects(warfFit)
plotRandomEffects(warfFit, covariates = c("wt", "age", "sex"))
plotVPC(warfFit, nsim = 500, bins = 8, seed = 1)

# FOCEI and SAEM ----
#
# The SAEM fits run all 500 iterations and end unconverged: the expected increase of -2 log L from
# the Monte Carlo error of the averaged estimate stays above `mcTol` (0.05).
compareMethods(warfRuns, "Warfarin", list(obj = warfObj, random = warfRandom, starts = warfStarts))


# Mavoglurant ----
#
# 120 subjects, an infusion of 10 min (20 min in one occasion) as a time switch, 78 subjects on two
# occasions. Two compartments; clearances scale with weight to the power 0.75, volumes with weight.
# Clearance, the central and peripheral volume and the intercompartmental clearance vary between
# subjects, the clearance also between the occasions of a subject, `(1 | ID:OCC)`. SAEM fits the
# model without the variability between occasions.
data(mavoglurant)
mavoData <- as.datalist(mavoglurant, splitBy = "condition",
                        keepCovariates = c("ID", "OCC", "dose", "duration", "AGE", "SEX", "WT", "HT"))
mavoNet <- eqnlist() |>
  assignCompartment(C1 = "central", volume = "V1") |>
  assignCompartment(C2 = "peripheral", volume = "V2") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("", "C1", "piecewise(dose / duration, time < duration, 0)", "infusion", rateCompartment = "unit") |>
  addReaction("C1", "", "CL * C1", "elimination", rateCompartment = "unit") |>
  addReaction("C1", "C2", "Q * C1", "distribution", rateCompartment = "unit") |>
  addReaction("C2", "C1", "Q * C2", "redistribution", rateCompartment = "unit")
mavoNet
# solver tolerances 1e-10: at the defaults the solver error, absolute in mg/L, moves -2 log L by
# about 1e-3, more than the fit resolves
mavoX <- Xs(odemodel(mavoNet, modelname = "pm_mavo", compile = FALSE, derivMode = dm),
            options = list(atol = 1e-10, rtol = 1e-10))
# the dose in mg, the concentration in ng/mL
mavoG <- Y(c(cp = "1000 * C1"), mavoX, modelname = "pm_mavo_g", compile = FALSE, derivMode = dm)
mavoE <- Y(c(cp = "sqrt(sigma_abs^2 + (sigma_rel * cp)^2)"), mavoG, attachInput = FALSE, modelname = "pm_mavo_e",
           compile = FALSE, derivMode = dm)
mavoTable <- attr(mavoData, "condition.grid")
mavoTable$ID <- as.character(mavoTable$ID)
mavoTrafo <- function(between) {
  tr <- branch(define(NULL, "x <- x", x = c("CL", "V1", "Q", "V2", "C1", "C2", "dose", "duration", "sigma_abs",
                                                "sigma_rel")),
               table = mavoTable)
  tr <- insert(tr, "dose <- d", d = dose)
  tr <- insert(tr, "duration <- u", u = duration)
  tr <- insert(tr, "x <- 0", x = c("C1", "C2"))
  tr <- insert(tr, "x <- x * (w / 70)^0.75", x = c("CL", "Q"), w = WT)
  tr <- insert(tr, "x <- x * w / 70", x = c("V1", "V2"), w = WT)
  tr <- insert(tr, "x <- exp(x)", x = c("CL", "V1", "Q", "V2", "sigma_abs", "sigma_rel"))
  if (between) tr <- insert(tr, CL ~ 1 + (1 | ID) + (1 | ID:OCC)) else tr <- insert(tr, CL ~ 1 + (1 | ID))
  insert(tr, V1 + Q + V2 ~ 1 + (1 | ID))
}
mavoP <- P(mavoTrafo(TRUE), method = "explicit", modelname = "pm_mavo_p", compile = FALSE,
           derivMode = dm)
mavoP1 <- P(mavoTrafo(FALSE), method = "explicit", modelname = "pm_mavo_p1", compile = FALSE,
            derivMode = dm)
compileSeconds["mavoglurant"] <- system.time(compile(mavoX, mavoG, mavoE, mavoP, mavoP1,
                                                     output = "pm_mavo_all"))[["elapsed"]]
mavoObj <- normL2(mavoData, mavoG * mavoX * mavoP, errmodel = mavoE)
mavoObj1 <- normL2(mavoData, mavoG * mavoX * mavoP1, errmodel = mavoE)
mavoRandom <- randomEffects(mavoP)
mavoRandom1 <- randomEffects(mavoP1)
mavoRandom
mavoInit <- c(CL = log(30), V1 = log(40), Q = log(30), V2 = log(150), sigma_abs = log(1), sigma_rel = log(0.2),
              omega_CL_CL = log(0.3),
              omega_V1_V1 = log(0.3), omega_Q_Q = log(0.3), omega_V2_V2 = log(0.3))
mavoStarts <- drawStarts(c(mavoInit, omega_CL_ID_OCC_CL_ID_OCC = log(0.15)))
mavoStarts1 <- mavoStarts[, names(mavoInit)]
mavoRuns <- cached("Mavoglurant dMod2", list(focei = runDMod2(mavoObj, mavoRandom, mavoStarts, "focei"),
                                              saem = runDMod2(mavoObj1, mavoRandom1, mavoStarts1, "saem")))
mavoFit <- cachedEm("Mavoglurant dMod2 focei best", mavoObj, mavoRandom,
                    mavoStarts[which.min(mavoRuns$focei$value), ], verbose = FALSE)
summary(mavoFit)
compileSeconds

# Diagnostics ----
plot(mavoFit)
plotResiduals(mavoFit, nsim = 200, seed = 1)
plotCombined(mavoFit, conditionsPerPage = 12)[[1]]
plotRandomEffects(mavoFit)
plotRandomEffects(mavoFit, covariates = c("WT", "AGE", "SEX", "HT"))
plotVPC(mavoFit, nsim = 200, bins = 8, pc = TRUE, strata = "dose", seed = 1)

# FOCEI and SAEM ----
compareMethods(mavoRuns, "Mavoglurant", list(obj = mavoObj1, random = mavoRandom1, starts = mavoStarts1))
