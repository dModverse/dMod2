# STAT5 dimerisation after Epo stimulation, from reactions to profile likelihoods:
# multi-start fit, an error model by restricted maximum likelihood, confidence
# intervals and a cross-check against the shipped PEtab form.
# Boehm et al. (2014) J Proteome Res, doi:10.1021/pr5006923.

library(dMod2)
library(ggplot2)

.modelname <- "boehm"
.outdir    <- file.path(tempdir(), .modelname)
if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
# generated sources and shared objects go to .outdir, see ?dMod2 for the options;
# compiling, fitting and profiling use the free cores of the local machine
.opt   <- options(dMod2.outdir = .outdir)
.cores <- detectFreeCores()
set.seed(5555)
.petabDir <- system.file("extdata", "petab_boehm", package = "dMod2")


# Reactions ---------------------------------------------------------------------
# BaF3-EpoR cells respond to erythropoietin: STAT5A and STAT5B are phosphorylated,
# form the dimers ApA, ApB and BpB, move into the nucleus and return. The Epo
# stimulus decays and enters the rates in closed form, so `time` appears in a rate.
# Nuclear species get their compartment before the reactions that produce them.
epo <- "1.25e-7*exp(-Epo_degradation_BaF3*time)"

reactions <- eqnlist() |>
  assignCompartment(nucpApA = "nuc", nucpApB = "nuc", nucpBpB = "nuc", volume = "0.45") |>
  addReaction("2*STAT5A", "pApA", paste0("k_phos*", epo, "*STAT5A^2"),
              "STAT5A-STAT5A phosphorylation", compartment = "cyt") |>
  addReaction("STAT5A + STAT5B", "pApB", paste0("k_phos*", epo, "*STAT5A*STAT5B"),
              "STAT5A-STAT5B phosphorylation", compartment = "cyt") |>
  addReaction("2*STAT5B", "pBpB", paste0("k_phos*", epo, "*STAT5B^2"),
              "STAT5B-STAT5B phosphorylation", compartment = "cyt") |>
  addReaction("pApA", "nucpApA", "k_imp_homo*pApA", "pApA import", compartment = "cyt") |>
  addReaction("pApB", "nucpApB", "k_imp_hetero*pApB", "pApB import", compartment = "cyt") |>
  addReaction("pBpB", "nucpBpB", "k_imp_homo*pBpB", "pBpB import", compartment = "cyt") |>
  addReaction("nucpApA", "2*STAT5A", "k_exp_homo*nucpApA", "pApA export", compartment = "nuc") |>
  addReaction("nucpApB", "STAT5A + STAT5B", "k_exp_hetero*nucpApB", "pApB export", compartment = "nuc") |>
  addReaction("nucpBpB", "2*STAT5B", "k_exp_homo*nucpBpB", "pBpB export", compartment = "nuc") |>
  setCompartmentVolume(cyt = "1.4")

# The volumes enter the right-hand side as the ratios that convert a flux from the
# compartment it is written in into that of the species it changes.
as.eqnvec(reactions)


# Prediction, observation and error model ---------------------------------------
# Only relative quantities were measured, mixed by the isotope ratio specC17. The
# error model has one standard deviation per observable, estimated with the rest.
myOptions <- list(atol = 1e-8, rtol = 1e-6, maxattempts = 100L, maxsteps = 1e6)

model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE)
x <- Xs(model, options = myOptions)

observables <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))"
)
errors <- eqnvec(
  pSTAT5A_rel = "sd_pSTAT5A_rel",
  pSTAT5B_rel = "sd_pSTAT5B_rel",
  rSTAT5A_rel = "sd_rSTAT5A_rel"
)

g <- Y(observables, x, modelname = "boehm_obs", compile = FALSE)
e <- Y(errors, g, modelname = "boehm_err", compile = FALSE)


# Parameter transformation ------------------------------------------------------
# Initial dimers are zero, the STAT5 pool is split by the measured ratio, and the
# nine estimated parameters are on log10.
innerpars <- getParameters(model, g, e)
estimated <- c("Epo_degradation_BaF3", "k_exp_hetero", "k_exp_homo", "k_imp_hetero",
               "k_imp_homo", "k_phos", "sd_pSTAT5A_rel", "sd_pSTAT5B_rel", "sd_rSTAT5A_rel")

p <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("x~0", x = c("pApA", "pApB", "pBpB", "nucpApA", "nucpApB", "nucpBpB")) |>
  insert("STAT5A ~ 207.6*ratio") |>
  insert("STAT5B ~ 207.6 - 207.6*ratio") |>
  insert("ratio ~ 0.693") |>
  insert("specC17 ~ 0.107") |>
  insert("x~exp10(x)", x = estimated) |>
  P(modelname = "boehm_trafo", condition = "Boehm2014", compile = FALSE)

# One call compiles ODE, sensitivities, observation, error model and trafo
# together into one shared object.
compile(g, x, p, e, output = .modelname, cores = .cores)
prd       <- g*x*p
outerpars <- getParameters(prd)


# Objective and multi-start fit -------------------------------------------------
# Two directions are practically non-identifiable, so a weak prior keeps the
# dynamic parameters finite. It stays off the standard deviations, whose REML
# estimate below solves the stationarity condition of the data term alone.
data(boehm)
data <- as.datalist(boehm)

pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
dyn    <- setdiff(outerpars, grep("^sd_", outerpars, value = TRUE))
obj    <- normL2(data, prd, e) + constraintL2(pouter[dyn], sigma = 4)

# Fifty trust-region fits from random starts, sorted by value: the lowest
# plateau holds the runs that reached the best optimum.
fits <- mstrust(obj, pouter, sd = 3, fits = 50, cores = .cores, iterlim = 5000) |>
  as.parframe()
plotValues(fits, tol = 0.1, value < 1e4)
best <- as.parvec(fits)


# Error model by restricted maximum likelihood ----------------------------------
# Maximum likelihood underestimates the standard deviations, because the fitted
# mean shrinks the residuals. `reml()` counts each data point with one minus its
# leverage, the stationarity condition of the restricted likelihood (Patterson
# and Thompson 1971), which corrects the variances for the degrees of freedom the
# mean uses.
remlfit <- reml(obj, best)
stopifnot(remlfit$converged)
rbind(ML = best[remlfit$errpars], REML = remlfit$argument[remlfit$errpars])
remlfit$dof     # data points minus the summed leverage, per observable
remlfit$rank    # effective number of mean parameters

# The summed leverage per observable adds up to the rank.
tapply(remlfit$leverage$leverage, remlfit$leverage$name, sum)


# Fit and uncertainty band ------------------------------------------------------
# One standard deviation of the REML error model around the prediction.
times      <- seq(0, 240, length.out = 200)
prediction <- prd(times, remlfit$argument, deriv = FALSE)

plot(prediction, data) +
  geom_ribbon(data = as.data.frame(prediction, errfn = e),
              aes(x = time, ymin = value - sigma, ymax = value + sigma),
              linetype = "dashed", alpha = 0.3) +
  labs(x = "time [min]", y = "relative signal [%]", colour = NULL)


# Profile likelihood ------------------------------------------------------------
# A profile re-optimises every other parameter while one is held away from the
# optimum, so the interval comes from the likelihood itself (Raue et al. 2009).
# `profileThreshold()` gives the threshold for `profile()` and `plotProfile()`;
# `confint()` computes the same one from `method`, `n` and `p`. `method = "F"`
# calibrates it for the finite sample, with the effective rank of the REML fit as
# the number of mean parameters.
nd     <- sum(vapply(data, nrow, 0L))
thr_95 <- profileThreshold(0.95, method = "F", n = nd, p = remlfit$rank)
thr_90 <- profileThreshold(0.90, method = "F", n = nd, p = remlfit$rank)

# The intervals are taken on the data term; the prior only keeps the fit finite.
profiles <- profile(
  obj, best, whichPar = names(best), limits = c(lower = -5, upper = 5),
  cores = .cores, delta = thr_95,
  stepControl = list(atol = 1e-3, rtol = 1e-3, limit = 200, stop = "data"),
  algoControl = list(reoptimize = TRUE),
  optControl  = list(iterlim = 20))

plotProfile(profiles, mode %in% c("data", "prior"),
            threshold = c("90%" = thr_90, "95%" = thr_95))
confint(profiles, level = 0.95, val.column = "data",
        method = "F", n = nd, p = remlfit$rank)

# For k_imp_homo and k_exp_hetero the data term stays level and only the prior
# rises: homodimer import is faster than phosphorylation, and the heterodimer
# barely returns within 240 minutes. Those intervals stay open on one side.


# Cross-check against the PEtab form of the same problem ------------------------
# The import rebuilds the model from SBML and the tables. At the published
# optimum, `petab$bestfit`, both forms give the same sum of squares. The explicit
# modelname keeps the import's files apart from the hand-built ones.
petab <- importPEtab(file.path(.petabDir, "Boehm.yaml"), backend = "cppDE",
                     cores = .cores, modelname = "boehm_petab", options = myOptions)
stopifnot(setequal(names(petab$bestfit), outerpars))
published <- petab$bestfit[outerpars]
chi2 <- c(hand  = attr(obj(published, deriv = FALSE), "chi2"),
          PEtab = attr(petab$obj(published, deriv = FALSE), "chi2"))
chi2
stopifnot(isTRUE(all.equal(unname(chi2[1]), unname(chi2[2]), tolerance = 1e-6)))

options(.opt)
