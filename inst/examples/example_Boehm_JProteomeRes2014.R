# STAT5 dimerisation after Epo stimulation: multi-start comparison of the gn and
# bfgs Hessian sources, reverse-mode gradient, PEtab cross-check.
# Boehm et al. (2014) J Proteome Res; PEtab copy in inst/extdata/petab_boehm.

library(dMod2)
library(ggplot2)

.modelname <- "boehm"
# generated sources, objects and libraries go here, not the working directory
.outdir    <- file.path(tempdir(), .modelname)
.fit       <- TRUE   # multi-start comparison

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
.petabDir <- system.file("extdata", "petab_boehm", package = "dMod2")


# Load data ---------------------------------------------------------------------
# Standard deviations are estimated: `sigma` is NA and the error model supplies it.
data(boehm)
mydataL <- as.datalist(boehm)

cat(sprintf("%d points, %d condition, %d observables\n", nrow(boehm),
            length(mydataL), length(unique(boehm$name))))


# Model -------------------------------------------------------------------------
# The Epo stimulus enters the rates in closed form, so `time` appears in a rate.
# Nuclear species are assigned their compartment before they are produced.
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

myOptions <- list(atol = 1e-8, rtol = 1e-6, maxattempts = 100L, maxsteps = 1e6)

# `derivMode = c("forward", "reverse")` adds the checkpointed value solve and
# backward sweep the reverse section needs.
model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE,
                  derivMode = c("forward", "reverse"), outdir = .outdir)
x <- Xs(model, options = myOptions)

# Only relative quantities were measured, mixed by the isotope ratio specC17.
observables <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))"
)
errorModels <- eqnvec(
  pSTAT5A_rel = "sd_pSTAT5A_rel",
  pSTAT5B_rel = "sd_pSTAT5B_rel",
  rSTAT5A_rel = "sd_rSTAT5A_rel"
)

g <- Y(observables, x, modelname = "boehm_obs", attachInput = FALSE,
       compile = FALSE, outdir = .outdir)
e <- Y(errorModels, g, modelname = "boehm_err", attachInput = FALSE,
       compile = FALSE, outdir = .outdir)


# Parameter transformation ------------------------------------------------------
# Initial dimers are zero, the STAT5 pool is split by the measured ratio, and
# every estimated parameter is on log10.
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
  P(modelname = "boehm_trafo", condition = "Boehm2014", compile = FALSE, outdir = .outdir)

# One call compiles ODE, sensitivities, observation, error and trafo together.
compile(g, x, p, e, output = .modelname, cores = 6)

prd       <- g*x*p
outerpars <- getParameters(prd)


# Objective ---------------------------------------------------------------------
# A weak prior on the dynamic parameters keeps the non-identifiable directions
# finite; bounds and the published optimum come from the PEtab parameter table.
.pars <- read.delim(file.path(.petabDir, "parameters_Boehm_JProteomeRes2014.tsv"))
.pars <- .pars[.pars$estimate == 1, ]
.lower <- setNames(log10(.pars$lowerBound), .pars$parameterId)[outerpars]
.upper <- setNames(log10(.pars$upperBound), .pars$parameterId)[outerpars]

pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
dyn    <- setdiff(outerpars, grep("^sd_", outerpars, value = TRUE))
obj    <- normL2(mydataL, prd, e) + constraintL2(pouter[dyn], sigma = 4)

# the published optimum, on the log10 scale the table declares
bestfit <- setNames(log10(.pars$nominalValue), .pars$parameterId)[outerpars]
stopifnot(setequal(names(bestfit), outerpars))


# Multi-start fit: gn vs bfgs vs a gn run handing over to bfgs ------------------
# All methods share one set of starts and one evaluation limit, so they differ
# only in their Hessian source; qnEval counts quasi-Newton evaluations.
if (.fit) {
  # A run is a Hessian source and, optionally, a fallback it hands over to.
  .methods <- list(gn        = list(hessianMethod = "gn"),
                   bfgs      = list(hessianMethod = "bfgs"),
                   `gn->bfgs` = list(hessianFallback = "bfgs"))
  set.seed(20260905)
  .starts  <- msParframe(pouter, n = 100, sd = 3)   # shared across methods

  runs <- lapply(.methods, function(a)
    do.call(mstrust, c(list(obj, center = .starts, fits = nrow(.starts), cores = 20,
                            rinit = 0.1, rmax = 10, iterlim = 1500,
                            parlower = .lower, parupper = .upper), a)))
  frames <- lapply(runs, as.parframe)

  # The best value over all methods is the reference; a start succeeds within
  # 0.1 of it, and `converged` separates a miss from a run that never stops.
  .best <- min(vapply(frames, function(f) min(f$value), 0.0))
  summary <- do.call(rbind, lapply(.methods, function(hm) {
    f  <- frames[[hm]]
    ok <- f$value <= .best + 0.1
    data.frame(method       = hm,
               fits         = nrow(f),
               converged    = sum(f$converged),
               success_rate = mean(ok),
               median_neval = stats::median(f$neval[ok]),
               total_neval  = sum(f$neval),
               qn_share     = sum(f$qnEval) / sum(f$neval),
               best_value   = min(f$value))
  }))
  print(summary, row.names = FALSE)

  lapply(.methods, function(hm) {
      outframe <- frames[[hm]]
      bestfit  <- as.parvec(outframe)
      print(plotValues(outframe, tol = 0.1, value < 1e4))
    })
}


# The gradient the other way round ----------------------------------------------
# Reverse integrates the states alone and sweeps one tape back. It agrees with
# forward up to O(tol) and returns no Hessian, see example_ReverseAD.R.
.fwd <- obj(bestfit, deriv = TRUE)
.rev <- obj(bestfit, deriv = TRUE, sweep = "reverse")

print(rbind(forward = .fwd$gradient,
            reverse = .rev$gradient[names(.fwd$gradient)]))
cat(sprintf("value    %.10g vs %.10g\ngradient max relative difference %.2e\n",
            .fwd$value, .rev$value,
            max(abs(.fwd$gradient - .rev$gradient[names(.fwd$gradient)])) /
              max(abs(.fwd$gradient))))
cat("the reverse objective has no Hessian:", is.null(.rev$hessian), "\n")

# Forward cost rises with n_theta and reverse cost does not, so which side wins
# is a property of the problem.
.reps <- 10
.tf <- system.time(for (i in seq_len(.reps)) obj(bestfit, deriv = TRUE))[["elapsed"]]
.tr <- system.time(for (i in seq_len(.reps))
                     obj(bestfit, deriv = TRUE, sweep = "reverse"))[["elapsed"]]
cat(sprintf("%d gradients at %d parameters: forward %.2fs, reverse %.2fs\n",
            .reps, length(outerpars), .tf, .tr))


# Fit and uncertainty band ------------------------------------------------------
times  <- seq(0, 240, length.out = 200)
prdout <- as.data.frame(prd(times, bestfit, deriv = FALSE), errfn = e)

plot(prd(times, bestfit, deriv = FALSE), mydataL) +
  geom_ribbon(data = prdout, aes(x = time, ymin = value - sigma, ymax = value + sigma),
              linetype = "dashed", alpha = 0.2) +
  labs(x = "time [min]", y = "relative signal [%]", colour = NULL, fill = NULL,
       title = "STAT5 dimerisation, Boehm 2014")


# Cross-check against the PEtab form of the same problem ------------------------
# The explicit modelname keeps the import's files apart from the hand-built ones,
# which the YAML basename would clash with on a case-insensitive filesystem.
petab <- importPEtab(file.path(.petabDir, "Boehm.yaml"),
                     backend = "cppDE", cores = 4, modelname = "boehm_petab",
                     outdir = .outdir,
                     options = myOptions)

stopifnot(setequal(names(petab$bestfit), outerpars))
chi2 <- c(hand  = attr(obj(bestfit, deriv = FALSE),       "chi2"),
          PEtab = attr(petab$obj(bestfit, deriv = FALSE), "chi2"))
chi2
