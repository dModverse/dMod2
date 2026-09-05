# -------------------------------------------------------------------------#
# STAT5 dimerisation after Epo stimulation
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Boehm et al. (2014) in dMod2: STAT5A/STAT5B are phosphorylated, form the
# dimers ApA, ApB and BpB, shuttle into the nucleus and back. 48 measurements,
# one condition, 9 estimated parameters. The multi-start section fits the model
# with the three trust-region Hessian sources gn, bfgs and hybrid over one
# shared set of starting points.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Fri 05 Sep 2026
#
# [Info]
# The data ship with the package as `boehm`, the PEtab form as
# inst/extdata/petab_boehm. The last section imports that one and checks the two
# against each other. Parameters are on log10, so the coordinates are already
# comparable and the trust region needs no parscale.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "boehm"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.fit       <- FALSE   # multi-start comparison, minutes on this model

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
.petabDir <- system.file("extdata", "petab_boehm", package = "dMod2")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Load data
#
# Relative quantities from mass spectrometry, one condition (Boehm2014). The
# standard deviations are estimated, so `sigma` is NA and the error model below
# supplies it.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
data(boehm)
mydataL <- as.datalist(boehm)

cat(sprintf("%d points, %d condition, %d observables\n", nrow(boehm),
            length(mydataL), length(unique(boehm$name))))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
#
# The Epo stimulus decays exponentially and enters the phosphorylation rates in
# closed form, so `time` appears directly in a rate. Nuclear species are
# assigned their compartment before they are produced.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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

myOptionsODE  <- list(atol = 1e-8, rtol = 1e-6, maxattemps = 100L, maxsteps = 1e6)
myOptionsSens <- myOptionsODE

model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, optionsOde = myOptionsODE, optionsSens = myOptionsSens)

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

g <- Y(observables, x, modelname = "boehm_obs", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)
e <- Y(errorModels, g, modelname = "boehm_err", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameter transformation
#
# Fix the initial dimers to zero, split the total STAT5 pool by the measured
# ratio, and put every estimated parameter on log10.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Objective
#
# Two directions are practically non-identifiable, so a weak prior on the
# dynamic parameters keeps them finite; it stays off the sd_*. Bounds and the
# published optimum come from the PEtab parameter table.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Multi-start fit: gn vs bfgs vs hybrid
#
# One set of 100 starting points, drawn once with msParframe() and handed to
# every method as the `center`, so the three fits differ only in their Hessian
# source. neval counts objective (gradient) evaluations.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
if (.fit) {
  .methods <- c("gn", "bfgs", "hybrid")
  .starts  <- msParframe(pouter, n = 100, sd = 3)   # shared across methods

  runs <- lapply(.methods, function(hm)
    mstrust(obj, center = .starts, fits = nrow(.starts), cores = 6,
            rinit = 0.1, rmax = 10, iterlim = 5000,
            parlower = .lower, parupper = .upper, hessianMethod = hm))
  names(runs) <- .methods
  frames <- lapply(runs, as.parframe)

  # Best value seen anywhere is the reference optimum; a start "succeeds" if it
  # lands within 0.1 of it.
  .best <- min(vapply(frames, function(f) min(f$value), 0.0))
  summary <- do.call(rbind, lapply(.methods, function(hm) {
    f  <- frames[[hm]]
    ok <- f$value <= .best + 0.1
    data.frame(method       = hm,
               fits         = nrow(f),
               success_rate = mean(ok),
               median_neval = stats::median(f$neval[ok]),
               total_neval  = sum(f$neval),
               best_value   = min(f$value))
  }))
  print(summary, row.names = FALSE)

  outframe <- frames[["hybrid"]]
  bestfit  <- as.parvec(outframe)
  print(plotValues(outframe, tol = 0.1, value < 1e4))
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and uncertainty band
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
times  <- seq(0, 240, length.out = 200)
prdout <- as.data.frame(prd(times, bestfit, deriv = FALSE), errfn = e)

plot(prd(times, bestfit, deriv = FALSE), mydataL) +
  geom_ribbon(data = prdout, aes(x = time, ymin = value - sigma, ymax = value + sigma),
              linetype = "dashed", alpha = 0.2) +
  labs(x = "time [min]", y = "relative signal [%]", colour = NULL, fill = NULL,
       title = "STAT5 dimerisation, Boehm 2014")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Cross-check against the PEtab form of the same problem
#
# The benchmark collection distributes Boehm as PEtab, and the package ships
# that copy. Importing it builds the model a second time, from SBML and the
# tables instead of from the reactions above. An explicit modelname keeps its
# sources and shared object apart from the hand-built "boehm" ones -- the YAML
# basename would default to "Boehm", which collides on a case-insensitive
# filesystem.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
petab <- importPEtab(file.path(.petabDir, "Boehm.yaml"),
                     backend = "cppDE", cores = 4, modelname = "boehm_petab",
                     outdir = .outdir,
                     optionsOde = myOptionsODE, optionsSens = myOptionsSens)

stopifnot(setequal(names(petab$bestfit), outerpars))
chi2 <- c(hand  = attr(obj(bestfit, deriv = FALSE),       "chi2"),
          PEtab = attr(petab$obj(bestfit, deriv = FALSE), "chi2"))
chi2
