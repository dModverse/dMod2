# -------------------------------------------------------------------------#
# Thermal abuse of a lithium-ion cell: symmetries of an Arrhenius model
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A literature model in which states sit inside exponentials and every
# parameter and initial value is positive. Four decomposition reactions heat a
# cell in an oven; only its temperature is measured, as in an oven or ARC test.
# symmetryDetection() finds seven scalings and one translation of the SEI
# thickness inside exp(-tsei/tref); symmetryReduction() removes all eight with
# a chart certified positive. The runaway is too sharp for a fit without
# multiple shooting, so the bench stops at the analysis.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Thu 24 Sep 2026
#
# [Info]
# Lumped form of Hatchard TD, MacNeil DD, Basu A, Dahn JR (2001). Thermal model
# of cylindrical and prismatic lithium-ion cells. J Electrochem Soc 148, A755,
# with the reaction set and parameters of Kim GH, Pesaran A, Spotnitz R (2007).
# A three-dimensional thermal abuse model for lithium-ion cells. J Power Sources
# 170, 476-489. Activation energies are given in K (E_a/R), the reaction orders
# are 1. Needs dMod2 from devel-symmetry.
#
# Detection and reduction run in the linear coordinates with the positive
# domain declared; the log10 parametrisation is inserted only after the
# reduction, so that its chart is certified on the domain the fit explores.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryBattery"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
#
# csei  metastable SEI fraction          Rsei  SEI decomposition
# cneg  lithium in the anode             Rne   anode-electrolyte reaction, slowed
# tsei  dimensionless SEI thickness            by the SEI it builds, exp(-tsei/tref)
# alpha cathode conversion               Rpe   cathode-electrolyte reaction
# ce    electrolyte                      Re    electrolyte decomposition
# T     cell temperature, heated by all four and cooled towards the oven at Ta
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
Rsei <- "Asei*exp(-Esei/T)*csei"
Rne  <- "Ane*exp(-tsei/tref)*exp(-Ene/T)*cneg"
Rpe  <- "Ape*exp(-Epe/T)*alpha*(1 - alpha)"
Re   <- "Ae*exp(-Ee/T)*ce"

f <- eqnvec(
  csei  = paste0("-", Rsei),
  cneg  = paste0("-", Rne),
  tsei  = Rne,
  alpha = Rpe,
  ce    = paste0("-", Re),
  T     = sprintf("(Hsei*Wc*%s + Hne*Wc*%s + Hpe*Wp*%s + He*We*%s - hA*(T - Ta))/rhocp",
                  Rsei, Rne, Rpe, Re))

observables <- eqnvec(y = "T")

# Kim et al. (2007), Table 1; rhocp and hA (per cell volume) for an 18650 cell
truth <- c(
  csei = 0.15, cneg = 0.75, tsei = 0.033, alpha = 0.04, ce = 1, T = 298,
  Asei = 1.667e15, Esei = 16247, Hsei = 257,
  Ane  = 2.5e13,   Ene  = 16247, Hne  = 1714, tref = 0.033,
  Ape  = 6.667e13, Epe  = 16791, Hpe  = 314,
  Ae   = 5.14e25,  Ee   = 32956, He   = 155,
  Wc = 6.104e5, Wp = 1.221e6, We = 4.07e5,
  rhocp = 2.14e6, hA = 1816, Ta = 428)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Symmetries
#
# Every coordinate is positive, the default. The seven scalings say that heat
# and mass enter only as products H*W*c; the eighth shifts tsei by eps*tref and
# scales Ane by exp(eps), so only Ane*exp(-tsei/tref) reaches the temperature.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
res <- symmetryDetection(f, observables, reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

red <- symmetryReduction(res)
summary(red)
red$trafo

# the reduction is complete: in its chart nothing is left to find
symmetryDetection(f, observables, trafo = red$trafo, verbose = FALSE)$identifiable

# Known material constants are the more physical gauge. They leave the same
# eight directions removed, with other outer parameters.
red.phys <- symmetryReduction(res, fixed = c("Wc", "Wp", "We", "rhocp", "hA"))
red.phys$trafo


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# The ODE is integrated as it stands. Both trafos go to log10 only after the
# reduction: the chart is certified positive, so every outer parameter of it
# can be fitted as exp10().
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptionsODE  <- list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)
myOptionsSens <- myOptionsODE

model <- odemodel(f, modelname = "battery_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "oven", optionsOde = myOptionsODE, optionsSens = myOptionsSens)
g <- Y(observables, f = x, attach.input = TRUE,
       modelname = "battery_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)

# every coordinate free; the oven temperature is set by the experiment
p.full <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("Ta~428") |>
  insert("x~exp10(x)", x = setdiff(.currentSymbols, "Ta")) |>
  P(modelname = "battery_pfull", compile = FALSE, outdir = .outdir)

p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  define("Ta~428") |>
  insert("x~exp10(x)", x = setdiff(.currentSymbols, "Ta")) |>
  P(modelname = "battery_pred", compile = FALSE, outdir = .outdir)

# One flow per direction, dz/deps = sgn*xi(z), in the log10 coordinates of the
# moved coordinates
flow <- lapply(seq_along(res$symmetries), function(i) {
  fl <- log10Transform(res$symmetries[[i]]$completeGenerator)
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
              deriv = FALSE, modelname = paste0("battery_flow", i), compile = FALSE,
              outdir = .outdir),
     condition = "flow",
     optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
})

do.call(compile, c(list(x, g, p.full, p.red), flow,
                   list(output = .modelname, cores = .cores)))

prd     <- g*x*p.full
prd.red <- g*x*p.red

truthL <- log10(truth[getParameters(prd)])
times  <- seq(0, 4000, len = 400)
pred0  <- as.data.frame(prd(times, truthL, deriv = FALSE))

ggplot(pred0, aes(time, value)) + geom_line() +
  facet_wrap(~ name, scales = "free_y") + theme_dMod(base_size = 9)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Every direction, both ways
#
# Along each flow the states move and T does not. The deviation is relative to
# the RMS of the unmoved temperature.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
eps <- seq(0, 1, len = 11)
obs <- pred0$name %in% names(observables)
rms <- sqrt(mean(pred0$value[obs]^2))
tl  <- setNames(log10(truth), paste0(names(truth), "_l10"))

for (i in seq_along(flow)) {
  mv  <- names(res$symmetries[[i]]$generator)
  mvL <- paste0(mv, "_l10")
  z0  <- c(tl[mvL], truth[setdiff(coords, mv)])
  for (sgn in c(1, -1)) {
    zi <- try(flow[[i]](eps, c(z0, sgn = sgn), deriv = FALSE)$flow, silent = TRUE)
    if (inherits(zi, "try-error")) {
      cat(sprintf("X%d  eps %+.2f  escapes\n", i, sgn * max(eps)))
      next
    }
    d <- sapply(seq_len(nrow(zi)), function(j) {
      q <- as.data.frame(prd(times, replace(truthL, mv, zi[j, mvL]), deriv = FALSE))
      max(abs(q$value[obs] - pred0$value[obs])) / rms
    })
    cat(sprintf("X%d  %-8s on %-30s  eps %+.2f  moved %.2f  deviates %.1e\n", i,
                res$symmetries[[i]]$type,
                paste(res$symmetries[[i]]$support, collapse = ","),
                sgn * max(zi[, "time"]),
                max(abs(10^(zi[nrow(zi), mvL] - tl[mvL]) - 1)), max(d)))
  }
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Reduced truth
#
# survivorMeaning names the invariant each outer parameter carries; evaluated
# at the truth it gives the reduced truth, whose prediction is the full one.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
meaning <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
outer   <- getParameters(prd.red)
vals    <- c(truth, sapply(meaning, function(m) eval(parse(text = m), as.list(truth))))
truth.red  <- vals[outer]
truthL.red <- log10(truth.red)

predR <- as.data.frame(prd.red(times, truthL.red, deriv = FALSE))
max(abs(predR$value[obs] - pred0$value[obs])) / rms
