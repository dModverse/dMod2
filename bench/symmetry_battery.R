# -------------------------------------------------------------------------#
# Thermal abuse of a lithium-ion cell: symmetries of an Arrhenius model
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A literature model with states inside exponentials and every parameter and
# initial value positive. Four decomposition reactions heat a cell in an oven;
# only its temperature is observed, as in an oven or ARC test.
# symmetryDetection() finds the scalings and a translation of the SEI
# thickness inside exp(-tsei/tref), symmetryReduction() removes them in a
# chart certified positive. Each direction is flown both ways against the
# prediction, and the reduced truth is checked against the full one. The
# script stops at the analysis and does not fit.
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
# are 1.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryBattery"
.outdir    <- file.path(tempdir(), .modelname)  # all generated sources and libraries
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# ---- Model ------------------------------------------------------------------

Rsei <- "Asei*exp(-Esei/T)*csei"                    # SEI decomposition
Rne  <- "Ane*exp(-tsei/tref)*exp(-Ene/T)*cneg"      # anode-electrolyte, slowed by the SEI
Rpe  <- "Ape*exp(-Epe/T)*alpha*(1 - alpha)"         # cathode-electrolyte
Re   <- "Ae*exp(-Ee/T)*ce"                          # electrolyte decomposition

f <- eqnvec(
  csei  = paste0("-", Rsei),                        # metastable SEI fraction
  cneg  = paste0("-", Rne),                         # lithium in the anode
  tsei  = Rne,                                      # dimensionless SEI thickness
  alpha = Rpe,                                      # cathode conversion
  ce    = paste0("-", Re),                          # electrolyte
  T     = sprintf("(Hsei*Wc*%s + Hne*Wc*%s + Hpe*Wp*%s + He*We*%s - hA*(T - Ta))/rhocp",
                  Rsei, Rne, Rpe, Re))              # cell temperature, oven at Ta

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


# ---- Symmetries -------------------------------------------------------------

# Detection and reduction run in the linear coordinates with every coordinate
# positive, the default. reconstruct = TRUE returns the exact generators.
res <- symmetryDetection(f, observables, reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

red <- symmetryReduction(res)
summary(red)
red$trafo

# a second detection in the reduced chart: TRUE when nothing is left
check.red <- symmetryDetection(f, observables, trafo = red$trafo, verbose = FALSE)
check.red$identifiable

# the same reduction gauged by the known material constants
red.phys <- symmetryReduction(res, fixed = c("Wc", "Wp", "We", "rhocp", "hA"))
red.phys$trafo


# ---- Prediction chain -------------------------------------------------------

# The ODE is integrated in linear coordinates. Both parameter trafos switch to
# log10 after the reduction, since its chart is certified positive.
myOptions <- list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)

model <- odemodel(f, modelname = "battery_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "oven", options = myOptions)
g <- Y(observables, f = x, attachInput = TRUE,
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

# one flow per direction, dz/deps = sgn*xi(z), in log10 of the moved coordinates
flow <- lapply(seq_along(res$symmetries), function(i) {
  fl <- log10Transform(res$symmetries[[i]]$completeGenerator)
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
              deriv = FALSE, modelname = paste0("battery_flow", i), compile = FALSE,
              outdir = .outdir),
     condition = "flow",
     options = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
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


# ---- Every direction, both ways ---------------------------------------------

# Along each flow the states move and the observable does not. The deviation
# is relative to the RMS of the unmoved observable; NA marks a flow that escapes.
eps <- seq(0, 1, len = 11)
obs <- pred0$name %in% names(observables)
rms <- sqrt(mean(pred0$value[obs]^2))
tl  <- setNames(log10(truth), paste0(names(truth), "_l10"))

runs <- expand.grid(i = seq_along(flow), sgn = c(1, -1))

flowcheck <- do.call(rbind, lapply(seq_len(nrow(runs)), function(k) {
  i   <- runs$i[k]
  sgn <- runs$sgn[k]
  sym <- res$symmetries[[i]]
  mv  <- names(sym$generator)
  mvL <- paste0(mv, "_l10")
  z0  <- c(tl[mvL], truth[setdiff(coords, mv)])
  row <- data.frame(direction = paste0("X", i), type = sym$type,
                    support = paste(sym$support, collapse = ","),
                    eps = NA_real_, moved = NA_real_, deviates = NA_real_)
  zi <- try(flow[[i]](eps, c(z0, sgn = sgn), deriv = FALSE)$flow, silent = TRUE)
  if (inherits(zi, "try-error")) return(row)
  d <- sapply(seq_len(nrow(zi)), function(j) {
    q <- as.data.frame(prd(times, replace(truthL, mv, zi[j, mvL]), deriv = FALSE))
    max(abs(q$value[obs] - pred0$value[obs])) / rms
  })
  row$eps      <- sgn * max(zi[, "time"])
  row$moved    <- max(abs(10^(zi[nrow(zi), mvL] - tl[mvL]) - 1))
  row$deviates <- max(d)
  row
}))
flowcheck


# ---- Reduced truth ----------------------------------------------------------

# survivorMeaning names the invariant behind each outer parameter; evaluated at
# the truth it gives the reduced truth, whose prediction should equal the full one.
meaning <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
outer   <- getParameters(prd.red)
vals    <- c(truth, sapply(meaning, function(m) eval(parse(text = m), as.list(truth))))
truth.red  <- vals[outer]
truthL.red <- log10(truth.red)

predR <- as.data.frame(prd.red(times, truthL.red, deriv = FALSE))
dev.red <- max(abs(predR$value[obs] - pred0$value[obs])) / rms
dev.red
