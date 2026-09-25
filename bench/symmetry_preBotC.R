# -------------------------------------------------------------------------#
# Pre-Boetzinger pacemaker neuron: symmetries on a mixed-sign domain
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A literature model without a global positivity constraint: the membrane
# potential, the reversal potentials and the half-activation voltages take
# either sign, and so do the slopes of the Boltzmann gates. Only the
# capacitance, the conductances, the time constants and the gates are positive.
# Voltage is measured. symmetryDetection() finds one scaling and two general
# directions, symmetryReduction() removes all three on the declared domain. The
# last part compares with exp10() in the trafo handed to the detection. A
# spiking model is not fitted here: without multiple shooting the spike phase
# makes the objective too rugged.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Thu 24 Sep 2026
#
# [Info]
# Model 1 of Butera RJ, Rinzel J, Smith JC (1999). Models of respiratory
# rhythm generation in the pre-Boetzinger complex. I. Bursting pacemaker
# neurons. J Neurophysiol 82, 382-397. Gates x_inf = 1/(1 + exp((V - thx)/sx)),
# time constants taux/cosh((V - thx)/(2*sx)); m and p are instantaneous.
# Units mV, ms, pF, nS. Needs dMod2 from devel-symmetry.
#
# The log10 parametrisation belongs after the reduction, not into the trafo of
# the detection: the reduction certifies its chart on the declared domain in
# linear coordinates, and exp10() then only reparametrises what is known to be
# positive. The comparison at the end reaches the same chart the long way.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryPreBotC"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
#
# INaP persistent sodium, gated by p (instantaneous) and h (slow inactivation)
# INa, IK fast spike currents, gated by m (instantaneous) and n
# IL, Iton leak and tonic excitatory drive
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
xinf <- function(x) sprintf("1/(1 + exp((V - th%s)/s%s))", x, x)
taux <- function(x) sprintf("tau%s/cosh((V - th%s)/(2*s%s))", x, x, x)

f <- eqnvec(
  V = sprintf(paste0("-(gNaP*%s*h*(V - ENa) + gNa*(%s)^3*(1 - n)*(V - ENa) + ",
                     "gK*n^4*(V - EK) + gL*(V - EL) + gton*(V - Esyn))/C"),
              xinf("p"), xinf("m")),
  n = sprintf("(%s - n)/(%s)", xinf("n"), taux("n")),
  h = sprintf("(%s - h)/(%s)", xinf("h"), taux("h")))

observables <- eqnvec(y = "V")

positive <- c("C", "gNaP", "gNa", "gK", "gL", "gton", "taun", "tauh", "n", "h")

# Butera et al. (1999), Table 1, gNaP = 2.8 nS (bursting) and gton = 0.3 nS
truth <- c(
  V = -60, n = 0.01, h = 0.6,
  C = 21, gNaP = 2.8, gNa = 28, gK = 11.2, gL = 2.8, gton = 0.3,
  ENa = 50, EK = -85, EL = -65, Esyn = 0,
  thm = -34, sm = -5, thn = -29, sn = -4, thp = -40, sp = -6, thh = -48, sh = 6,
  taun = 10, tauh = 10000)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Symmetries on the declared domain
#
# V alone cannot tell the capacitance from the conductances (X1), and the two
# ohmic currents IL and Iton act as one: only gL + gton and gL*EL + gton*Esyn
# enter (X2, X3).
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
res <- symmetryDetection(f, observables, positive = positive, reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

red <- symmetryReduction(res, positive = positive)
summary(red)
red$trafo

symmetryDetection(f, observables, trafo = red$trafo, positive = positive,
                  verbose = FALSE)$identifiable


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# Positive coordinates and positive outer parameters of the reduction go to
# log10, everything real-valued stays linear.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# spikes last about 2 ms; the check of the reduced chart needs tight solves
myOptionsODE  <- list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)
myOptionsSens <- myOptionsODE

model <- odemodel(f, modelname = "preBotC_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "clamp", optionsOde = myOptionsODE, optionsSens = myOptionsSens)
g <- Y(observables, f = x, attach.input = TRUE,
       modelname = "preBotC_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)

carrier  <- unlist(lapply(red$blocks, `[[`, "carrierDomain"))
logRed   <- c(positive, names(carrier)[carrier == "positive"])

p.full <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~exp10(x)", x = intersect(.currentSymbols, positive)) |>
  P(modelname = "preBotC_pfull", compile = FALSE, outdir = .outdir)

p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  insert("x~exp10(x)", x = intersect(.currentSymbols, logRed)) |>
  P(modelname = "preBotC_pred", compile = FALSE, outdir = .outdir)

compile(x, g, p.full, p.red, output = .modelname, cores = .cores)

prd     <- g*x*p.full
prd.red <- g*x*p.red

toOuter <- function(v, logs) replace(v, names(v) %in% logs, log10(v[names(v) %in% logs]))
truthO  <- toOuter(truth[getParameters(prd)], positive)

times <- seq(0, 10000, by = 0.5)
pred0 <- as.data.frame(prd(times, truthO, deriv = FALSE))
ggplot(pred0, aes(time, value)) + geom_line() +
  facet_wrap(~ name, scales = "free_y", ncol = 1) + theme_dMod(base_size = 9)

# reduced truth from the invariants each outer parameter carries
meaning <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
truth.red <- sapply(getParameters(prd.red), function(p)
  if (p %in% names(meaning)) eval(parse(text = meaning[[p]]), as.list(truth))
  else truth[[p]])
truthO.red <- toOuter(truth.red, logRed)

obs   <- pred0$name == "y"
predR <- as.data.frame(prd.red(times, truthO.red, deriv = FALSE))
max(abs(predR$value[obs] - pred0$value[obs]))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Comparison: exp10() in the trafo of the detection
#
# The detection accepts parameters that enter only as exp10(k_l10) and reports
# in k_l10: the same three directions, the scaling now a translation. The
# reduction works in X = 10^k_l10, where they are rational again, and maps the
# chart back: the same q_k, with log() entries. Nothing is gained over reducing
# linearly and inserting exp10() afterwards, as above.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
posPars <- setdiff(positive, c("n", "h"))
lt <- as.eqnvec(setNames(paste0("exp10(", posPars, "_l10)"), posPars))

resL <- symmetryDetection(f, observables, trafo = lt, positive = c("n", "h"),
                          reconstruct = TRUE)
redL <- symmetryReduction(resL, positive = c("n", "h"))
redL
