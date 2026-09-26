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
# directions, symmetryReduction() removes all three on the declared domain.
# Every direction is then flown both ways, as in bench/symmetry_AB.R: the
# parameters move along the orbit, the currents they set move with them, and
# the voltage does not. The last part compares with exp10() in the trafo handed
# to the detection. A spiking model is not fitted here; see
# bench/multipleShooting_preBotC.R.
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

# The currents of the model, to watch what a direction moves; the voltage is
# the observable.
currents <- eqnvec(
  INaP = sprintf("gNaP*%s*h*(V - ENa)", xinf("p")),
  INa  = sprintf("gNa*(%s)^3*(1 - n)*(V - ENa)", xinf("m")),
  IK   = "gK*n^4*(V - EK)",
  IL   = "gL*(V - EL)",
  Iton = "gton*(V - Esyn)")
gI <- Y(c(observables, currents), f = x, modelname = "preBotC_currents",
        compile = FALSE, outdir = .outdir)

# One flow per direction, dz/deps = sgn*xi(z), in the chart of the prediction:
# the positive coordinates a direction moves on log10, the others, which take
# either sign, linear. The sign runs the same compiled model both ways.
flowOf <- function(gen) {
  pos <- intersect(names(gen), positive)
  lin <- setdiff(names(gen), pos)
  eq <- if (length(pos)) unclass(log10Transform(gen[pos])) else character(0)
  if (length(lin)) {
    el <- unclass(gen[lin])
    if (length(pos))
      el <- cOde::replaceSymbols(pos, sprintf("(10^(%s_l10))", pos), el)
    eq <- c(eq, setNames(el, lin))
  }
  as.eqnvec(setNames(paste0("sgn*(", eq, ")"), names(eq)))
}
flow <- lapply(seq_along(res$symmetries), function(i)
  Xf(odemodel(flowOf(as.eqnvec(res$symmetries[[i]]$completeGenerator)), deriv = FALSE,
              modelname = paste0("preBotC_flow", i), compile = FALSE, outdir = .outdir),
     condition = "flow", optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)))

do.call(compile, c(list(x, g, gI, p.full, p.red), flow,
                   list(output = .modelname, cores = .cores)))

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
# Every direction, flown
#
# From the truth both ways, over a range that moves the parameters by up to
# half a decade or 5 mV to first order. At every point of the orbit the model
# is solved over one burst: the verdict is the largest deviation of the
# voltage relative to its RMS, the plot shows voltage and currents coloured by
# eps. X1 scales the capacitance and every conductance together, so every
# current scales and V stays; X2 and X3 trade the leak against the tonic drive.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
prdI  <- gI*x*p.full
tw    <- seq(2000, 3000, by = 0.5)
predI0 <- as.data.frame(prdI(tw, truthO, deriv = FALSE))
obsI  <- predI0$name == "y"
rmsV  <- sqrt(mean(predI0$value[obsI]^2))

orbit <- function(i, n = 21) {
  mv   <- names(res$symmetries[[i]]$completeGenerator)
  pos  <- intersect(mv, positive); lin <- setdiff(mv, pos)
  posL <- sprintf("%s_l10", pos)
  z0   <- c(setNames(truthO[pos], posL), truth[lin], truth[setdiff(coords, mv)])
  # first-order speed of the flow at the truth, per decade or per 5 mV
  gen  <- flowOf(as.eqnvec(res$symmetries[[i]]$completeGenerator))
  spd  <- abs(sapply(names(gen), function(k)
    eval(parse(text = gsub("sgn", "1", gen[[k]])), as.list(c(z0, sgn = 1)))))
  spd  <- spd / ifelse(names(spd) %in% lin, 10, 1)
  emax <- 0.5 / max(spd)
  eps  <- seq(0, emax, length.out = (n + 1) / 2)
  zf <- flow[[i]](eps, c(z0, sgn =  1), deriv = FALSE)$flow
  zb <- flow[[i]](eps, c(z0, sgn = -1), deriv = FALSE)$flow
  z  <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
  z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])
  theta <- lapply(seq_len(nrow(z)), function(j) {
    th <- truthO
    th[pos] <- z[j, posL]
    th[lin] <- z[j, lin]
    th
  })
  list(eps = z[, "time"], z = z, theta = theta, moved = mv)
}

orbits <- lapply(seq_along(flow), orbit)
for (i in seq_along(orbits)) {
  o <- orbits[[i]]
  d <- sapply(o$theta, function(th) {
    q <- as.data.frame(prdI(tw, th, deriv = FALSE))
    max(abs(q$value[q$name == "y"] - predI0$value[obsI])) / rmsV
  })
  moved <- max(abs(sapply(o$theta, function(th) max(abs(th - truthO)))))
  cat(sprintf("X%d  %-8s on %-26s  eps [%+.3f, %+.3f]  moved %.2f  V deviates %.2e\n", i,
              res$symmetries[[i]]$type, paste(o$moved, collapse = ","),
              min(o$eps), max(o$eps), moved, max(d)))
}

# the currents move, the voltage does not
flowPlot <- function(i) {
  o <- orbits[[i]]
  long <- do.call(rbind, Map(function(th, e) transform(as.data.frame(prdI(tw, th, deriv = FALSE)),
                                                      eps = e), o$theta, o$eps))
  long$name <- factor(ifelse(long$name == "y", "V [observed]", as.character(long$name)),
                      c("V [observed]", names(currents)))
  ggplot(long, aes(time, value, group = eps, colour = eps)) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ name, scales = "free_y", ncol = 2) +
    scale_color_dMod_div(name = expression(epsilon)) +
    labs(title = sprintf("pre-Boetzinger, X%d (%s)", i, res$symmetries[[i]]$type),
         x = "time [ms]", y = NULL) +
    theme_dMod(base_size = 9)
}
flowPlot(1)
flowPlot(2)
flowPlot(3)


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
