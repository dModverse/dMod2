# -------------------------------------------------------------------------#
# Pre-Boetzinger pacemaker neuron: symmetries on a mixed-sign domain
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A literature model without a global positivity constraint: potentials,
# half-activation voltages and gate slopes take either sign; capacitance,
# conductances, time constants and gates are positive. Voltage is observed.
# symmetryDetection() finds the directions the voltage cannot resolve and
# symmetryReduction() removes them on the declared domain. Each direction is
# then flown both ways: parameters and currents move, the voltage does not.
# The last part runs the detection with exp10() in its trafo for comparison.
# Fitting this model is the subject of bench/multipleShooting_preBotC.R.
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
# Units mV, ms, pF, nS.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryPreBotC"
# generated sources, objects and shared libraries go to tempdir()
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
pinf  <- "1/(1 + exp((V - thp)/sp))"
minf  <- "1/(1 + exp((V - thm)/sm))"
ninf  <- "1/(1 + exp((V - thn)/sn))"
hinf  <- "1/(1 + exp((V - thh)/sh))"
tauN  <- "taun/cosh((V - thn)/(2*sn))"
tauH  <- "tauh/cosh((V - thh)/(2*sh))"

currents <- eqnvec(
  INaP = sprintf("gNaP*%s*h*(V - ENa)", pinf),
  INa  = sprintf("gNa*(%s)^3*(1 - n)*(V - ENa)", minf),
  IK   = "gK*n^4*(V - EK)",
  IL   = "gL*(V - EL)",
  Iton = "gton*(V - Esyn)")

f <- eqnvec(
  V = sprintf("-(%s)/C", paste(unclass(currents), collapse = " + ")),
  n = sprintf("(%s - n)/(%s)", ninf, tauN),
  h = sprintf("(%s - h)/(%s)", hinf, tauH))

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
# Detection with reconstructed generators, reduction on the declared positive
# set, and a second detection on the reduced chart, which must be identifiable.
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
# Positive coordinates and the outer parameters of the reduction certified
# positive go to log10, all others stay linear.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# spikes last about 2 ms; the check of the reduced chart needs tight solves
myOptions <- list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)

model <- odemodel(f, modelname = "preBotC_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "clamp", options = myOptions)
g <- Y(observables, f = x, attach.input = TRUE,
       modelname = "preBotC_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)

carrier <- unlist(lapply(red$blocks, `[[`, "carrierDomain"))
logRed  <- c(positive, names(carrier)[carrier == "positive"])

p.full <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~exp10(x)", x = intersect(.currentSymbols, positive)) |>
  P(modelname = "preBotC_pfull", compile = FALSE, outdir = .outdir)

p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  insert("x~exp10(x)", x = intersect(.currentSymbols, logRed)) |>
  P(modelname = "preBotC_pred", compile = FALSE, outdir = .outdir)

# observable plus the currents, to show what a direction moves
gI <- Y(c(observables, currents), f = x, modelname = "preBotC_currents",
        compile = FALSE, outdir = .outdir)

# Flow dz/deps = sgn*xi(z) per direction: positive coordinates on log10, the
# others linear. The sign runs the same compiled model both ways.
flowEq <- lapply(res$symmetries, function(s) {
  gen <- as.eqnvec(s$completeGenerator)
  pos <- intersect(names(gen), positive)
  lin <- setdiff(names(gen), pos)
  eqPos <- character(0)
  eqLin <- character(0)
  if (length(pos)) eqPos <- unclass(log10Transform(gen[pos]))
  if (length(lin)) eqLin <- setNames(unclass(gen[lin]), lin)
  if (length(lin) && length(pos))
    eqLin[] <- cOde::replaceSymbols(pos, sprintf("(10^(%s_l10))", pos), eqLin)
  eq <- c(eqPos, eqLin)
  as.eqnvec(setNames(paste0("sgn*(", eq, ")"), names(eq)))
})

flow <- lapply(seq_along(flowEq), function(i)
  Xf(odemodel(flowEq[[i]], deriv = FALSE, modelname = paste0("preBotC_flow", i),
              compile = FALSE, outdir = .outdir),
     condition = "flow", options = myOptions))

do.call(compile, c(list(x, g, gI, p.full, p.red), flow,
                   list(output = .modelname, cores = .cores)))

prd     <- g*x*p.full
prd.red <- g*x*p.red

truthO <- truth[getParameters(prd)]
isLog  <- names(truthO) %in% positive
truthO[isLog] <- log10(truthO[isLog])

times <- seq(0, 10000, by = 0.5)
pred0 <- as.data.frame(prd(times, truthO, deriv = FALSE))
ggplot(pred0, aes(time, value)) + geom_line() +
  facet_wrap(~ name, scales = "free_y", ncol = 1) + theme_dMod(base_size = 9)

# reduced truth: each outer parameter of the reduction evaluated as its invariant
meaning   <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
parsRed   <- getParameters(prd.red)
truth.red <- setNames(truth[parsRed], parsRed)
fromMeaning <- intersect(parsRed, names(meaning))
truth.red[fromMeaning] <- sapply(meaning[fromMeaning], function(m)
  eval(parse(text = m), as.list(truth)))

truthO.red <- truth.red
isLogRed   <- names(truthO.red) %in% logRed
truthO.red[isLogRed] <- log10(truthO.red[isLogRed])

obs   <- pred0$name == "y"
predR <- as.data.frame(prd.red(times, truthO.red, deriv = FALSE))
max(abs(predR$value[obs] - pred0$value[obs]))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Every direction, flown
#
# From the truth both ways, over a range that moves the parameters by up to
# half a decade or 5 mV to first order. At each orbit point the model is solved
# over one burst; the verdict is the largest voltage deviation relative to its RMS.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
prdI   <- gI*x*p.full
tw     <- seq(2000, 3000, by = 0.5)
predI0 <- as.data.frame(prdI(tw, truthO, deriv = FALSE))
obsI   <- predI0$name == "y"
rmsV   <- sqrt(mean(predI0$value[obsI]^2))

# 21 points per orbit, eps symmetric around the truth
orbits <- list()
for (i in seq_along(flow)) {
  mv   <- names(res$symmetries[[i]]$completeGenerator)
  pos  <- intersect(mv, positive)
  lin  <- setdiff(mv, pos)
  posL <- sprintf("%s_l10", pos)
  z0   <- c(setNames(truthO[pos], posL), truth[lin], truth[setdiff(coords, mv)])

  # first-order speed at the truth, per decade (log10) or per 10 mV (linear)
  spd <- abs(sapply(names(flowEq[[i]]), function(k)
    eval(parse(text = flowEq[[i]][[k]]), as.list(c(z0, sgn = 1)))))
  spd <- spd / ifelse(names(spd) %in% lin, 10, 1)
  eps <- seq(0, 0.5 / max(spd), length.out = 11)

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
  orbits[[i]] <- list(eps = z[, "time"], z = z, theta = theta, moved = mv)
}

# voltage and currents along each orbit, stacked over eps
orbitPred <- lapply(orbits, function(o)
  do.call(rbind, Map(function(th, e) transform(as.data.frame(prdI(tw, th, deriv = FALSE)), eps = e),
                     o$theta, o$eps)))

verdict <- data.frame(
  direction  = paste0("X", seq_along(orbits)),
  type       = sapply(res$symmetries, `[[`, "type"),
  moves      = sapply(orbits, function(o) paste(o$moved, collapse = ",")),
  epsMin     = sapply(orbits, function(o) min(o$eps)),
  epsMax     = sapply(orbits, function(o) max(o$eps)),
  maxMove    = sapply(orbits, function(o) max(abs(sapply(o$theta, function(th) th - truthO)))),
  devV       = sapply(orbitPred, function(q) {
    qy <- q[q$name == "y", ]
    max(abs(qy$value - predI0$value[obsI][match(qy$time, predI0$time[obsI])])) / rmsV
  }))
verdict

orbitPred <- lapply(orbitPred, function(q) {
  q$name <- factor(ifelse(q$name == "y", "V [observed]", as.character(q$name)),
                   c("V [observed]", names(currents)))
  q
})

# the currents move, the voltage does not
ggplot(orbitPred[[1]], aes(time, value, group = eps, colour = eps)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = sprintf("pre-Boetzinger, X1 (%s)", res$symmetries[[1]]$type),
       x = "time [ms]", y = NULL) +
  theme_dMod(base_size = 9)

ggplot(orbitPred[[2]], aes(time, value, group = eps, colour = eps)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = sprintf("pre-Boetzinger, X2 (%s)", res$symmetries[[2]]$type),
       x = "time [ms]", y = NULL) +
  theme_dMod(base_size = 9)

ggplot(orbitPred[[3]], aes(time, value, group = eps, colour = eps)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = sprintf("pre-Boetzinger, X3 (%s)", res$symmetries[[3]]$type),
       x = "time [ms]", y = NULL) +
  theme_dMod(base_size = 9)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Comparison: exp10() in the trafo of the detection
#
# Positive parameters enter as exp10(k_l10); detection and reduction then work
# in k_l10. Compare redL$trafo with red$trafo above followed by exp10().
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
posPars <- setdiff(positive, c("n", "h"))
lt <- as.eqnvec(setNames(paste0("exp10(", posPars, "_l10)"), posPars))

resL <- symmetryDetection(f, observables, trafo = lt, positive = c("n", "h"),
                          reconstruct = TRUE)
redL <- symmetryReduction(resL, positive = c("n", "h"))
redL
