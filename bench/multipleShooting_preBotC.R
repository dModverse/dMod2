# -------------------------------------------------------------------------#
# Pre-Boetzinger pacemaker neuron: a burst of spikes fitted by multiple shooting
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The model of bench/symmetry_preBotC.R, in the chart its symmetry reduction
# certifies, fitted to one simulated burst of its voltage: thirteen spikes of
# about 2 ms in 700 ms, sampled at 10 kHz with 0.25 mV noise, as a patch-clamp
# recording would give it. A change of a parameter by a per cent moves the late
# spikes of the burst by more than their width, so single shooting finds the
# optimum only from next to it. Multiple shooting with nodes = "transitions"
# puts a node just before every spike: each spike starts where the data have
# it, the slow stretch behind it depends on the parameters smoothly, and a
# wrong timing shows as a gap.
#
# 1. At the true parameters both reach the optimum.
# 2. From ten starts at 0.02, 0.05 and 0.1 decades from the truth (positive
#    parameters in log10; potentials, half-activations and slopes by a fifth of
#    that relatively), single shooting and multiple shooting. A start is a hit
#    when the single-shooting objective at its end lies within 1 of the
#    optimum.
# 3. A start single shooting misses and multiple shooting hits, against the
#    data.
#
# On 26 Sep 2026, half an hour on twenty cores: single shooting hit from none
# of the 30 starts, multiple shooting from 4, 4 and 1 of 10 at 0.02, 0.05 and
# 0.1 decades. Its misses end where single shooting ends, with spikes missing
# or added.
#
# With MS_RESULTS set, the numbers and the trajectories of the figure are
# written there; dev/optimisation/prebotc.rds, which the Optimisation vignette
# reads, is that file.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Sat 26 Sep 2026
#
# [Info]
# Model 1 of Butera RJ, Rinzel J, Smith JC (1999). J Neurophysiol 82, 382-397,
# as in bench/symmetry_preBotC.R. The ODE model is built with
# includeTimeZero = FALSE, which multiple shooting needs: a segment starts at
# its own node, not at zero. The solver runs at 1e-10: at 1e-8 the objective
# scatters by more than a step predicts far from the optimum. Needs cppDE
# 0.10.3, whose first step resolves a late start under that tolerance.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "msPreBotC"
.outdir    <- Sys.getenv("MS_OUTDIR", file.path(tempdir(), .modelname))
.cores     <- detectFreeCores()
.fits      <- as.integer(Sys.getenv("MS_FITS", "10"))

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model and reduced chart, as in bench/symmetry_preBotC.R
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

truth <- c(
  V = -60, n = 0.01, h = 0.6,
  C = 21, gNaP = 2.8, gNa = 28, gK = 11.2, gL = 2.8, gton = 0.3,
  ENa = 50, EK = -85, EL = -65, Esyn = 0,
  thm = -34, sm = -5, thn = -29, sn = -4, thp = -40, sp = -6, thh = -48, sh = 6,
  taun = 10, tauh = 10000)

res <- symmetryDetection(f, observables, positive = positive, reconstruct = TRUE)
red <- symmetryReduction(res, positive = positive)
red$trafo

# 1e-10: far from the optimum, where misplaced spikes make the residuals
# large, a looser solve scatters the objective by more than the descent a step
# predicts. The cap of 10^6 steps ends a trial point in a regime the data never
# show, which trust() then rejects, instead of integrating it for minutes.
opts <- list(atol = 1e-10, rtol = 1e-10, maxsteps = 1e6)

model <- odemodel(f, modelname = "msPreBotC_ode", compile = FALSE,
                  outdir = .outdir, includeTimeZero = FALSE)
x <- Xs(model, condition = "clamp", options = opts)
g <- Y(observables, f = x, modelname = "msPreBotC_obs", compile = FALSE,
       outdir = .outdir)

innerpars <- getParameters(g, x)
carrier <- unlist(lapply(red$blocks, `[[`, "carrierDomain"))
logRed  <- c(positive, names(carrier)[carrier == "positive"])

p <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  insert("x~exp10(x)", x = intersect(.currentSymbols, logRed)) |>
  P(modelname = "msPreBotC_p", compile = FALSE, outdir = .outdir)

compile(x, g, p, output = .modelname, cores = .cores)
prd <- g * x * p

toOuter <- function(v, logs) replace(v, names(v) %in% logs, log10(v[names(v) %in% logs]))
meaning <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
truth.red <- sapply(getParameters(prd), function(q)
  if (q %in% names(meaning)) eval(parse(text = meaning[[q]]), as.list(truth))
  else truth[[q]])
truthO <- toOuter(truth.red, logRed)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Data: one burst, the voltage at 10 kHz for 700 ms, 0.25 mV noise
#
# The window starts at 2100 ms of the trajectory from the published initial
# state, shortly before the second burst; its state there is the initial state
# of the fit, n and h in log10 like every positive quantity.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
full <- (x * p)(seq(0, 2100, by = 1), truthO, deriv = FALSE)[[1]]
x0 <- full[nrow(full), c("V", "n", "h")]
truthB <- truthO
truthB[c("V", "n", "h")] <- c(x0[["V"]], log10(x0[["n"]]), log10(x0[["h"]]))

times <- seq(0, 700, by = 0.1)
set.seed(5)
Vtrue <- prd(times, truthB, deriv = FALSE)[[1]][, "y"]
data <- as.datalist(data.frame(name = "y", time = times,
                               value = Vtrue + rnorm(length(times), 0, 0.25),
                               sigma = 0.25, condition = "clamp"))
spikes <- function(v) sum(diff(v > -20) == 1)
spikes(Vtrue)

charts <- c(n = "log10", h = "log10")
obj.ss <- normL2(data, prd)
obj.ms <- normL2(data, prd, multipleShootingControl = list(nodes = "transitions",
                                                           charts = charts))
ss <- function(th) obj.ss(th, deriv = FALSE)$value

# the nodes the data give: one just before every spike
attr(dMod2:::.shootObjOf(obj.ms), "spec")$tau$clamp

fitOne <- function(obj, s, iterlim, cores = .cores)
  suppressWarnings(trust(obj, s, rinit = 0.1, rmax = 10, iterlim = iterlim,
                         cores = cores))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. At the true parameters
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
atTruth <- list(
  "single shooting"   = fitOne(obj.ss, truthB, 500),
  "multiple shooting" = fitOne(obj.ms, truthB, 300))
optimum <- min(vapply(atTruth, function(r) ss(r$argument), 0.0))
truth.table <- data.frame(
  run    = names(atTruth),
  ss     = vapply(atTruth, function(r) ss(r$argument), 0.0),
  gaps   = vapply(atTruth, function(r)
    if (is.null(r$multipleShooting)) NA_real_ else r$multipleShooting$violation, 0.0),
  row.names = NULL)
truth.table
c(truth = ss(truthB), optimum = optimum)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2. From perturbed starts
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
pars <- names(truthB)
fixedInit <- intersect(c("V", "n", "h"), pars)
free <- setdiff(pars, fixedInit)
lg <- free[free %in% logRed]; ln <- setdiff(free, lg)
scales <- c(0.02, 0.05, 0.1)
starts <- lapply(setNames(scales, scales), function(sc) {
  set.seed(round(1000 * sc))
  lapply(seq_len(.fits), function(i) {
    s <- truthB
    s[lg] <- s[lg] + rnorm(length(lg), 0, sc)
    s[ln] <- s[ln] * (1 + rnorm(length(ln), 0, sc / 5))
    s
  })
})

jobs <- expand.grid(method = c("single shooting", "multiple shooting"),
                    scale = scales, start = seq_len(.fits), stringsAsFactors = FALSE)
# one start per core, each fit on one
out <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), mc.cores = .cores,
                                         mc.preschedule = FALSE, function(j) {
  jb <- jobs[j, ]
  s <- starts[[as.character(jb$scale)]][[jb$start]]
  t0 <- Sys.time()
  r <- try(if (jb$method == "single shooting") fitOne(obj.ss, s, 500, 1L)
           else fitOne(obj.ms, s, 300, 1L), silent = TRUE)
  sec <- as.numeric(Sys.time() - t0, units = "secs")
  th <- if (inherits(r, "try-error")) NULL else r$argument
  data.frame(jb, ss = if (is.null(th)) NA_real_ else ss(th),
             spikes = if (is.null(th)) NA_integer_
                      else spikes(prd(times, th, deriv = FALSE)[[1]][, "y"]),
             sec = sec, theta = I(list(th)))
}))
out$hit <- !is.na(out$ss) & out$ss < optimum + 1
out$fits <- 1L
start.table <- aggregate(cbind(hits = hit, fits = fits) ~ method + scale, out, sum)
start.table$medianSec <- aggregate(sec ~ method + scale, out, function(z) round(median(z)))$sec
start.table


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. One start, both methods, against the data
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
pick <- merge(subset(out, method == "multiple shooting" & hit, c(scale, start)),
              subset(out, method == "single shooting" & !hit, c(scale, start)))
if (!nrow(pick))   # no such start: the best multiple-shooting fit
  pick <- subset(out, method == "multiple shooting")[
    which.min(subset(out, method == "multiple shooting")$ss), c("scale", "start")]
pick <- pick[order(-pick$scale), ][1, ]
thOf <- function(m) out$theta[[which(out$method == m & out$scale == pick$scale &
                                     out$start == pick$start)]]
tf <- seq(0, 700, by = 0.05)
panels <- c("start", "single shooting", "multiple shooting")
traj <- do.call(rbind, Map(function(th, lab) data.frame(
  time = tf, value = prd(tf, th, deriv = FALSE)[[1]][, "y"], run = lab),
  list(starts[[as.character(pick$scale)]][[pick$start]], thOf("single shooting"),
       thOf("multiple shooting")), panels))
traj$run <- factor(traj$run, panels)
shown <- data$clamp[seq(1, nrow(data$clamp), by = 5), ]

ggplot(traj, aes(time, value)) +
  geom_point(data = shown, size = 0.15, colour = "grey60") +
  geom_line(linewidth = 0.3, colour = "#0072B2") +
  facet_wrap(~run, ncol = 1) +
  labs(x = "time [ms]", y = "V [mV]",
       title = sprintf("start %d at %.2f decades", pick$start, pick$scale)) +
  theme_dMod(base_size = 9)

# The curves are stored to 0.01 mV, far below the noise, which keeps the file
# small.
if (nzchar(Sys.getenv("MS_RESULTS")))
  saveRDS(list(truth = truth.table, optimum = optimum, starts = start.table,
               perStart = out[, setdiff(names(out), c("theta", "fits"))], pick = pick,
               trajectories = transform(traj, value = round(value, 2)),
               data = transform(shown, value = round(value, 2)), spikes = spikes(Vtrue),
               nodes = attr(dMod2:::.shootObjOf(obj.ms), "spec")$tau$clamp),
          Sys.getenv("MS_RESULTS"))
