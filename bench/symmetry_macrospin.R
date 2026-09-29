# -------------------------------------------------------------------------#
# Macrospin ringdown: what a polar Kerr signal fixes of damping and gamma
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Time-resolved MOKE on a thin film with perpendicular anisotropy: a pump
# pulse tilts the magnetization, which then precesses back to equilibrium
# under the Landau-Lifshitz-Gilbert equation. The polar Kerr signal reads
# m_z = cos(theta). symmetryDetection() works on the LLG in spherical
# coordinates, sin() and cos() of the states included. It finds that a field
# along the easy axis hides the precession from the Kerr signal, so that
# damping alpha and gyromagnetic ratio gamma cannot be told apart, and that a
# tilted field separates them. The easy-axis direction is flown, then both
# experiments are fitted and profiled.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Mon 28 Sep 2026
#
# [Info]
# Needs dMod2 installed from devel-symmetry. Units: time in ns, fields in mT,
# gamma in rad/(ns mT), so that gamma = 0.176 is the free-electron value.
# phi0 = 0 fixes the in-plane direction the pump tilts into.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryMacrospin"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
#
# Energy density per Ms in field units,
#   e = Hk/2 sin(th)^2 - Hz cos(th) - Hx sin(th) cos(ph),
# with uniaxial anisotropy field Hk along z and the applied field (Hx, 0, Hz).
# LLG in (th, ph):
#   th' = -gam/(1 + al^2) (al de/dth + de/dph / sin(th))
#   ph' =  gam/(1 + al^2) (de/dth / sin(th) - al de/dph / sin(th)^2)
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
dEdth <- "(Hk*sin(th)*cos(th) + Hz*sin(th) - Hx*cos(th)*cos(ph))"
dEdph <- "(Hx*sin(th)*sin(ph))"

f <- eqnvec(
  th = sprintf("-gam/(1 + al^2)*(al*%s + %s/sin(th))", dEdth, dEdph),
  ph = sprintf("gam/(1 + al^2)*(%s/sin(th) - al*%s/sin(th)^2)", dEdth, dEdph))

observables <- eqnvec(y = "s*cos(th)")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Symmetries
#
# Everything free: only frequencies gam*H are seen, so gamma scales against
# every field. The applied field is the experimental control and therefore
# known; Hk stays unknown.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
res.free <- symmetryDetection(f, observables, reconstruct = TRUE)
summary(res.free)

# Field along the easy axis: ph drops out of th', the Kerr signal sees only
# the relaxation rate gam*al/(1 + al^2)*(Hk*cos(th) + Hz) and no precession.
# Two directions: the azimuth, which nothing observes, and the trade between
# alpha and gamma along that rate.
res.easy <- symmetryDetection(f, observables, trafo = eqnvec(Hx = "0"),
                              fixed = "Hz", reconstruct = TRUE)
summary(res.easy)

# Tilted field: the precession projects onto m_z and alpha, gamma and Hk
# separate.
res.tilt <- symmetryDetection(f, observables, trafo = eqnvec(ph = "0"),
                              fixed = c("Hz", "Hx"), reconstruct = TRUE)
summary(res.tilt)

# the rate the easy-axis experiment fixes is invariant along its direction
dir.easy <- Filter(function(d) all(c("al", "gam") %in% d$support), res.easy$symmetries)[[1]]
dir.easy$generator


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# Two conditions, one per field: "easy" (Hx = 0) and "tilt". The tilt angle
# th0 after the pump is unknown per condition, ph0 = 0, fields known. The
# in-plane components mx, my are predicted but never observed.
#
# The LLG is a lightly damped oscillator, eigenvalues gam*H*(-al +- i), and not
# stiff, so the explicit Tsit5 fits it; BDF is slower at every tolerance. A
# start that makes the precession orders of magnitude too fast fails at
# maxsteps within milliseconds instead of integrating for minutes.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptionsODE  <- list(atol = 1e-8, rtol = 1e-8, maxsteps = 1e5)
myOptionsSens <- myOptionsODE

model <- odemodel(f, modelname = "macrospin_ode", compile = FALSE, outdir = .outdir,
                  method = "tsit5")
x <- Xs(model, optionsOde = myOptionsODE, optionsSens = myOptionsSens)

g <- Y(c(observables, eqnvec(mx = "sin(th)*cos(ph)", my = "sin(th)*sin(ph)")),
       f = x, modelname = "macrospin_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)
fields    <- data.frame(Hz = c(100, 30), Hx = c(0, 100), row.names = c("easy", "tilt"))

# positive physical parameters on log10, the tilt angle linear
p <- Reduce(`+`, lapply(rownames(fields), function(cn) {
  eqnvec() |>
    define("x~x", x = innerpars) |>
    define("x~y", x = c("Hz", "Hx"), y = unlist(fields[cn, ])) |>
    define("ph~0") |>
    define(paste0("th~th0_", cn)) |>
    insert("x~10^x", x = c("gam", "al", "Hk", "s")) |>
    P(condition = cn, modelname = paste0("macrospin_p_", cn), compile = FALSE,
      outdir = .outdir)
}))

# One flow along the easy-axis direction, dz/deps = sgn*xi(z), in the log10
# coordinates of the moved parameters.
fl <- log10Transform(dir.easy$completeGenerator)
flow <- Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
                    deriv = FALSE, modelname = "macrospin_flow", compile = FALSE,
                    outdir = .outdir),
           condition = "flow", optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))

# One call compiles ODE, observation, trafos and the flow together.
compile(x, g, p, flow, output = .modelname, cores = .cores)

prd <- g*x*p


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Ringdown in both fields
#
# On the easy axis m_z only relaxes; in the tilted field it rings.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth  <- c(gam = log10(0.176), al = log10(0.05), Hk = log10(300), s = 0,
            th0_easy = 0.6, th0_tilt = 0.6)
times  <- seq(0, 1.5, by = 0.005)

plot(prd(times, truth, deriv = FALSE))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The easy-axis direction, flown
#
# Along the flow the Kerr signal stays put while the precession in the plane,
# which the Kerr signal does not see, runs faster or slower.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mv   <- names(dir.easy$generator)                   # al, gam
mvL  <- paste0(mv, "_l10")
eps  <- seq(0, 2, len = 30)

z0 <- setNames(truth[mv], mvL)
zf <- flow(eps, c(z0, sgn =  1), deriv = FALSE)$flow
zb <- flow(eps, c(z0, sgn = -1), deriv = FALSE)$flow

z <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])

# the orbit in (alpha, gamma): the curve gam*al/(1 + al^2) = const
orbit <- data.frame(eps = z[, "time"], al = 10^z[, "al_l10"], gam = 10^z[, "gam_l10"])
ggplot(orbit, aes(al, gam, colour = eps)) + geom_path() + geom_point() +
  scale_x_log10() + scale_y_log10() +
  scale_color_dMod_div(name = expression(epsilon)) + theme_dMod(base_size = 9)
range(orbit$gam*orbit$al/(1 + orbit$al^2))

pred0 <- as.data.frame(prd(times, truth, deriv = FALSE))
easy0 <- pred0[pred0$condition == "easy", ]
pred  <- lapply(seq_len(nrow(z)), function(j) {
  th <- replace(truth, mv, z[j, mvL])
  d  <- as.data.frame(prd(times, th, deriv = FALSE))
  transform(d[d$condition == "easy", ], eps = z[j, "time"])
})

# verdict: max deviation of y relative to its RMS, and how far mx moved
obs <- easy0$name == "y"
rms <- sqrt(mean(easy0$value[obs]^2))
sapply(pred, function(d) max(abs(d$value[obs] - easy0$value[obs])) / rms)
sapply(pred, function(d) max(abs(d$value[easy0$name == "mx"] - easy0$value[easy0$name == "mx"])))

long <- do.call(rbind, pred)
long$name <- factor(long$name, c("y", "mx", "my"), c("y [observed]", "mx", "my"))
ggplot(long, aes(time, value, group = eps, colour = eps)) +
  geom_line() + facet_wrap(~ name, ncol = 1, scales = "free_y") +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = "easy axis, alpha-gamma direction", x = "time [ns]", y = NULL) +
  theme_dMod(base_size = 9)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Simulate data
#
# Kerr signal every 10 ps with noise sd 0.005, in both fields.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.sigma <- 0.005
set.seed(1111)
sim <- as.data.frame(prd(seq(0, 1.5, by = 0.01), truth, deriv = FALSE))
sim <- sim[sim$name == "y", ]
sim$value <- sim$value + rnorm(nrow(sim), 0, .sigma)
sim$sigma <- .sigma

mydata <- as.datalist(sim[, c("name", "time", "value", "sigma", "condition")],
                      split.by = "condition")
plot(prd(times, truth, deriv = FALSE), mydata)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and profile: easy axis alone, then both fields
#
# Easy axis alone, th0_tilt does not enter and is left out. Its profiles of
# al and gam are flat along the direction above; with the tilted field they
# close, narrowly, since the precession frequency is fitted to high precision.
# The ringing makes the joint fit multimodal, so mstrust is needed.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
stepControl <- list(stepsize = 1e-5, min = 1e-6, max = Inf,
                    atol = 1e-2, rtol = 1e-2, limit = 1e3)
optControl  <- list(rinit = 0.1, rmax = 10, iterlim = 2e3)
algoControl <- list(reoptimize = TRUE)

# starts drawn uniformly from what is known before the experiment: g between 1
# and 3.4, alpha between 1e-3 and 0.3, Hk between 10 mT and 1 T, the Kerr scale
# within a factor 2, the tilt below 70 degrees
lower <- c(gam = log10(0.088), al = -3, Hk = 1, s = -0.3, th0_easy = 0.1, th0_tilt = 0.1)
upper <- c(gam = log10(0.3), al = log10(0.3), Hk = 3, s = 0.3, th0_easy = 1.2,
           th0_tilt = 1.2)

pars.easy <- names(truth)[names(truth) != "th0_tilt"]
obj.easy  <- normL2(mydata["easy"], prd, times = times)
fit.easy  <- mstrust(obj.easy, 0*truth[pars.easy], samplefun = "runif",
                     min = lower[pars.easy], max = upper[pars.easy],
                     fixed = truth["th0_tilt"], rinit = 0.1, rmax = 10, fits = 100,
                     iterlim = 5e3, cores = .cores, name = "easy")
plotValues(as.parframe(fit.easy), tol = 0.1)
best.easy <- as.parvec(as.parframe(fit.easy))

prof.easy <- profile(obj.easy, best.easy, whichPar = names(best.easy),
                     fixed = truth["th0_tilt"], method = "integrate",
                     algoControl = algoControl, stepControl = stepControl,
                     optControl = optControl, limits = c(lower = -3, upper = 3),
                     cores = .cores)
plotProfile(prof.easy)

obj.both  <- normL2(mydata, prd, times = times)
fit.both  <- mstrust(obj.both, 0*truth, samplefun = "runif", min = lower[names(truth)],
                     max = upper[names(truth)], rinit = 0.1, rmax = 10, fits = 100,
                     iterlim = 1e3, cores = .cores, name = "both")
plotValues(as.parframe(fit.both), tol = 0.1)
best.both <- as.parvec(as.parframe(fit.both))
plot(prd(times, best.both, deriv = FALSE), mydata)

prof.both <- profile(obj.both, best.both, whichPar = names(best.both),
                     method = "integrate", algoControl = algoControl,
                     stepControl = stepControl, optControl = optControl,
                     limits = c(lower = -2, upper = 2), cores = .cores)
plotProfile(prof.both)
plotProfile(list(easy = prof.easy, both = prof.both))
