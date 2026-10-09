# -------------------------------------------------------------------------#
# Macrospin ringdown and MTJ switching: what Kerr and TMR signals fix
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Landau-Lifshitz-Gilbert macrospin in spherical coordinates, observed through
# m_z = cos(theta). Part one is time-resolved polar MOKE on a film with
# perpendicular anisotropy: symmetryDetection() in an easy-axis and a tilted
# field, the easy-axis direction flown, both experiments fitted and profiled.
# Part two adds spin-transfer torque: the film as free layer of a
# perpendicular magnetic tunnel junction, switched by a current pulse and read
# through the TMR, with symmetryDetection(), symmetryReduction() and flows of
# the easy-axis directions.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Mon 28 Sep 2026
#
# [Info]
# Units: time in ns, fields in mT, gamma in rad/(ns mT), current density u in
# MA/cm^2, torques eta*u and xi*u in mT.
# The ODE models use includeTimeZero = FALSE, which multiple shooting needs.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryMacrospin"
# generated sources, objects and shared libraries go to tempdir()
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
#
# LLG in (th, ph) built from the energy gradient: uniaxial anisotropy field Hk
# along z and an applied field (Hx, 0, Hz).
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
# All parameters free, then the applied field known: along the easy axis and
# tilted, with the pump direction ph = 0.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
res.free <- symmetryDetection(f, observables, reconstruct = TRUE)
summary(res.free)

res.easy <- symmetryDetection(f, observables, trafo = eqnvec(Hx = "0"),
                              fixed = "Hz", reconstruct = TRUE)
summary(res.easy)

res.tilt <- symmetryDetection(f, observables, trafo = eqnvec(ph = "0"),
                              fixed = c("Hz", "Hx"), reconstruct = TRUE)
summary(res.tilt)

# the easy-axis direction that moves both alpha and gamma
dir.easy <- Filter(function(d) all(c("al", "gam") %in% d$support), res.easy$symmetries)[[1]]
dir.easy$generator


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# Two conditions, "easy" (Hx = 0) and "tilt". Fields known, ph0 = 0, the tilt
# th0 after the pump unknown per condition. mx and my are predicted only.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptions <- list(atol = 1e-8, rtol = 1e-8, maxsteps = 1e5)

model <- odemodel(f, modelname = "macrospin_ode", compile = FALSE, outdir = .outdir,
                  method = "tsit5", includeTimeZero = FALSE)
x <- Xs(model, options = myOptions)

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

# flow dz/deps = sgn*xi(z) along the easy-axis direction, in log10 coordinates
fl <- log10Transform(dir.easy$completeGenerator)
flow <- Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
                    deriv = FALSE, modelname = "macrospin_flow", compile = FALSE,
                    outdir = .outdir),
           condition = "flow", options = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))

compile(x, g, p, flow, output = .modelname, cores = .cores)

prd <- g*x*p


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Ringdown in both fields
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth  <- c(gam = log10(0.176), al = log10(0.05), Hk = log10(300), s = 0,
            th0_easy = 0.6, th0_tilt = 0.6)
times  <- seq(0, 1.5, by = 0.005)

plot(prd(times, truth, deriv = FALSE))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The easy-axis direction, flown
#
# The flow is run both ways from the truth; the observed y should stay put
# while the unobserved mx, my change.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mv   <- names(dir.easy$generator)
mvL  <- paste0(mv, "_l10")
eps  <- seq(0, 2, len = 30)

z0 <- setNames(truth[mv], mvL)
zf <- flow(eps, c(z0, sgn =  1), deriv = FALSE)$flow
zb <- flow(eps, c(z0, sgn = -1), deriv = FALSE)$flow

z <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])

# the orbit in (alpha, gamma), and the relaxation rate along it
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

# max deviation of y relative to its RMS, and the max change of mx
obs <- easy0$name == "y"
rms <- sqrt(mean(easy0$value[obs]^2))
dev.y  <- sapply(pred, function(d) max(abs(d$value[obs] - easy0$value[obs])) / rms)
dev.mx <- sapply(pred, function(d) max(abs(d$value[easy0$name == "mx"] - easy0$value[easy0$name == "mx"])))
dev.y
dev.mx

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
                      splitBy = "condition")
plot(prd(times, truth, deriv = FALSE), mydata)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and profile: easy axis alone, then both fields
#
# The fit in both fields uses multiple shooting with ph on an angle chart,
# whose gaps count modulo 2 pi. Profiles evaluate the single-shooting objective.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
stepControl <- list(stepsize = 1e-5, min = 1e-6, max = Inf,
                    atol = 1e-2, rtol = 1e-2, limit = 1e3)
optControl  <- list(rinit = 0.1, rmax = 10, iterlim = 2e3)
algoControl <- list(reoptimize = TRUE)

# starts uniform around the truth: gam, al and Hk within a factor 2,
# the Kerr scale within 25 %, the tilt within 0.3 rad
spread <- c(gam = 0.3, al = 0.3, Hk = 0.3, s = 0.1, th0_easy = 0.3, th0_tilt = 0.3)
lower  <- truth - spread[names(truth)]
upper  <- truth + spread[names(truth)]

# easy axis alone: th0_tilt does not enter and is held fixed
pars.easy <- names(truth)[names(truth) != "th0_tilt"]
msControl <- list(charts = c(ph = "angle"))
obj.easy  <- normL2(mydata["easy"], prd, times = times)

fit.easy  <- mstrust(obj.easy, 0*truth[pars.easy], samplefun = "runif",
                     min = lower[pars.easy], max = upper[pars.easy],
                     fixed = truth["th0_tilt"], rinit = 0.1, rmax = 10, fits = 200,
                     iterlim = 5e3, cores = .cores, name = "easy")
plotValues(as.parframe(fit.easy), tol = 0.1)
best.easy <- as.parvec(as.parframe(fit.easy))

prof.easy <- profile(obj.easy, best.easy, whichPar = names(best.easy),
                     fixed = truth["th0_tilt"], method = "integrate",
                     algoControl = algoControl, stepControl = stepControl,
                     optControl = optControl, limits = c(lower = -3, upper = 3),
                     cores = .cores)
plotProfile(prof.easy)

obj.both  <- normL2(mydata, prd, times = times, multipleShootingControl = msControl)
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MTJ switching: model
#
# Reference layer along z. The field-like torque enters as a field xi*u along
# z, the damping-like torque as eta*u*(m x p). The current u is a state, set to
# zero by an event at the end of the pulse.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
dEdth.stt <- "(Hk*sin(th)*cos(th) + (Hz + xi*u)*sin(th) - Hx*cos(th)*cos(ph))"

f.stt <- eqnvec(
  th = sprintf("-gam/(1 + al^2)*(al*%s + %s/sin(th) + eta*u*sin(th))", dEdth.stt, dEdph),
  ph = sprintf("gam/(1 + al^2)*(%s/sin(th) - al*%s/sin(th)^2 - al*eta*u)", dEdth.stt, dEdph),
  u  = "0")

.tpulse <- 20
ev.stt  <- addEvent(eventlist(), var = "u", time = .tpulse, value = 0, method = "replace")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MTJ switching: symmetries
#
# Field along the easy axis, three supercritical currents, tilt th0 unknown
# per condition; then the reduction of the directions found.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.col <- c("col_m6", "col_m8", "col_m10")
cond.col <- data.frame(u = c(-6, -8, -10), th = paste0("th0_", .col), row.names = .col)

res.stt <- symmetryDetection(f.stt, observables, trafo = eqnvec(Hx = "0", ph = "0"),
                             events = ev.stt, conditions = cond.col, fixed = "Hz",
                             reconstruct = TRUE)
summary(res.stt)

red.stt <- symmetryReduction(res.stt)
print(red.stt)

# adds two subcritical currents of both signs in a tilted field
.tilt <- c("tilt_m1", "tilt_p1")
cond.all <- rbind(
  data.frame(Hx = 0, cond.col),
  data.frame(Hx = 100, u = c(-1, 1), th = paste0("th0_", .tilt), row.names = .tilt))

res.stt.all <- symmetryDetection(f.stt, observables, trafo = eqnvec(ph = "0"),
                                 events = ev.stt, conditions = cond.all, fixed = "Hz",
                                 reconstruct = TRUE)
summary(res.stt.all)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MTJ switching: prediction chain and flows
#
# One condition per row of cond.all, Hz = 20 mT in all of them. One flow per
# direction of the easy-axis analysis, in linear coordinates.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
model.stt <- odemodel(f.stt, events = ev.stt, modelname = "mtj_ode", compile = FALSE,
                      outdir = .outdir, method = "tsit5")
x.stt <- Xs(model.stt, options = myOptions)

g.stt <- Y(c(observables, eqnvec(mx = "sin(th)*cos(ph)", my = "sin(th)*sin(ph)")),
           f = x.stt, modelname = "mtj_obs", compile = FALSE, outdir = .outdir)

innerpars.stt <- getParameters(g.stt, x.stt)
design <- data.frame(Hz = 20, cond.all[, c("Hx", "u")])

p.stt <- Reduce(`+`, lapply(rownames(design), function(cn) {
  eqnvec() |>
    define("x~x", x = innerpars.stt) |>
    define("x~y", x = c("Hz", "Hx", "u"), y = unlist(design[cn, ])) |>
    define("ph~0") |>
    define(paste0("th~th0_", cn)) |>
    insert("x~10^x", x = c("gam", "al", "Hk", "eta", "xi", "s")) |>
    P(condition = cn, modelname = paste0("mtj_p_", cn), compile = FALSE,
      outdir = .outdir)
}))

dirs.stt  <- res.stt$symmetries
flows.stt <- lapply(seq_along(dirs.stt), function(i) {
  gen <- dirs.stt[[i]]$completeGenerator
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", gen, ")"), names(gen))),
              deriv = FALSE, modelname = paste0("mtj_flow", i), compile = FALSE,
              outdir = .outdir),
     condition = paste0("flow", i),
     options = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
})

do.call(compile, c(list(x.stt, g.stt, p.stt), flows.stt,
                   list(output = "symmetryMTJ", cores = .cores)))

prd.stt <- g.stt*x.stt*p.stt


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MTJ switching: traces
#
# Gilbert damping 0.01, damping-like efficiency 1 mT per MA/cm^2, field-like
# ratio 0.2. The critical current uc is printed below.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth.lin <- c(gam = 0.176, al = 0.01, Hk = 300, eta = 1, xi = 0.2, s = 1, Hz = 20,
               setNames(rep(0.1, nrow(design)), paste0("th0_", rownames(design))))
.logpars   <- c("gam", "al", "Hk", "eta", "xi", "s")
truth.stt <- c(log10(truth.lin[.logpars]), truth.lin[grepl("^th0_", names(truth.lin))])
times.stt <- seq(0, 30, by = 0.01)

uc <- with(as.list(truth.lin), -al*(Hk + Hz)/(al*xi + eta))
uc

plot(prd.stt(times.stt, truth.stt, deriv = FALSE))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MTJ switching: the easy-axis directions, flown
#
# Each flow is run both ways from the truth. The TMR traces of the collinear
# conditions should stay put; uc and the moved parameters are tracked.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
eps <- seq(0, 1, len = 20)

orbits.stt <- lapply(seq_along(flows.stt), function(i) {
  fp <- setdiff(getParameters(flows.stt[[i]]), "sgn")
  zf <- flows.stt[[i]](eps, c(truth.lin[fp], sgn =  1), deriv = FALSE)[[paste0("flow", i)]]
  zb <- flows.stt[[i]](eps, c(truth.lin[fp], sgn = -1), deriv = FALSE)[[paste0("flow", i)]]
  z  <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
  z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])
  z
})

# range of uc along every orbit
uc.orbits <- lapply(orbits.stt, function(z) {
  mv <- setdiff(colnames(z), "time")
  range(apply(z[, mv, drop = FALSE], 1, function(v)
    with(as.list(replace(truth.lin, mv, v)), -al*(Hk + Hz)/(al*xi + eta))))
})
uc.orbits

# range of the moved parameters along every orbit
range.orbits <- lapply(orbits.stt, function(z) apply(z[, setdiff(colnames(z), "time"), drop = FALSE], 2, range))
range.orbits

# max deviation of y in the collinear conditions relative to its RMS
pred0.stt <- as.data.frame(prd.stt(times.stt, truth.stt, deriv = FALSE))
col0 <- pred0.stt[pred0.stt$condition %in% .col & pred0.stt$name == "y", ]
rms  <- sqrt(mean(col0$value^2))
dev.stt <- lapply(orbits.stt, function(z) {
  mv <- setdiff(colnames(z), "time")
  sapply(seq_len(nrow(z)), function(j) {
    lin <- replace(truth.lin, mv, z[j, mv])
    th  <- c(log10(lin[.logpars]), lin[grepl("^th0_", names(lin))])
    d   <- as.data.frame(prd.stt(times.stt, th, deriv = FALSE))
    d   <- d[d$condition %in% .col & d$name == "y", ]
    max(abs(d$value - col0$value)) / rms
  })
})
dev.stt
