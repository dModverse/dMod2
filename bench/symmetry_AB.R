# -------------------------------------------------------------------------#
# A <-> B: what a symmetry does to a profile likelihood
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The smallest model with a non-scaling symmetry: A converts to B and back,
# only s*B is observed. symmetryDetection() finds two directions. Both are
# flown both ways: the states move, the observable does not. Then they are
# profiled: s is fixed as a gauge, the curved direction stays flat until its
# orbit leaves the model's domain, and symmetryReduction() closes that profile.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 23 Sep 2026
#
# [Info]
# Needs dMod2 installed from devel-symmetry. The observable is y = log2(s*B),
# the log of the rational observable the engine analyses. The reduction names
# two walls: towards {k1 = 0} the orbit runs off to A -> Inf; {k2 = 0} is a
# boundary it meets at finite A and k1, and there the flat profile ends. Which
# one binds depends on the true parameters, not on the model.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryAB"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model and its symmetries
#
# The engine works over Q, so it gets the rational observable s*B; the log
# only enters the prediction chain below.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
reactions <- eqnlist() |>
  addReaction("A", "B", "k1*A", "A to B") |>
  addReaction("B", "A", "k2*B", "B to A")

observables <- eqnvec(y = "s*B")

res <- symmetryDetection(reactions, observables,
                         method = "observability", reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

# fixed = "s" is the gauge the fit fixes too. The zero limits further down
# need reportZeroCompatibility, which is off by default.
red <- symmetryReduction(res, fixed = "s", reportZeroCompatibility = TRUE)
summary(red)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# The ODE is integrated as it stands and only the parameters live on log10,
# which keeps every rate positive and puts {k1 = 0} at -Inf.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptionsODE  <- list(atol = 1e-10, rtol = 1e-10, maxsteps = 1e7)
myOptionsSens <- myOptionsODE

model <- odemodel(reactions, modelname = "AB_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "C1", optionsOde = myOptionsODE, optionsSens = myOptionsSens)

# attach.input carries the states along for the flow picture; the verdict and
# the objective read only the observable
g <- Y(eqnvec(y = "log2(s*B)"), f = x, attach.input = TRUE,
       modelname = "AB_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)

# The flow check needs every coordinate free, in the _l10 coordinates the flow
# is integrated in, so output and input meet without conversion.
p.flow <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x ~ exp10(x_l10)", x = .currentSymbols) |>
  P(modelname = "AB_pflow", compile = FALSE, outdir = .outdir)

# s gauged: the log2 offset is 0
p.lin <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("s~1") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "AB_plin", compile = FALSE, outdir = .outdir)

# red$trafo is stated in linear coordinates, so it goes in before exp10. That
# needs a positive carrier (red$blocks[[1]]$carrierDomain); a "[real-valued]"
# one would have to stay linear.
p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  define("s~1") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "AB_pred", compile = FALSE, outdir = .outdir)

# One flow per direction, dz/deps = sgn*xi(z): the sign runs the same compiled
# model both ways. log10Transform() renames only the moved coordinates, the
# rest enter linearly and keep their nominal value.
flow <- lapply(seq_along(res$symmetries), function(i) {
  fl <- log10Transform(res$symmetries[[i]]$completeGenerator)
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
              deriv = FALSE, modelname = paste0("AB_flow", i), compile = FALSE,
              outdir = .outdir),
     condition = "flow",
     optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
})

# One call compiles ODE, observation, the three trafos and the flows together.
do.call(compile, c(list(x, g, p.flow, p.lin, p.red), flow,
                   list(output = .modelname, cores = .cores)))

prd     <- g*x*p.flow   # every coordinate free
prd.lin <- g*x*p.lin    # s gauged
prd.red <- g*x*p.red    # s gauged, curved direction removed


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Zero limits
#
# Which coordinate the symmetry can switch off without another one running to
# infinity. The reduction states the condition, the truth decides it: with
# k1*A > k2*B, k2 goes and the back reaction with it. Swap the pair round
# (A = 0.2, B = 1, k1 = 0.05, k2 = 0.3) and A = 0 is the reachable zero
# instead. k1 never is, since k1*(A + B) is invariant and strictly positive.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth  <- c(A = 1, B = 0.2, k1 = 0.3, k2 = 0.05, s = 1)
truthL <- setNames(log10(truth), paste0(names(truth), "_l10"))

red$zeroCompatibility[, c("coordinates", "verdict", "limit", "condition", "at")]
cond <- setNames(red$zeroCompatibility$condition, red$zeroCompatibility$coordinates)
sapply(cond[c("A", "k2")], function(cc) eval(parse(text = cc), as.list(truth)))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# One direction, flown
#
# The flow runs along completeGenerator, so it goes on for every
# eps both ways: towards {k1 = 0} and towards the boundary {k2 = 0}, which it
# approaches without reaching.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
times <- seq(0, 10, len = 300)
pred0 <- as.data.frame(prd(times, truthL, deriv = FALSE))   # unmoved reference

direction <- 2   # click through 1 .. length(flow)
eps <- seq(0, 10, len = 50)

mv  <- names(res$symmetries[[direction]]$generator)
mvL <- paste0(mv, "_l10")
z0  <- c(truthL[mvL], truth[setdiff(coords, mv)])

zf <- flow[[direction]](eps, c(z0, sgn =  1), deriv = FALSE)$flow
zb <- flow[[direction]](eps, c(z0, sgn = -1), deriv = FALSE)$flow

# one orbit through the truth, indexed by signed eps taken from the returned
# time column, so a run the solver cut short stays labelled correctly
z <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])

ggplot(wide2long(z), aes(time, value)) +
  geom_line() + geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ name, scales = "free_y") +
  labs(x = expression(epsilon), y = NULL) + theme_dMod(base_size = 9)

theta <- lapply(seq_len(nrow(z)), function(j) replace(truthL, mvL, z[j, mvL]))
pred  <- lapply(theta, function(th) as.data.frame(prd(times, th, deriv = FALSE)))

# verdict: max deviation of the observable relative to the RMS of the unmoved
# curve
obs <- pred0$name %in% names(observables)
rms <- ave(pred0$value[obs], pred0$name[obs], FUN = function(v) sqrt(mean(v^2)))

sapply(pred, function(d) max(abs(d$value[obs] - pred0$value[obs]) / rms))

# the states move, the observable does not
long <- do.call(rbind, Map(function(d, e) transform(d, eps = e), pred, z[, "time"]))
lev  <- levels(long$name)
levels(long$name) <- ifelse(lev %in% names(observables), paste0(lev, " [observed]"), lev)
long$name <- factor(long$name, levels(long$name)[order(lev %in% names(observables))])

ggplot(long, aes(time, value, group = eps, colour = eps)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = paste0("A <-> B, X", direction), x = "time", y = NULL) +
  theme_dMod(base_size = 9)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Every direction, both ways
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
for (i in seq_along(flow)) {
  mv  <- names(res$symmetries[[i]]$generator)
  mvL <- paste0(mv, "_l10")
  z0  <- c(truthL[mvL], truth[setdiff(coords, mv)])
  for (sgn in c(1, -1)) {
    zi <- try(flow[[i]](eps, c(z0, sgn = sgn), deriv = FALSE)$flow, silent = TRUE)
    if (inherits(zi, "try-error")) {
      cat(sprintf("X%d  %-8s  eps %+.2f  escapes\n", i,
                  res$symmetries[[i]]$type, sgn * max(eps)))
      next
    }
    d <- sapply(seq_len(nrow(zi)), function(j) {
      q <- as.data.frame(prd(times, replace(truthL, mvL, zi[j, mvL]), deriv = FALSE))
      max(abs(q$value[obs] - pred0$value[obs]) / rms)
    })
    cat(sprintf("X%d  %-8s on %-12s  eps %+.2f  moved %.2f  deviates %.2e\n", i,
                res$symmetries[[i]]$type,
                paste(res$symmetries[[i]]$support, collapse = ","),
                sgn * max(zi[, "time"]),
                max(abs(10^(zi[nrow(zi), mvL] - truthL[mvL]) - 1)), max(d)))
  }
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Simulate data
#
# 12 replicates per time point with noise sd 0.3 on the log2 observable,
# pooled by reduceReplicates().
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth.fit <- log10(truth[getParameters(prd.lin)])
.sigma    <- 0.3
.nrep     <- 12

set.seed(1111)
sim <- as.data.frame(prd.lin(seq(0, 10, by = 0.8), truth.fit, deriv = FALSE))
sim <- sim[sim$name == "y", ]

mydataL <- data.frame(
  name = "y", time = rep(sim$time, each = .nrep),
  value = rep(sim$value, each = .nrep) + rnorm(.nrep * nrow(sim), 0, .sigma),
  sigma = NA, condition = "C1") |>
  reduceReplicates() |>
  as.datalist(split.by = "condition")

plot(prd.lin(times, truth.fit, deriv = FALSE), mydataL)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and profile, full model
#
# With s gauged one direction is left, and its profile is flat until the
# orbit reaches k2 = 0.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
stepControl <- list(stepsize = 1e-8, min = 1e-8, max = Inf,
                    atol = 1e-3, rtol = 1e-3, limit = 1e3)
optControl  <- list(rinit = 0.1, rmax = 10, iterlim = 2e3)
algoControl <- list(reoptimize = TRUE)

pouter <- structure(rep(-1, length(getParameters(prd.lin))),
                    names = getParameters(prd.lin))
obj    <- normL2(mydataL, prd.lin)
fit    <- mstrust(obj, pouter, rinit = 0.1, rmax = 10, sd = 4, fits = 100,
                  iterlim = 1e3, cores = .cores, studyname = "full")

plotValues(as.parframe(fit), tol = 0.1, value < 1e4)
bestfit <- as.parvec(as.parframe(fit))
plot(prd.lin(times, bestfit, deriv = FALSE), mydataL)

prof <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate",
                algoControl = algoControl,
                stepControl = stepControl, optControl = optControl,
                limits = c(lower = -5, upper = 5), cores = .cores)
plotProfilesAndPaths(prof, names(bestfit))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and profile, reduced model
#
# The reduction removes a direction, not a degree of freedom: the optima
# agree, only the reduced profiles are finite.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
pouter.red <- structure(rep(-1, length(getParameters(prd.red))),
                        names = getParameters(prd.red))
obj.red    <- normL2(mydataL, prd.red)
fit.red    <- mstrust(obj.red, pouter.red, rinit = 0.1, rmax = 10, sd = 4, fits = 100,
                      iterlim = 1e3, cores = .cores, studyname = "reduced")

plotValues(as.parframe(fit.red), tol = 0.1, value < 1e4)
bestfit.red <- as.parvec(as.parframe(fit.red))
plot(prd.red(times, bestfit.red, deriv = FALSE), mydataL)

prof.red <- profile(obj.red, bestfit.red, whichPar = names(bestfit.red),
                    method = "integrate", algoControl = algoControl,
                    stepControl = stepControl, optControl = optControl,
                    limits = c(lower = -10, upper = 10), cores = .cores)
plotProfilesAndPaths(prof.red, names(bestfit.red))

cat(sprintf("optimum full %.6f, reduced %.6f, difference %.2e\n",
            obj(bestfit)$value, obj.red(bestfit.red)$value,
            abs(obj.red(bestfit.red)$value - obj(bestfit)$value)))
