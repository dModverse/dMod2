# -------------------------------------------------------------------------#
# A <-> B: what a symmetry does to a profile likelihood
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The smallest model with a non-scaling symmetry: A converts to B and back,
# only s*B is observed. symmetryDetection() finds its directions, each is
# flown both ways to show that the states move and the observable does not.
# The full and the reduced model from symmetryReduction() are then fitted and
# profiled side by side.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 23 Sep 2026
#
# [Info]
# Step through top to bottom. Compiled code goes to tempdir(), nothing else
# is written to disk.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "symmetryAB"
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model and its symmetries
#
# The engine works over Q, so it gets the rational observable; the log
# enters only in the prediction chain.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
reactions <- eqnlist() |>
  addReaction("A", "B", "k1*A", "A to B") |>
  addReaction("B", "A", "k2*B", "B to A")

observables <- eqnvec(y = "s*B")

res <- symmetryDetection(reactions, observables, reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

# s is the gauge the fit fixes as well; the zero limits below need
# reportZeroCompatibility
red <- symmetryReduction(res, fixed = "s", reportZeroCompatibility = TRUE)
summary(red)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Prediction chain
#
# The ODE is integrated in linear states, the parameters live on log10.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptions <- list(atol = 1e-10, rtol = 1e-10, maxsteps = 1e7)

model <- odemodel(reactions, modelname = "AB_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, condition = "C1", options = myOptions)

# attachInput adds the states to the output for the flow plots; the
# objective reads only the observable
g <- Y(eqnvec(y = "log2(s*B)"), f = x, attachInput = TRUE,
       modelname = "AB_obs", compile = FALSE, outdir = .outdir)

innerpars <- getParameters(g, x)

# every coordinate free, named _l10 like the coordinates of the flow
p.flow <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x ~ exp10(x_l10)", x = .currentSymbols) |>
  P(modelname = "AB_pflow", compile = FALSE, outdir = .outdir)

# s gauged
p.lin <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("s~1") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "AB_plin", compile = FALSE, outdir = .outdir)

# red$trafo is stated in linear coordinates and goes in before exp10, which
# requires a positive carrier domain (red$blocks[[1]]$carrierDomain)
red$blocks[[1]]$carrierDomain
p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  define("s~1") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "AB_pred", compile = FALSE, outdir = .outdir)

# One flow dz/deps = sgn*xi(z) per direction, so one compiled model runs both
# ways. log10Transform() renames the moved coordinates; the others stay linear.
flow <- lapply(seq_along(res$symmetries), function(i) {
  fl <- log10Transform(res$symmetries[[i]]$completeGenerator)
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
              deriv = FALSE, modelname = paste0("AB_flow", i), compile = FALSE,
              outdir = .outdir),
     condition = "flow",
     options = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
})

do.call(compile, c(list(x, g, p.flow, p.lin, p.red), flow,
                   list(output = .modelname, cores = .cores)))

prd     <- g*x*p.flow   # every coordinate free
prd.lin <- g*x*p.lin    # s gauged
prd.red <- g*x*p.red    # s gauged, symmetry directions removed


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Zero limits
#
# Which coordinates the symmetry can drive to zero, and under which condition
# on the parameters. The condition is evaluated at the true parameters.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth  <- c(A = 1, B = 0.2, k1 = 0.3, k2 = 0.05, s = 1)
# truth <- c(A = 0.2, B = 1, k1 = 0.05, k2 = 0.3, s = 1)   # the other regime
truthL <- setNames(log10(truth), paste0(names(truth), "_l10"))

red$zeroCompatibility[, c("coordinates", "verdict", "limit", "condition", "at")]
cond <- setNames(red$zeroCompatibility$condition, red$zeroCompatibility$coordinates)
sapply(cond[c("A", "k2")], function(cc) eval(parse(text = cc), as.list(truth)))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# One direction, flown
#
# The flow follows completeGenerator from the truth, both ways in eps.
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

# one orbit through the truth, signed eps taken from the returned time column
# so that a run the solver stopped early is labelled correctly
z <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])

ggplot(wide2long(z), aes(time, value)) +
  geom_line() + geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ name, scales = "free_y") +
  labs(x = expression(epsilon), y = NULL) + theme_dMod(base_size = 9)

theta <- lapply(seq_len(nrow(z)), function(j) replace(truthL, mvL, z[j, mvL]))
pred  <- lapply(theta, function(th) as.data.frame(prd(times, th, deriv = FALSE)))

# deviation of the observable along the orbit, relative to the RMS of the
# unmoved curve
obs <- pred0$name %in% names(observables)
rms <- ave(pred0$value[obs], pred0$name[obs], FUN = function(v) sqrt(mean(v^2)))

deviation <- sapply(pred, function(d) max(abs(d$value[obs] - pred0$value[obs]) / rms))
deviation

# states and observable along the orbit, observed panels last
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
#
# One row per direction and sign: how far the moved coordinates got and the
# largest relative deviation of the observable. NA marks a flow that failed.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
runs <- expand.grid(direction = seq_along(flow), sgn = c(1, -1))

sweep <- do.call(rbind, lapply(seq_len(nrow(runs)), function(r) {
  i   <- runs$direction[r]
  sgn <- runs$sgn[r]
  mv  <- names(res$symmetries[[i]]$generator)
  mvL <- paste0(mv, "_l10")
  z0  <- c(truthL[mvL], truth[setdiff(coords, mv)])
  out <- data.frame(direction = i, type = res$symmetries[[i]]$type,
                    support = paste(res$symmetries[[i]]$support, collapse = ","),
                    eps = sgn * max(eps), moved = NA, deviates = NA)
  zi <- try(flow[[i]](eps, c(z0, sgn = sgn), deriv = FALSE)$flow, silent = TRUE)
  if (inherits(zi, "try-error")) return(out)
  d <- sapply(seq_len(nrow(zi)), function(j) {
    q <- as.data.frame(prd(times, replace(truthL, mvL, zi[j, mvL]), deriv = FALSE))
    max(abs(q$value[obs] - pred0$value[obs]) / rms)
  })
  out$eps      <- sgn * max(zi[, "time"])
  out$moved    <- max(abs(10^(zi[nrow(zi), mvL] - truthL[mvL]) - 1))
  out$deviates <- max(d)
  out
}))
sweep


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Simulate data
#
# Replicates per time point with Gaussian noise on the log2 observable,
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
  as.datalist(splitBy = "condition")

plot(prd.lin(times, truth.fit, deriv = FALSE), mydataL)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Fit and profile, full model
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
# Same data and settings on the reduced parameters; the optima should agree.
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

value.full <- obj(bestfit)$value
value.red  <- obj.red(bestfit.red)$value
c(full = value.full, reduced = value.red, difference = abs(value.red - value.full))
