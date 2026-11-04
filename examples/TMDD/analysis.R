#
# The full TMDD model of an antibody binding a membrane-bound target, with a peripheral compartment:
# free drug binds the target, the complex dissociates or is internalised, the target is
# synthesised and degraded. The true parameters are the nominal ones of a fully human
# antibody and its membrane-bound target in Gibiansky and Gibiansky (2010), Table 1. No free TMDD
# data set exists, so the analysis simulates one: five intravenous doses, 45 subjects, eight
# samples each, free drug above a limit of quantification of 0.1 nmol/L. FOCEI, quadrature and,
# where nlmixr2 is installed, nlmixr2 FOCEI and SAEM fit it as multistarts of 100 fits.

library(dMod2)
dir.create("build", showWarnings = FALSE)
# generated sources and shared objects in build/, reused when nothing changed
options(dMod2.outdir = "build")
cores <- max(1L, min(20L, parallel::detectCores() - 2L))
nstart <- 100L
withNlmixr2 <- requireNamespace("nlmixr2", quietly = TRUE)

# fits are saved in results/cache.rds under a key and read from there when present; delete the
# file to fit everything again
cached <- function(key, expr) {
  file <- "results/cache.rds"
  dir.create("results", showWarnings = FALSE)
  store <- if (file.exists(file)) readRDS(file) else list()
  if (is.null(store[[key]])) {
    store[[key]] <- expr
    saveRDS(store, file)
  }
  store[[key]]
}
# a fit of em() cached without its model, which an evaluation at its estimates restores
cachedEm <- function(key, obj, random, start, control = list(), ...) {
  model <- c("design", "prdfn", "errfn", "data")
  f <- cached(key, { f <- em(obj, random, start, control = control, ...); f[model] <- NULL; f })
  f[model] <- em(obj, random, f$argument, control = modifyList(control, list(trustControl = list(iterlim = 0))),
                 verbose = FALSE)[model]
  f
}

dm <- c("forward", "reverse", "forward-reverse")

# Model ----
#
# Concentrations in nmol/L, time in days; the rates of distribution and linear elimination are
# fluxes in nmol per day in a frame of volume 1, binding and turnover concentration rates in the
# central compartment.
net <- eqnlist() |>
  assignCompartment(C = "central", volume = "V") |>
  assignCompartment(Cp = "peripheral", volume = "Vp") |>
  assignCompartment(R = "central") |>
  assignCompartment(RC = "central") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("C", "", "CL * C", "elimination", rateCompartment = "unit") |>
  addReaction("C", "Cp", "Q * C", "distribution", rateCompartment = "unit") |>
  addReaction("Cp", "C", "Q * Cp", "redistribution", rateCompartment = "unit") |>
  addReaction("C + R", "RC", "kon * C * R", "binding") |>
  addReaction("RC", "C + R", "koff * RC", "dissociation") |>
  addReaction("RC", "", "kint * RC", "internalisation") |>
  addReaction("", "R", "ksyn", "synthesis", compartment = "central") |>
  addReaction("R", "", "kdeg * R", "degradation")
net

# The target before the dose is at the steady state of the drug-free system, R = ksyn / kdeg;
# with R0 as parameter the synthesis is kdeg R0.
steady <- steadyStates(eqnlist() |> addReaction("", "R", "ksyn", "synthesis") |>
                         addReaction("R", "", "kdeg * R", "degradation"))
steady

# Identifiability ----
#
# With the dose known and the target at its steady state, R(0) = R0 and ksyn = kdeg R0, the free
# drug alone determines all nine parameters; total target is not measured.
start <- eqnvec(C = "dose / V", Cp = "0", R = "R0", RC = "0", ksyn = "kdeg * R0")
parameters <- c("CL", "V", "Q", "Vp", "kon", "koff", "kint", "kdeg", "R0")
identifiability <- symmetryDetection(as.eqnvec(net), eqnvec(lcp = "log(C)"), fixed = "dose",
                                     trafo = c(start, eqnvec(setNames(parameters, parameters)), eqnvec(dose = "dose")))
identifiability

# Design and simulation ----
#
# Five dose groups of 9 subjects, samples from one hour to eight weeks; clearance, central volume
# and target level vary between subjects with 20% CV, the readings with 15% on the log scale. The
# seed gives a typical draw: its share of censored readings is the median of seeds 1 to 50.
truth <- log(c(CL = 0.3, V = 3, Q = 0.2, Vp = 3, kon = 5, koff = 0.25, kint = 15, kdeg = 15, R0 = 0.1))
omega <- sqrt(log(1 + 0.2^2))
sigma <- sqrt(log(1 + 0.15^2))
doses <- c(10, 30, 100, 300, 1000)
times <- c(1 / 24, 1, 3, 7, 14, 21, 35, 56)
lloq <- 0.1
subjects <- sprintf("S%02d", 1:45)
table <- data.frame(ID = subjects, dose = rep(doses, each = 9), row.names = subjects)

# The free drug in log10 coordinates, C = 10^C_l10: it stays positive where binding takes it far
# below the absolute tolerance of the solver, and log(C) is linear in the state.
ode <- log10Transform(as.eqnvec(net), states = "C")
ode
# solver tolerances 1e-10: at the defaults the solver error moves -2 log L by about 1e-4, more than
# FOCEI resolves, and the fit stalls short of its optimum
x <- Xs(odemodel(ode, modelname = "tmdd", compile = FALSE, deriv2 = TRUE, derivMode = dm),
        options = list(atol = 1e-10, rtol = 1e-10))
g <- Y(c(lcp = "log(10) * C_l10"), x, modelname = "tmdd_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(lcp = "sigma"), g, attachInput = FALSE, modelname = "tmdd_e", compile = FALSE, deriv2 = TRUE, derivMode = dm)
trafo <- branch(define(NULL, "x <- x", x = c("CL", "V", "Q", "Vp", "kon", "koff", "kint", "ksyn", "kdeg",
                                             "C_l10", "Cp", "R", "RC", "sigma")), table = table)
trafo <- insert(trafo, "x <- y", x = c("C_l10", "Cp", "RC"), y = c("log10(dose / V)", "0", "0"))
trafo <- insert(trafo, paste("R <-", steady[["R"]]))
trafo <- insert(trafo, "ksyn <- kdeg * R0")
trafo <- insert(trafo, "dose <- d", d = dose)
trafo <- insert(trafo, "x <- exp(x)", x = c(parameters, "sigma"))
trafo <- insert(trafo, CL + V + R0 ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "tmdd_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(x, g, e, p, output = "tmdd_all")
random <- randomEffects(p)

simulateData <- function(seed) {
  set.seed(seed)
  eta <- matrix(rnorm(45 * 3, 0, omega), 45, dimnames = list(subjects, c("CL", "V", "R0")))
  pars <- c(truth, sigma = log(sigma),
            setNames(as.vector(eta), paste0("eta_", rep(colnames(eta), each = 45), "_", subjects)))
  pred <- (g * x * p)(times = times, pars = pars, deriv = FALSE)
  do.call(rbind, lapply(subjects, function(s) {
    v <- pred[[s]][match(times, pred[[s]][, "time"]), "lcp"] + rnorm(length(times), 0, sigma)
    data.frame(name = "lcp", time = times, value = pmax(v, log(lloq)), sigma = NA_real_, lloq = log(lloq),
               condition = s, ID = s, dose = table[s, "dose"])
  }))
}
censored <- vapply(1:50, function(s) mean(simulateData(s)$value <= log(lloq)), 0)
seed <- which.min(abs(censored - stats::median(censored)))
c(seed = seed, censored = censored[seed], median = stats::median(censored))
sim <- simulateData(seed)
data <- as.datalist(sim, splitBy = "condition", keepCovariates = c("ID", "dose"))
plotData(data) + ggplot2::aes(colour = factor(dose)) +
  ggplot2::labs(x = "time [days]", y = "ln free drug [nmol/L]", colour = "dose [nmol]")

# Fits ----
#
# The starts: the truth on the log scale and 99 draws around it, N(0, 0.5) per parameter. Every
# fit runs on one core. Up to level 5 the nodes resolve -2 log L only to the level error the
# summary of quadrature prints, and a rebuild of the nodes moves it by about 1e-3. A quadrature
# fit converges where a step no longer lowers -2 log L on the rebuilt nodes; one still moving
# after ten rebuilds ends unconverged.
obj <- normL2(data, g * x * p, errmodel = e)
# second derivatives forward over forward; quadrature up to level 5, which the certificate passes
ctrl <- list(derivMode = "forward-forward", quadrature = list(maxLevel = 5))
# -2 log L of FOCEI at given estimates, the definition under which all fits are compared
focei <- function(obj, random, pars)
  em(obj, random, pars, control = c(ctrl, list(trustControl = list(iterlim = 0))), verbose = FALSE)$value
init <- c(truth, sigma = log(sigma), omega_CL_CL = log(omega), omega_V_V = log(omega), omega_R0_R0 = log(omega))
set.seed(1)
starts <- rbind(init, t(vapply(seq_len(nstart - 1L), function(i) init + rnorm(length(init), 0, 0.5), init)))
runs <- cached("dMod2", do.call(rbind, unlist(lapply(c(focei = "focei", quadrature = "quadrature"), function(method)
  parallel::mclapply(seq_len(nrow(starts)), function(i) {
    options(dMod2.cores = 1L)
    t <- system.time(f <- tryCatch(em(obj, random, starts[i, ], method = method, control = ctrl, verbose = FALSE),
                                   error = function(e) NULL))[["elapsed"]]
    data.frame(method = method, start = i, value = if (is.null(f)) NA_real_ else f$value,
               focei = if (is.null(f)) NA_real_ else if (method == "focei") f$value else focei(obj, random, f$argument),
               converged = !is.null(f) && isTRUE(f$converged), seconds = t)
  }, mc.cores = cores, mc.preschedule = FALSE)), recursive = FALSE)))
best <- tapply(runs$value, runs$method, min, na.rm = TRUE)
data.frame(best = best, hits = tapply(runs$value - best[runs$method] < 0.1, runs$method, sum, na.rm = TRUE),
           converged = tapply(runs$converged, runs$method, sum), failed = tapply(is.na(runs$value), runs$method, sum),
           medianSeconds = tapply(runs$seconds, runs$method, stats::median))
fit <- cachedEm("dMod2 focei best", obj, random, starts[which.min(runs$value[runs$method == "focei"]), ],
                control = ctrl, verbose = FALSE)
quad <- cachedEm("dMod2 quadrature best", obj, random, starts[which.min(runs$value[runs$method == "quadrature"]), ],
                 method = "quadrature", verbose = FALSE,
                 control = modifyList(ctrl, list(quadrature = list(certificate = list(points = 256)))))
summary(fit)
quad
round(cbind(truth = init, focei = fit$argument[names(init)], quadrature = quad$argument[names(init)],
            se = sqrt(pmax(diag(vcov(fit)), 0))[names(init)]), 3)

# Diagnostics ----
#
# CWRES linearise around the modes, NPDE come from simulations of the fitted model and take the
# censored readings by their simulated fraction.
plotCombined(fit, conditionsPerPage = 15)[[1]]
plotResiduals(fit, type = c("CWRES", "NPDE"), nsim = 500, seed = 1)
attr(residuals(fit, type = "NPDE", nsim = 500, seed = 1), "tests")
plotRandomEffects(fit)
plotRandomEffects(fit, covariates = "dose")
plotVPC(fit, nsim = 500, bins = c(0, (head(times, -1) + tail(times, -1)) / 2, 60), pc = TRUE, logScale = TRUE,
        strata = "dose", seed = 1)

# nlmixr2 FOCEI and SAEM ----
#
# The same model with the censored readings by M3 (CENS = 1, DV the limit) and the same starts.
# Every fit, quadrature included, counts by the -2 log L of FOCEI at its estimates.
if (withNlmixr2) {
  ui <- nlmixr2est::nlmixr2(function() {
    ini({
      lcl <- -1.2; lv <- 1.1; lq <- -1.6; lvp <- 1.1; lkon <- 1.6; lkoff <- -1.4; lkint <- 2.7; lkdeg <- 2.7
      lr0 <- -2.3
      eta.cl ~ 0.04; eta.v ~ 0.04; eta.r0 ~ 0.04
      add.sd <- 0.15
    })
    model({
      cl <- exp(lcl + eta.cl); v <- exp(lv + eta.v); q <- exp(lq); vp <- exp(lvp)
      kon <- exp(lkon); koff <- exp(lkoff); kint <- exp(lkint); kdeg <- exp(lkdeg); r0 <- exp(lr0 + eta.r0)
      cf <- central / v
      d/dt(central) <- -cl * cf - q * cf + q * peripheral / vp - (kon * cf * r - koff * rc) * v
      d/dt(peripheral) <- q * cf - q * peripheral / vp
      r(0) <- r0
      d/dt(r) <- kdeg * r0 - kdeg * r - kon * cf * r + koff * rc
      d/dt(rc) <- kon * cf * r - koff * rc - kint * rc
      lcp <- log(cf)
      lcp ~ add(add.sd)
    })
  })
  ndata <- rbind(data.frame(ID = match(table$ID, subjects), TIME = 0, AMT = table$dose, EVID = 1L, DV = NA_real_,
                            CENS = 0L, CMT = "central"),
                 data.frame(ID = match(sim$ID, subjects), TIME = sim$time, AMT = 0, EVID = 0L, DV = sim$value,
                            CENS = as.integer(sim$value <= log(lloq)), CMT = "central"))
  ndata <- ndata[order(ndata$ID, ndata$TIME, -ndata$EVID), ]
  toIni <- function(s)
    c(list(lcl = s[["CL"]], lv = s[["V"]], lq = s[["Q"]], lvp = s[["Vp"]], lkon = s[["kon"]], lkoff = s[["koff"]],
           lkint = s[["kint"]], lkdeg = s[["kdeg"]], lr0 = s[["R0"]], add.sd = exp(s[["sigma"]])),
      unname(Map(function(eta, o) stats::as.formula(sprintf("%s ~ %.10g", eta, exp(2 * s[[o]]))),
                 c("eta.cl", "eta.v", "eta.r0"), c("omega_CL_CL", "omega_V_V", "omega_R0_R0"))))
  fromNlmixr2 <- function(f) {
    th <- f$theta; om <- diag(f$omega)
    c(CL = th[["lcl"]], V = th[["lv"]], Q = th[["lq"]], Vp = th[["lvp"]], kon = th[["lkon"]], koff = th[["lkoff"]],
      kint = th[["lkint"]], kdeg = th[["lkdeg"]], R0 = th[["lr0"]], sigma = log(th[["add.sd"]]),
      omega_CL_CL = log(sqrt(om[["eta.cl"]])), omega_V_V = log(sqrt(om[["eta.v"]])),
      omega_R0_R0 = log(sqrt(om[["eta.r0"]])))
  }
  # rxode2::ini() evaluates the random effects in `envir`, where it needs lotri()
  nfit <- function(s, est) {
    out <- NULL
    control <- if (est == "focei") nlmixr2est::foceiControl(print = 0) else nlmixr2est::saemControl(print = 0)
    utils::capture.output(out <- tryCatch(suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr2(do.call(rxode2::ini, c(list(ui), toIni(s), envir = asNamespace("rxode2"))), ndata, est = est,
                          control = control))),
      error = function(e) NULL))
    out
  }
  rxode2::setRxThreads(1L)
  nruns <- cached("nlmixr2", {
    invisible(nfit(starts[1, ], "focei"))
    do.call(rbind, unlist(lapply(c("focei", "saem"), function(est)
      parallel::mclapply(seq_len(nrow(starts)), function(i) {
        rxode2::setRxThreads(1L)
        t <- system.time(f <- nfit(starts[i, ], est))[["elapsed"]]
        data.frame(method = paste("nlmixr2", est), start = i, seconds = t,
                   value = if (is.null(f)) NA_real_ else focei(obj, random, fromNlmixr2(f)))
      }, mc.cores = cores, mc.preschedule = FALSE)), recursive = FALSE))
  })
  all <- rbind(data.frame(method = runs$method, start = runs$start, value = runs$focei, seconds = runs$seconds), nruns)
  all$method <- factor(all$method, c("focei", "quadrature", "nlmixr2 focei", "nlmixr2 saem"),
                       c("dMod2 FOCEI", "dMod2 quadrature", "nlmixr2 FOCEI", "nlmixr2 SAEM"))
  all$delta <- all$value - min(all$value, na.rm = TRUE)
  all$rank <- stats::ave(rank(all$value, ties.method = "first", na.last = TRUE), all$method,
                         FUN = function(r) rank(r, na.last = TRUE))
  print(ggplot2::ggplot(all[is.finite(all$value), ], ggplot2::aes(rank, delta + 1e-3, color = method,
                                                                  shape = method)) +
          ggplot2::geom_point() + ggplot2::scale_y_log10(labels = scales::label_log()) +
          ggplot2::scale_shape_manual(values = c(19, 17, 1, 2)) + scale_color_dMod() + theme_dMod() +
          ggplot2::labs(x = "fit, ordered", y = "-2 log L (FOCEI) above the best + 0.001", color = NULL,
                        shape = NULL))
  print(do.call(rbind, lapply(split(all, all$method), function(r)
    data.frame(method = r$method[1], best = min(r$value, na.rm = TRUE),
               hits = sum(r$value - min(all$value, na.rm = TRUE) < 0.1, na.rm = TRUE),
               failed = sum(is.na(r$value)), medianSeconds = stats::median(r$seconds), totalSeconds = sum(r$seconds)))))
  nbest <- cached("nlmixr2 best", lapply(c(focei = "nlmixr2 focei", saem = "nlmixr2 saem"), function(m)
    fromNlmixr2(nfit(starts[nruns$start[nruns$method == m][which.min(nruns$value[nruns$method == m])], ],
                     sub("nlmixr2 ", "", m)))))
  print(round(cbind(truth = init, dMod2 = fit$argument[names(init)], nlmixr2 = nbest$focei[names(init)],
                    nlmixr2Saem = nbest$saem[names(init)]), 3))
}
