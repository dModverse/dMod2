## Full TMDD with a peripheral compartment, the design and truth of demo("TMDD") (the membrane-bound
## target of Gibiansky and Gibiansky 2010, Table 1), over `nsets` simulated data sets: per data set
## FOCEI and quadrature of em() and nlmixr2 FOCEI and SAEM, each the best of a multistart of
## `nstart` fits from the same starts, one core per fit. Per method and parameter the bias, the
## RMSE and the coverage of the 95% Wald intervals (nlmixr2 reports no standard errors of Omega).
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "tmddBench"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)
nsets <- 20
nstart <- 100
cores <- max(1L, min(20L, parallel::detectCores() - 2L))
dm <- c("forward", "reverse", "forward-reverse")

## Model and design ----------------------------------------------------------------------
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
  addReaction("RC", "", "kint * RC", "internalization") |>
  addReaction("", "R", "ksyn", "synthesis", compartment = "central") |>
  addReaction("R", "", "kdeg * R", "degradation")
truth <- log(c(CL = 0.3, V = 3, Q = 0.2, Vp = 3, kon = 5, koff = 0.25, kint = 15, kdeg = 15, R0 = 0.1))
omega <- sqrt(log(1 + 0.2^2))
sigma <- sqrt(log(1 + 0.15^2))
init <- c(truth, sigma = log(sigma), omega_CL_CL = log(omega), omega_V_V = log(omega), omega_R0_R0 = log(omega))
doses <- c(10, 30, 100, 300, 1000)
times <- c(1 / 24, 1, 3, 7, 14, 21, 35, 56)
lloq <- 0.1
subjects <- sprintf("S%02d", 1:45)
table <- data.frame(ID = subjects, dose = rep(doses, each = 9), row.names = subjects)
# the free drug in log10 coordinates and the solver tolerances 1e-10, as in the demo
x <- Xs(odemodel(log10Transform(as.eqnvec(net), states = "C"), modelname = "tmddBench", compile = FALSE,
                 deriv2 = TRUE, derivMode = dm), options = list(atol = 1e-10, rtol = 1e-10))
g <- Y(c(lcp = "log(10) * C_l10"), x, modelname = "tmddBench_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(lcp = "sigma"), g, attachInput = FALSE, modelname = "tmddBench_e", compile = FALSE, deriv2 = TRUE, derivMode = dm)
trafo <- branch(define(NULL, "x <- x", x = c("CL", "V", "Q", "Vp", "kon", "koff", "kint", "ksyn", "kdeg",
                                             "C_l10", "Cp", "R", "RC", "sigma")), table = table)
trafo <- insert(trafo, "x <- y", x = c("C_l10", "Cp", "RC"), y = c("log10(dose / V)", "0", "0"))
trafo <- insert(trafo, "R <- ksyn / kdeg")
trafo <- insert(trafo, "ksyn <- kdeg * R0")
trafo <- insert(trafo, "dose <- d", d = dose)
trafo <- insert(trafo, "x <- exp(x)", x = c(names(truth), "sigma"))
trafo <- insert(trafo, CL + V + R0 ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "tmddBench_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(x, g, e, p, output = "tmddBench")
random <- randomEffects(p)
# second derivatives forward over forward; quadrature up to level 5
ctrl <- list(derivMode = "forward-forward", quadrature = list(maxLevel = 5))

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
toIni <- function(s)
  c(list(lcl = s[["CL"]], lv = s[["V"]], lq = s[["Q"]], lvp = s[["Vp"]], lkon = s[["kon"]], lkoff = s[["koff"]],
         lkint = s[["kint"]], lkdeg = s[["kdeg"]], lr0 = s[["R0"]], add.sd = exp(s[["sigma"]])),
    unname(Map(function(eta, o) stats::as.formula(sprintf("%s ~ %.10g", eta, exp(2 * s[[o]]))),
               c("eta.cl", "eta.v", "eta.r0"), c("omega_CL_CL", "omega_V_V", "omega_R0_R0"))))
thetaNames <- c(CL = "lcl", V = "lv", Q = "lq", Vp = "lvp", kon = "lkon", koff = "lkoff", kint = "lkint",
                kdeg = "lkdeg", R0 = "lr0")

## Data sets -----------------------------------------------------------------------------
results <- list()
for (k in seq_len(nsets)) {
  set.seed(1000 + k)
  eta <- matrix(rnorm(45 * 3, 0, omega), 45, dimnames = list(subjects, c("CL", "V", "R0")))
  pars <- c(truth, sigma = log(sigma), setNames(as.vector(eta), paste0("eta_", rep(colnames(eta), each = 45), "_", subjects)))
  pred <- (g * x * p)(times = times, pars = pars, deriv = FALSE)
  sim <- do.call(rbind, lapply(subjects, function(s) {
    v <- pred[[s]][match(times, pred[[s]][, "time"]), "lcp"] + rnorm(length(times), 0, sigma)
    data.frame(name = "lcp", time = times, value = pmax(v, log(lloq)), sigma = NA_real_, lloq = log(lloq),
               condition = s, ID = s, dose = table[s, "dose"])
  }))
  obj <- normL2(as.datalist(sim, splitBy = "condition", keepCovariates = c("ID", "dose")), g * x * p, errmodel = e)
  starts <- rbind(init, t(vapply(seq_len(nstart - 1), function(i) init + rnorm(length(init), 0, 0.5), init)))
  for (method in c("focei", "quadrature")) {
    values <- unlist(parallel::mclapply(seq_len(nstart), function(i) {
      options(dMod2.cores = 1L)
      f <- tryCatch(em(obj, random, starts[i, ], method = method, control = ctrl, verbose = FALSE),
                    error = function(e) NULL)
      if (is.null(f)) NA_real_ else f$value
    }, mc.cores = cores))
    best <- em(obj, random, starts[which.min(values), ], method = method, control = ctrl, verbose = FALSE)
    se <- sqrt(pmax(diag(vcov(best)), 0))
    results[[length(results) + 1]] <- data.frame(set = k, method = paste("dMod2", method), parameter = names(init),
                                                 estimate = unname(best$argument[names(init)]), se = unname(se[names(init)]))
  }
  ndata <- rbind(data.frame(ID = match(table$ID, subjects), TIME = 0, AMT = table$dose, EVID = 1L, DV = NA_real_,
                            CENS = 0L, CMT = "central"),
                 data.frame(ID = match(sim$ID, subjects), TIME = sim$time, AMT = 0, EVID = 0L, DV = sim$value,
                            CENS = as.integer(sim$value <= log(lloq)), CMT = "central"))
  ndata <- ndata[order(ndata$ID, ndata$TIME, -ndata$EVID), ]
  for (est in c("focei", "saem")) {
    control <- if (est == "focei") nlmixr2est::foceiControl(print = 0) else nlmixr2est::saemControl(print = 0)
    # rxode2::ini() evaluates the random effects in `envir`, where it needs lotri()
    nfit <- function(s) {
      out <- NULL
      utils::capture.output(out <- tryCatch(suppressMessages(suppressWarnings(
        nlmixr2est::nlmixr2(do.call(rxode2::ini, c(list(ui), toIni(s), envir = asNamespace("rxode2"))), ndata, est = est,
                            control = control))),
        error = function(e) NULL))
      out
    }
    rxode2::setRxThreads(1L)
    invisible(nfit(starts[1, ]))
    values <- unlist(parallel::mclapply(seq_len(nstart), function(i) {
      rxode2::setRxThreads(1L)
      f <- nfit(starts[i, ])
      if (is.null(f)) NA_real_ else -2 * as.numeric(suppressMessages(stats::logLik(f)))
    }, mc.cores = cores))
    f <- nfit(starts[which.min(values), ])
    th <- f$theta; om <- diag(f$omega); se <- setNames(f$parFixedDf$SE, rownames(f$parFixedDf))
    estimate <- c(th[thetaNames], log(th[["add.sd"]]), log(sqrt(om[c("eta.cl", "eta.v", "eta.r0")])))
    results[[length(results) + 1]] <- data.frame(set = k, method = paste("nlmixr2", est), parameter = names(init),
                                                 estimate = unname(estimate),
                                                 se = c(unname(se[thetaNames]), se[["add.sd"]] / th[["add.sd"]], NA, NA, NA))
  }
}

## Bias, RMSE and coverage ---------------------------------------------------------------
res <- do.call(rbind, results)
res$error <- res$estimate - init[res$parameter]
res$covered <- abs(res$error) < 1.96 * res$se
summary <- do.call(rbind, lapply(split(res, list(res$method, res$parameter), drop = TRUE), function(r)
  data.frame(method = r$method[1], parameter = r$parameter[1], bias = mean(r$error), rmse = sqrt(mean(r$error^2)),
             coverage = mean(r$covered))))
summary[order(summary$parameter, summary$method), ]
ggplot2::ggplot(res, ggplot2::aes(method, error)) + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_boxplot() + ggplot2::facet_wrap(~ parameter, scales = "free_y") + theme_dMod() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

## Single fits on 10 cores ---------------------------------------------------------------
# On the data set of demo("TMDD"), seed 10: one FOCEI fit per package from five starts halfway out,
# the optimum moved by 0.3 on the log scale in every parameter with random signs, on `threads`
# cores: dMod2 over the conditions, nlmixr2 over the subjects. Per start the minimum seconds over
# `reps` runs, and the -2 log L of both by FOCEI of dMod2 above the optimum.
threads <- 10L
reps <- 3
set.seed(10)
eta <- matrix(rnorm(45 * 3, 0, omega), 45, dimnames = list(subjects, c("CL", "V", "R0")))
pars <- c(truth, sigma = log(sigma), setNames(as.vector(eta), paste0("eta_", rep(colnames(eta), each = 45), "_", subjects)))
pred <- (g * x * p)(times = times, pars = pars, deriv = FALSE)
sim <- do.call(rbind, lapply(subjects, function(s) {
  v <- pred[[s]][match(times, pred[[s]][, "time"]), "lcp"] + rnorm(length(times), 0, sigma)
  data.frame(name = "lcp", time = times, value = pmax(v, log(lloq)), sigma = NA_real_, lloq = log(lloq),
             condition = s, ID = s, dose = table[s, "dose"])
}))
obj <- normL2(as.datalist(sim, splitBy = "condition", keepCovariates = c("ID", "dose")), g * x * p, errmodel = e)
ndata <- rbind(data.frame(ID = match(table$ID, subjects), TIME = 0, AMT = table$dose, EVID = 1L, DV = NA_real_,
                          CENS = 0L, CMT = "central"),
               data.frame(ID = match(sim$ID, subjects), TIME = sim$time, AMT = 0, EVID = 0L, DV = sim$value,
                          CENS = as.integer(sim$value <= log(lloq)), CMT = "central"))
ndata <- ndata[order(ndata$ID, ndata$TIME, -ndata$EVID), ]
nfit <- function(s) {
  out <- NULL
  utils::capture.output(out <- tryCatch(suppressMessages(suppressWarnings(
    nlmixr2est::nlmixr2(do.call(rxode2::ini, c(list(ui), toIni(s), envir = asNamespace("rxode2"))), ndata, est = "focei",
                        control = nlmixr2est::foceiControl(print = 0)))),
    error = function(e) NULL))
  out
}
fromNlmixr2 <- function(f) {
  th <- f$theta; om <- diag(f$omega)
  setNames(c(th[thetaNames], log(th[["add.sd"]]), log(sqrt(om[c("eta.cl", "eta.v", "eta.r0")]))), names(init))
}
focei <- function(pars) em(obj, random, pars, control = c(ctrl, list(trustControl = list(iterlim = 0))), verbose = FALSE)$value
optimum <- em(obj, random, init, control = ctrl, verbose = FALSE)
set.seed(2)
halfway <- t(vapply(1:5, function(i) optimum$argument[names(init)] + 0.3 * sample(c(-1, 1), length(init), TRUE), init))
rxode2::setRxThreads(threads)
invisible(nfit(init))
multi <- do.call(rbind, lapply(seq_len(nrow(halfway)), function(i) {
  d <- replicate(reps, {
    t <- system.time(f <- em(obj, random, halfway[i, ], control = ctrl, cores = threads, verbose = FALSE))[["elapsed"]]
    c(t, f$value)
  })
  n <- replicate(reps, {
    t <- system.time(f <- nfit(halfway[i, ]))[["elapsed"]]
    c(t, if (is.null(f)) NA_real_ else focei(fromNlmixr2(f)))
  })
  data.frame(start = i, dMod2 = min(d[1, ]), nlmixr2 = min(n[1, ]), dMod2Above = d[2, 1] - optimum$value,
             nlmixr2Above = n[2, 1] - optimum$value)
}))
rxode2::setRxThreads(1L)
multi
