## Mavoglurant, two compartments, an infusion of 10 min as a time switch, weight allometric with
## exponents 0.75 and 1, random effects on CL, V1, Q and V2 per subject and on CL per occasion:
## FOCEI of em() against nlmixr2 FOCEI, both with the variability between occasions, and SAEM of
## both without it, each a multistart of `nstart` fits from the same starts, one core per fit. Per
## method the best -2 log L (FOCEI by dMod2 at the estimates of every fit), the starts within 0.1 of
## it, the failed fits, and the time of the fit from the best start, the minimum over `reps` runs.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "mavoglurantBench"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)
reps <- 3
nstart <- 100
cores <- max(1L, min(20L, parallel::detectCores() - 2L))
dm <- c("forward", "reverse", "forward-reverse")

## dMod2 ---------------------------------------------------------------------------------
data(mavoglurant)
data <- as.datalist(mavoglurant, splitBy = "condition",
                    keepCovariates = c("ID", "OCC", "dose", "duration", "AGE", "SEX", "WT", "HT"))
net <- eqnlist() |>
  assignCompartment(C1 = "central", volume = "V1") |>
  assignCompartment(C2 = "peripheral", volume = "V2") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("", "C1", "piecewise(dose / duration, time < duration, 0)", "infusion", rateCompartment = "unit") |>
  addReaction("C1", "", "CL * C1", "elimination", rateCompartment = "unit") |>
  addReaction("C1", "C2", "Q * C1", "distribution", rateCompartment = "unit") |>
  addReaction("C2", "C1", "Q * C2", "redistribution", rateCompartment = "unit")
# solver tolerances 1e-10, as in the demo
x <- Xs(odemodel(net, modelname = "mavoBench", compile = FALSE, deriv2 = TRUE, derivMode = dm),
        options = list(atol = 1e-10, rtol = 1e-10))
g <- Y(c(cp = "1000 * C1"), x, modelname = "mavoBench_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(cp = "sqrt(sigma_abs^2 + (sigma_rel * cp)^2)"), g, attachInput = FALSE, modelname = "mavoBench_e",
       compile = FALSE, deriv2 = TRUE, derivMode = dm)
table <- attr(data, "condition.grid")
table$ID <- as.character(table$ID)
trafo <- branch(define(NULL, "x <- x", x = c("CL", "V1", "Q", "V2", "C1", "C2", "dose", "duration", "sigma_abs",
                                             "sigma_rel")), table = table)
trafo <- insert(trafo, "dose <- d", d = dose)
trafo <- insert(trafo, "duration <- u", u = duration)
trafo <- insert(trafo, "x <- 0", x = c("C1", "C2"))
trafo <- insert(trafo, "x <- x * (w / 70)^0.75", x = c("CL", "Q"), w = WT)
trafo <- insert(trafo, "x <- x * w / 70", x = c("V1", "V2"), w = WT)
trafo <- insert(trafo, "x <- exp(x)", x = c("CL", "V1", "Q", "V2", "sigma_abs", "sigma_rel"))
trafo1 <- insert(trafo, CL ~ 1 + (1 | ID))
trafo1 <- insert(trafo1, V1 + Q + V2 ~ 1 + (1 | ID))
trafo <- insert(trafo, CL ~ 1 + (1 | ID) + (1 | ID:OCC))
trafo <- insert(trafo, V1 + Q + V2 ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "mavoBench_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
p1 <- P(trafo1, method = "explicit", modelname = "mavoBench_p1", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compileSeconds <- system.time(compile(x, g, e, p, p1, output = "mavoBench"))[["elapsed"]]
obj <- normL2(data, g * x * p, errmodel = e)
obj1 <- normL2(data, g * x * p1, errmodel = e)
random <- randomEffects(p)
random1 <- randomEffects(p1)
init <- c(CL = log(30), V1 = log(40), Q = log(30), V2 = log(150), sigma_abs = log(1), sigma_rel = log(0.2),
          omega_CL_CL = log(0.3),
          omega_V1_V1 = log(0.3), omega_Q_Q = log(0.3), omega_V2_V2 = log(0.3), omega_CL_ID_OCC_CL_ID_OCC = log(0.15))
init1 <- init[-length(init)]
focei <- function(pars) em(obj, random, pars, control = list(trustControl = list(iterlim = 0)), verbose = FALSE)$value
focei1 <- function(pars) em(obj1, random1, pars, control = list(trustControl = list(iterlim = 0)), verbose = FALSE)$value
set.seed(1)
starts <- rbind(init, t(vapply(seq_len(nstart - 1), function(i) init + rnorm(length(init), 0, 0.5), init)))
starts1 <- starts[, names(init1)]

dmod <- rbind(
  do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
    options(dMod2.cores = 1L)
    f <- tryCatch(em(obj, random, starts[i, ], verbose = FALSE), error = function(e) NULL)
    data.frame(package = "dMod2", method = "focei", start = i, value = if (is.null(f)) NA_real_ else f$value)
  }, mc.cores = cores)),
  do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
    options(dMod2.cores = 1L)
    set.seed(i)
    f <- tryCatch(em(obj1, random1, starts1[i, ], method = "saem", verbose = FALSE), error = function(e) NULL)
    data.frame(package = "dMod2", method = "saem", start = i, value = if (is.null(f)) NA_real_ else focei1(f$argument[names(init1)]))
  }, mc.cores = cores)))

## nlmixr2 -------------------------------------------------------------------------------
# The occasions of a subject on one time axis, the second 1000 h after the first, the dose of an occasion a reset
# of the compartments (EVID 4); the variability between occasions as `iov.cl ~ ... | occ`.
ui <- nlmixr2est::nlmixr2(function() {
  ini({
    tcl <- log(30); tv1 <- log(40); tq <- log(30); tv2 <- log(150)
    eta.cl ~ 0.09; eta.v1 ~ 0.09; eta.q ~ 0.09; eta.v2 ~ 0.09
    iov.cl ~ 0.0225 | occ
    add.sd <- 1; prop.sd <- 0.2
  })
  model({
    cl <- exp(tcl + eta.cl + iov.cl) * (WT / 70)^0.75
    v1 <- exp(tv1 + eta.v1) * WT / 70
    q <- exp(tq + eta.q) * (WT / 70)^0.75
    v2 <- exp(tv2 + eta.v2) * WT / 70
    d/dt(central) <- -cl / v1 * central - q / v1 * central + q / v2 * peripheral
    d/dt(peripheral) <- q / v1 * central - q / v2 * peripheral
    cp <- 1000 * central / v1
    cp ~ add(add.sd) + prop(prop.sd)
  })
})
ui1 <- nlmixr2est::nlmixr2(function() {
  ini({
    tcl <- log(30); tv1 <- log(40); tq <- log(30); tv2 <- log(150)
    eta.cl ~ 0.09; eta.v1 ~ 0.09; eta.q ~ 0.09; eta.v2 ~ 0.09
    add.sd <- 1; prop.sd <- 0.2
  })
  model({
    cl <- exp(tcl + eta.cl) * (WT / 70)^0.75
    v1 <- exp(tv1 + eta.v1) * WT / 70
    q <- exp(tq + eta.q) * (WT / 70)^0.75
    v2 <- exp(tv2 + eta.v2) * WT / 70
    d/dt(central) <- -cl / v1 * central - q / v1 * central + q / v2 * peripheral
    d/dt(peripheral) <- q / v1 * central - q / v2 * peripheral
    cp <- 1000 * central / v1
    cp ~ add(add.sd) + prop(prop.sd)
  })
})
m <- nlmixr2data::mavoglurant
m <- m[!(m$EVID == 0 & m$MDV == 1), ]
ndata <- data.frame(ID = m$ID, TIME = m$TIME + 1000 * (m$OCC - 1), AMT = m$AMT, RATE = m$RATE,
                    EVID = ifelse(m$EVID == 1 & m$OCC == 2, 4L, m$EVID), DV = m$DV, CMT = "central", occ = m$OCC, WT = m$WT)
toIni <- function(s, between) c(list(tcl = s[["CL"]], tv1 = s[["V1"]], tq = s[["Q"]], tv2 = s[["V2"]],
                                     add.sd = exp(s[["sigma_abs"]]), prop.sd = exp(s[["sigma_rel"]])),
  unname(Map(function(eta, o) stats::as.formula(sprintf("%s ~ %.10g", eta, exp(2 * s[[o]]))),
             c("eta.cl", "eta.v1", "eta.q", "eta.v2", if (between) "iov.cl"),
             c("omega_CL_CL", "omega_V1_V1", "omega_Q_Q", "omega_V2_V2", if (between) "omega_CL_ID_OCC_CL_ID_OCC"))))
fromNlmixr2 <- function(f) {
  th <- f$theta
  # one matrix per level with occasions
  om <- if (is.list(f$omega)) unlist(lapply(unname(f$omega), diag)) else diag(f$omega)
  out <- c(CL = th[["tcl"]], V1 = th[["tv1"]], Q = th[["tq"]], V2 = th[["tv2"]], sigma_abs = log(th[["add.sd"]]),
           sigma_rel = log(th[["prop.sd"]]),
           omega_CL_CL = log(sqrt(om[["eta.cl"]])), omega_V1_V1 = log(sqrt(om[["eta.v1"]])),
           omega_Q_Q = log(sqrt(om[["eta.q"]])), omega_V2_V2 = log(sqrt(om[["eta.v2"]])))
  if (!is.na(om["iov.cl"])) out["omega_CL_ID_OCC_CL_ID_OCC"] <- log(sqrt(om[["iov.cl"]]))
  out
}
# rxode2::ini() evaluates the random effects in `envir`, where it needs lotri()
nfit <- function(s, est) {
  out <- NULL
  control <- if (est == "focei") nlmixr2est::foceiControl(print = 0) else nlmixr2est::saemControl(print = 0)
  model <- if (est == "focei") ui else ui1
  utils::capture.output(out <- tryCatch(suppressMessages(suppressWarnings(
    nlmixr2est::nlmixr2(do.call(rxode2::ini, c(list(model), toIni(s, est == "focei"), envir = asNamespace("rxode2"))), ndata,
                        est = est, control = control))), error = function(e) NULL))
  out
}
rxode2::setRxThreads(1L)
invisible(nfit(starts[1, ], "focei"))
invisible(nfit(starts[1, ], "saem"))
nlm <- rbind(
  do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
    rxode2::setRxThreads(1L)
    f <- nfit(starts[i, ], "focei")
    data.frame(package = "nlmixr2", method = "focei", start = i, value = if (is.null(f)) NA_real_ else focei(fromNlmixr2(f)))
  }, mc.cores = cores)),
  do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
    rxode2::setRxThreads(1L)
    set.seed(i)
    f <- nfit(starts1[i, ], "saem")
    data.frame(package = "nlmixr2", method = "saem", start = i, value = if (is.null(f)) NA_real_ else focei1(fromNlmixr2(f)))
  }, mc.cores = cores)))

## Comparison ----------------------------------------------------------------------------
runs <- rbind(dmod, nlm)
runs$fit <- paste(runs$package, runs$method)
best <- tapply(runs$value, runs$method, min, na.rm = TRUE)
runs <- runs[order(runs$fit, runs$value), ]
runs$rank <- stats::ave(seq_along(runs$value), runs$fit, FUN = seq_along)
ggplot2::ggplot(runs[is.finite(runs$value), ], ggplot2::aes(rank, value - best[method] + 1e-3, color = fit)) +
  ggplot2::geom_point() + ggplot2::scale_y_log10() + scale_color_dMod2() + theme_dMod2()
bestStart <- function(r, m) r$start[r$method == m][which.min(r$value[r$method == m])]
seconds <- c(
  dMod2focei = min(vapply(seq_len(reps), function(i)
    system.time(em(obj, random, starts[bestStart(dmod, "focei"), ], verbose = FALSE))[["elapsed"]], 0)),
  dMod2saem = min(vapply(seq_len(reps), function(i)
    system.time(em(obj1, random1, starts1[bestStart(dmod, "saem"), ], method = "saem", verbose = FALSE))[["elapsed"]], 0)),
  nlmixr2focei = min(vapply(seq_len(reps), function(i) system.time(nfit(starts[bestStart(nlm, "focei"), ], "focei"))[["elapsed"]], 0)),
  nlmixr2saem = min(vapply(seq_len(reps), function(i) system.time(nfit(starts1[bestStart(nlm, "saem"), ], "saem"))[["elapsed"]], 0)))
summary <- do.call(rbind, lapply(split(runs, runs$fit), function(r)
  data.frame(fit = r$fit[1], best = min(r$value, na.rm = TRUE), hits = sum(r$value - best[r$method[1]] < 0.1, na.rm = TRUE),
             failed = sum(is.na(r$value)))))
summary$seconds <- seconds[gsub(" ", "", summary$fit)]
summary
compileSeconds

## Single fits on 10 cores ---------------------------------------------------------------
# One FOCEI fit per package from five starts halfway out, the optimum moved by 0.3 on the log
# scale in every parameter with random signs, on `threads` cores: dMod2 over the conditions,
# nlmixr2 over the subjects. Per start the minimum seconds over `reps` runs, and the -2 log L of
# both by FOCEI of dMod2 above the optimum.
threads <- 10L
optimum <- em(obj, random, init, verbose = FALSE)
set.seed(2)
halfway <- t(vapply(1:5, function(i) optimum$argument[names(init)] + 0.3 * sample(c(-1, 1), length(init), TRUE), init))
rxode2::setRxThreads(threads)
multi <- do.call(rbind, lapply(seq_len(nrow(halfway)), function(i) {
  d <- replicate(reps, {
    t <- system.time(f <- em(obj, random, halfway[i, ], cores = threads, verbose = FALSE))[["elapsed"]]
    c(t, f$value)
  })
  n <- replicate(reps, {
    t <- system.time(f <- nfit(halfway[i, ], "focei"))[["elapsed"]]
    c(t, if (is.null(f)) NA_real_ else focei(fromNlmixr2(f)))
  })
  data.frame(start = i, dMod2 = min(d[1, ]), nlmixr2 = min(n[1, ]), dMod2Above = d[2, 1] - optimum$value,
             nlmixr2Above = n[2, 1] - optimum$value)
}))
rxode2::setRxThreads(1L)
multi
