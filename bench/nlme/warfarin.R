## Warfarin PK/PD, one oral compartment and a turnover of the prothrombin complex activity whose
## synthesis the concentration inhibits: FOCEI and SAEM of em() against nlmixr2 FOCEI and SAEM,
## each a multistart of `nstart` fits from the same starts, one core per fit. Per method the best
## -2 log L (FOCEI by dMod2 at the estimates of every fit), the starts within 0.1 of it, the
## failed fits, and the time of the fit from the best start, the minimum over `reps` runs.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "warfarinBench"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)
reps <- 3
nstart <- 100
cores <- max(1L, min(20L, parallel::detectCores() - 2L))
dm <- c("forward", "reverse", "forward-reverse")

## dMod2 ---------------------------------------------------------------------------------
data(warfarin)
data <- as.datalist(warfarin, splitBy = "condition", keepCovariates = c("dose", "wt", "age", "sex"))
net <- eqnlist() |>
  assignCompartment(Ag = "depot") |>
  assignCompartment(Cc = "central", volume = "V") |>
  assignCompartment(PCA = "effect") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("Ag", "Cc", "Ka * Ag", "absorption") |>
  addReaction("Cc", "", "Cl * Cc", "elimination", rateCompartment = "unit") |>
  addReaction("", "PCA", "kout * pca0 * IC50 / (IC50 + Cc)", "synthesis") |>
  addReaction("PCA", "", "kout * PCA", "degradation")
x <- Xs(odemodel(net, modelname = "warfBench", compile = FALSE, deriv2 = TRUE, derivMode = dm))
g <- Y(c(cp = "Cc", pca = "PCA"), x, modelname = "warfBench_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(cp = "sqrt(sigma_cp_abs^2 + (sigma_cp_rel * cp)^2)", pca = "sigma_pca"), g, attachInput = FALSE,
       modelname = "warfBench_e", compile = FALSE, deriv2 = TRUE, derivMode = dm)
table <- attr(data, "condition.grid")
table$ID <- rownames(table)
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "IC50", "kout", "pca0", "Ag", "Cc", "PCA", "sigma_cp_abs",
                                             "sigma_cp_rel", "sigma_pca")), table = table)
trafo <- insert(trafo, "Ag <- d", d = dose)
trafo <- insert(trafo, "Cc <- 0")
trafo <- insert(trafo, "PCA <- pca0")
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "IC50", "kout", "pca0", "sigma_cp_abs", "sigma_cp_rel",
                                           "sigma_pca"))
trafo <- insert(trafo, V + Cl + pca0 + kout ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "warfBench_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compileSeconds <- system.time(compile(x, g, e, p, output = "warfBench"))[["elapsed"]]
obj <- normL2(data, g * x * p, errmodel = e)
random <- randomEffects(p)
init <- c(Ka = 0, V = 2, Cl = -2, IC50 = 0, kout = -3, pca0 = 4.6, sigma_cp_abs = log(0.1), sigma_cp_rel = log(0.1),
          sigma_pca = log(5),
          omega_V_V = log(0.3), omega_Cl_Cl = log(0.3), omega_pca0_pca0 = log(0.3), omega_kout_kout = log(0.3))
focei <- function(pars) em(obj, random, pars, control = list(trustControl = list(iterlim = 0)), verbose = FALSE)$value
set.seed(1)
starts <- rbind(init, t(vapply(seq_len(nstart - 1), function(i) init + rnorm(length(init), 0, 0.5), init)))

dmod <- do.call(rbind, lapply(c("focei", "saem"), function(m) do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
  options(dMod2.cores = 1L)
  set.seed(i)
  f <- tryCatch(em(obj, random, starts[i, ], method = m, verbose = FALSE), error = function(e) NULL)
  data.frame(package = "dMod2", method = m, start = i,
             value = if (is.null(f)) NA_real_ else focei(f$argument[names(init)]))
}, mc.cores = cores))))

## nlmixr2 -------------------------------------------------------------------------------
# The same model on nlmixr2data::warfarin without the four concentrations of 0.
ui <- nlmixr2est::nlmixr2(function() {
  ini({
    tka <- 0; tv <- 2; tcl <- -2; tic50 <- 0; tkout <- -3; tpca0 <- 4.6
    eta.v ~ 0.09; eta.cl ~ 0.09; eta.pca0 ~ 0.09; eta.kout ~ 0.09
    add.cp <- 0.1; prop.cp <- 0.1; add.pca <- 5
  })
  model({
    ka <- exp(tka); v <- exp(tv + eta.v); cl <- exp(tcl + eta.cl); ic50 <- exp(tic50)
    kout <- exp(tkout + eta.kout); pca0 <- exp(tpca0 + eta.pca0)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    eff(0) <- pca0
    d/dt(eff) <- kout * pca0 * ic50 / (ic50 + center / v) - kout * eff
    cp <- center / v
    pca <- eff
    cp ~ add(add.cp) + prop(prop.cp)
    pca ~ add(add.pca)
  })
})
ndata <- nlmixr2data::warfarin
ndata <- ndata[!(ndata$evid == 0 & ndata$dvid == "cp" & ndata$dv == 0), ]
toIni <- function(s) c(list(tka = s[["Ka"]], tv = s[["V"]], tcl = s[["Cl"]], tic50 = s[["IC50"]], tkout = s[["kout"]],
                            tpca0 = s[["pca0"]], add.cp = exp(s[["sigma_cp_abs"]]), prop.cp = exp(s[["sigma_cp_rel"]]),
                            add.pca = exp(s[["sigma_pca"]])),
                       unname(Map(function(eta, o) stats::as.formula(sprintf("%s ~ %.10g", eta, exp(2 * s[[o]]))),
                                  c("eta.v", "eta.cl", "eta.pca0", "eta.kout"),
                                  c("omega_V_V", "omega_Cl_Cl", "omega_pca0_pca0", "omega_kout_kout"))))
fromNlmixr2 <- function(f) {
  th <- f$theta; om <- diag(f$omega)
  c(Ka = th[["tka"]], V = th[["tv"]], Cl = th[["tcl"]], IC50 = th[["tic50"]], kout = th[["tkout"]], pca0 = th[["tpca0"]],
    sigma_cp_abs = log(th[["add.cp"]]), sigma_cp_rel = log(th[["prop.cp"]]), sigma_pca = log(th[["add.pca"]]), omega_V_V = log(sqrt(om[["eta.v"]])),
    omega_Cl_Cl = log(sqrt(om[["eta.cl"]])), omega_pca0_pca0 = log(sqrt(om[["eta.pca0"]])),
    omega_kout_kout = log(sqrt(om[["eta.kout"]])))
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
invisible(nfit(starts[1, ], "focei"))
nlm <- do.call(rbind, lapply(c("focei", "saem"), function(m) do.call(rbind, parallel::mclapply(seq_len(nstart), function(i) {
  rxode2::setRxThreads(1L)
  set.seed(i)
  f <- nfit(starts[i, ], m)
  data.frame(package = "nlmixr2", method = m, start = i, value = if (is.null(f)) NA_real_ else focei(fromNlmixr2(f)))
}, mc.cores = cores))))

## Comparison ----------------------------------------------------------------------------
runs <- rbind(dmod, nlm)
runs$fit <- paste(runs$package, runs$method)
best <- min(runs$value, na.rm = TRUE)
runs <- runs[order(runs$fit, runs$value), ]
runs$rank <- stats::ave(seq_along(runs$value), runs$fit, FUN = seq_along)
ggplot2::ggplot(runs[is.finite(runs$value), ], ggplot2::aes(rank, value - best + 1e-3, color = fit)) +
  ggplot2::geom_point() + ggplot2::scale_y_log10() + scale_color_dMod() + theme_dMod()
seconds <- c(
  dMod2focei = min(vapply(seq_len(reps), function(i) system.time(em(obj, random, starts[which.min(dmod$value[dmod$method == "focei"]), ],
                                                                     verbose = FALSE))[["elapsed"]], 0)),
  dMod2saem = min(vapply(seq_len(reps), function(i) system.time(em(obj, random, starts[which.min(dmod$value[dmod$method == "saem"]), ],
                                                                    method = "saem", verbose = FALSE))[["elapsed"]], 0)),
  nlmixr2focei = min(vapply(seq_len(reps), function(i) system.time(nfit(starts[which.min(nlm$value[nlm$method == "focei"]), ],
                                                                        "focei"))[["elapsed"]], 0)),
  nlmixr2saem = min(vapply(seq_len(reps), function(i) system.time(nfit(starts[which.min(nlm$value[nlm$method == "saem"]), ],
                                                                       "saem"))[["elapsed"]], 0)))
summary <- do.call(rbind, lapply(split(runs, runs$fit), function(r)
  data.frame(fit = r$fit[1], best = min(r$value, na.rm = TRUE), hits = sum(r$value - best < 0.1, na.rm = TRUE),
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
