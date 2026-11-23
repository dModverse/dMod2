## Theophylline, one-compartment model with first-order absorption and log-normal Ka, V and Cl:
## the four methods of em() against nlmixr2 FOCEI and SAEM on the same data, the time (minimum
## over `reps` runs), -2 log L with all constants and the estimates on the log scale.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "theophBench"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)
reps <- 3

## dMod2 ---------------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
data <- as.datalist(data.frame(name = "y", time = Theoph$Time, value = Theoph$conc, sigma = NA_real_,
                               condition = Theoph$Subject))
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "", "Ka * Ag", "absorption")
reactions <- addReaction(reactions, "", "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "", "Cl/V * Cc", "elimination")
x   <- Xs(odemodel(reactions, modelname = "theophBench", compile = FALSE, derivMode = "forward-forward"))
g   <- Y(c(y = "Cc"), x, modelname = "theophBench_obs", compile = FALSE, derivMode = "forward-forward")
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, compile = FALSE, derivMode = "forward-forward",
         modelname = "theophBench_err")
subjects <- sort(unique(Theoph$Subject))
covariates <- data.frame(ID = subjects, dose = tapply(Theoph$Dose * Theoph$Wt, Theoph$Subject, `[`, 1)[subjects],
                         row.names = subjects)
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma_add")), table = covariates)
trafo <- insert(trafo, "Ag <- dose", dose = dose)
trafo <- insert(trafo, "Cc <- 0")
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
trafo <- insert(trafo, Ka + V + Cl ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", compile = FALSE, modelname = "theophBench_p", derivMode = "forward-forward")
compile(x, g, err, p, output = "theophBench")
obj  <- normL2(data, g * x * p, errmodel = err)
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5))

methods <- c("focei", "quadrature", "laplaceEM", "saem")
fits <- list(); seconds <- list()
for (m in methods) {
  times <- vapply(seq_len(reps), function(i) {
    set.seed(i)
    system.time(fits[[m]] <<- em(obj, randomEffects(p), init = init, method = m, verbose = FALSE))[["elapsed"]]
  }, 0)
  seconds[[m]] <- min(times)
}

## nlmixr2 -------------------------------------------------------------------------------
# The same model; theo_sd doses AMT = Dose * Wt into the depot, the observations are cp.
one.cmt <- function() {
  ini({
    tka <- 0; tv <- 3; tcl <- 0.7
    eta.ka ~ 0.09; eta.v ~ 0.09; eta.cl ~ 0.09
    add.sd <- 0.5
  })
  model({
    ka <- exp(tka + eta.ka)
    v  <- exp(tv + eta.v)
    cl <- exp(tcl + eta.cl)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}
theo <- nlmixr2data::theo_sd
nlm <- list()
for (est in c("focei", "saem")) {
  control <- if (est == "focei") nlmixr2est::foceiControl(print = 0) else nlmixr2est::saemControl(print = 0)
  times <- vapply(seq_len(reps), function(i) {
    set.seed(i)
    system.time(nlm[[est]] <<- nlmixr2est::nlmixr2(one.cmt, theo, est = est, control = control))[["elapsed"]]
  }, 0)
  seconds[[paste0("nlmixr2 ", est)]] <- min(times)
}

## Comparison ----------------------------------------------------------------------------
est <- c("Ka", "V", "Cl", "sigma_add", "omega_Ka_Ka", "omega_V_V", "omega_Cl_Cl")
fromNlmixr2 <- function(f) {
  th <- f$theta; om <- diag(f$omega)
  c(Ka = th[["tka"]], V = th[["tv"]], Cl = th[["tcl"]], sigma_add = log(th[["add.sd"]]),
    omega_Ka_Ka = log(sqrt(om[["eta.ka"]])), omega_V_V = log(sqrt(om[["eta.v"]])),
    omega_Cl_Cl = log(sqrt(om[["eta.cl"]])))
}
table <- rbind(t(vapply(fits, function(f) c(f$argument[est], minus2logL = f$value), numeric(8))),
               t(vapply(nlm, function(f) c(fromNlmixr2(f), minus2logL = -2 * as.numeric(stats::logLik(f))),
                        numeric(8))))
rownames(table) <- c(methods, paste("nlmixr2", names(nlm)))
round(cbind(table, seconds = unlist(seconds)[rownames(table)]), 4)
