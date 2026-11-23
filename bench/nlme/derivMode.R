## FOCEI with the second derivatives of the volume term forward over forward and forward over
## reverse (`control$derivMode`): the time (minimum over `reps` runs), -2 log L, the iterations
## and the standard errors, on Theophylline (seven parameters), on STAT5 dimerization in ten cell
## lines with four parameters varying between them (Boehm et al. 2014), and on an alignment with
## a random scale per gel.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "derivModeBench"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)
reps <- 3
modes <- c("forward-forward", "forward-reverse")
dm <- c("forward", "reverse", "forward-forward", "forward-reverse")

## Theophylline ---------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
data <- as.datalist(data.frame(name = "y", time = Theoph$Time, value = Theoph$conc, sigma = NA_real_,
                               condition = Theoph$Subject))
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "", "Ka * Ag", "absorption")
reactions <- addReaction(reactions, "", "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "", "Cl/V * Cc", "elimination")
x   <- Xs(odemodel(reactions, modelname = "dmTheoph", compile = FALSE, derivMode = dm))
g   <- Y(c(y = "Cc"), x, modelname = "dmTheoph_obs", compile = FALSE, derivMode = dm)
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, compile = FALSE, derivMode = dm,
         modelname = "dmTheoph_err")
subjects <- sort(unique(Theoph$Subject))
covariates <- data.frame(ID = subjects, dose = tapply(Theoph$Dose * Theoph$Wt, Theoph$Subject, `[`, 1)[subjects],
                         row.names = subjects)
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma_add")), table = covariates)
trafo <- insert(trafo, "Ag <- dose", dose = dose)
trafo <- insert(trafo, "Cc <- 0")
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
trafo <- insert(trafo, Ka + V + Cl ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", compile = FALSE, modelname = "dmTheoph_p", derivMode = dm)
compile(x, g, err, p, output = "dmTheoph")
obj  <- normL2(data, g * x * p, errmodel = err)
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5))

theoph <- lapply(setNames(nm = modes), function(mode) {
  fit <- NULL
  secs <- min(vapply(seq_len(reps), function(i)
    system.time(fit <<- em(obj, randomEffects(p), init, control = list(derivMode = mode),
                           verbose = FALSE))[["elapsed"]], 0))
  list(fit = fit, seconds = secs)
})
data.frame(mode = modes, seconds = sapply(theoph, `[[`, "seconds"),
           value = sapply(theoph, function(r) r$fit$value), iterations = sapply(theoph, function(r) r$fit$iterations))
sapply(theoph, function(r) summary(r$fit)$population[, "se"])

## STAT5 in ten cell lines ----------------------------------------------------------------
# The reactions of demo("cellLines"), normal random effects on the STAT5 pool, the split of the pool and two
# transport rates over the lines.
epo <- "1.25e-7*exp(-Epo_degradation_BaF3*time)"
stat5 <- eqnlist() |>
  assignCompartment(nucpApA = "nuc", nucpApB = "nuc", nucpBpB = "nuc", volume = "0.45") |>
  addReaction("2*STAT5A", "pApA", paste0("k_phos*", epo, "*STAT5A^2"), "phosphorylation AA", compartment = "cyt") |>
  addReaction("STAT5A + STAT5B", "pApB", paste0("k_phos*", epo, "*STAT5A*STAT5B"), "phosphorylation AB",
              compartment = "cyt") |>
  addReaction("2*STAT5B", "pBpB", paste0("k_phos*", epo, "*STAT5B^2"), "phosphorylation BB", compartment = "cyt") |>
  addReaction("pApA", "nucpApA", "k_imp_homo*pApA", "import AA", compartment = "cyt") |>
  addReaction("pApB", "nucpApB", "k_imp_hetero*pApB", "import AB", compartment = "cyt") |>
  addReaction("pBpB", "nucpBpB", "k_imp_homo*pBpB", "import BB", compartment = "cyt") |>
  addReaction("nucpApA", "2*STAT5A", "k_exp_homo*nucpApA", "export AA", compartment = "nuc") |>
  addReaction("nucpApB", "STAT5A + STAT5B", "k_exp_hetero*nucpApB", "export AB", compartment = "nuc") |>
  addReaction("nucpBpB", "2*STAT5B", "k_exp_homo*nucpBpB", "export BB", compartment = "nuc") |>
  setCompartmentVolume(cyt = "1.4")
sx <- Xs(odemodel(stat5, modelname = "dmStat5", compile = FALSE, derivMode = dm),
         options = list(atol = 1e-8, rtol = 1e-6, maxattempts = 100L, maxsteps = 1e6))
obs5 <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))")
sg <- Y(obs5, sx, modelname = "dmStat5_obs", compile = FALSE, derivMode = dm)
serr <- Y(c(pSTAT5A_rel = "sd_pA", pSTAT5B_rel = "sd_pB", rSTAT5A_rel = "sd_rA"), sg, attachInput = FALSE,
          compile = FALSE, derivMode = dm, modelname = "dmStat5_err")
lines <- paste0("L", 1:10)
rates <- c("Epo_degradation_BaF3", "k_exp_hetero", "k_exp_homo", "k_imp_hetero", "k_imp_homo", "k_phos",
           "sd_pA", "sd_pB", "sd_rA", "STAT5tot")
str5 <- branch(define(NULL, "x <- x", x = getParameters(sx, sg, serr)), table = data.frame(line = lines, row.names = lines))
str5 <- define(str5, "x <- 0", x = c("pApA", "pApB", "pBpB", "nucpApA", "nucpApB", "nucpBpB"))
str5 <- insert(str5, "STAT5A <- STAT5tot*ratio")
str5 <- insert(str5, "STAT5B <- STAT5tot - STAT5tot*ratio")
str5 <- insert(str5, "specC17 <- 0.107")
str5 <- insert(str5, "x <- exp10(x)", x = rates)
str5 <- insert(str5, "ratio <- 1/(1 + exp(-ratio))")
str5 <- insert(str5, x ~ 1 + (1 | line), x = c("STAT5tot", "ratio", "k_exp_homo", "k_imp_hetero"))
sp <- P(str5, method = "explicit", compile = FALSE, modelname = "dmStat5_p", derivMode = dm)
compile(sx, sg, serr, sp, output = "dmStat5")
published <- c(log10(c(Epo_degradation_BaF3 = 0.02698, k_exp_hetero = 1.0007e-05, k_exp_homo = 0.006170,
                       k_imp_hetero = 0.01637, k_imp_homo = 97749, k_phos = 15767, sd_pA = 3.853,
                       sd_pB = 6.591, sd_rA = 3.153, STAT5tot = 207.6)), ratio = qlogis(0.693))
set.seed(2)
etas5 <- unlist(lapply(c("STAT5tot", "ratio", "k_exp_homo", "k_imp_hetero"), function(v)
  setNames(rnorm(10, 0, 0.2), paste0("eta_", v, "_", lines))))
times5 <- c(0, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 160, 200, 240)
sim5 <- (sg * sx * sp)(times5, c(published, etas5), deriv = FALSE)
noise <- 10^published[c("sd_pA", "sd_pB", "sd_rA")]
data5 <- as.datalist(do.call(rbind, lapply(names(sim5), function(cn) do.call(rbind, Map(function(o, s)
  data.frame(name = o, time = times5, value = sim5[[cn]][, o] + rnorm(length(times5), 0, s), sigma = NA_real_,
             condition = cn), names(obs5), noise)))))
obj5 <- normL2(data5, sg * sx * sp, errmodel = serr)
fixed5 <- published["STAT5tot"]
init5 <- c(round(published[setdiff(names(published), "STAT5tot")], 1))

stat5fits <- lapply(setNames(nm = modes), function(mode) {
  fit <- NULL
  secs <- min(vapply(seq_len(reps), function(i)
    system.time(fit <<- em(obj5, randomEffects(sp), init5, fixed = fixed5, control = list(derivMode = mode),
                           verbose = FALSE))[["elapsed"]], 0))
  list(fit = fit, seconds = secs)
})
data.frame(mode = modes, seconds = sapply(stat5fits, `[[`, "seconds"),
           value = sapply(stat5fits, function(r) r$fit$value), iterations = sapply(stat5fits, function(r) r$fit$iterations))
sapply(stat5fits, function(r) summary(r$fit)$population[, "se"])

## An alignment with a random scale per gel ------------------------------------------------
# 50 gels with two targets each, the scale per gel and target plus a normal random effect per gel.
set.seed(1)
targets <- sprintf("T%02d", 1:10); tps <- sprintf("t%02d", 1:8)
truth <- matrix(rlnorm(80, 0, 1), 10, dimnames = list(targets, tps))
gels <- sprintf("G%02d", 1:50)
pairs <- data.frame(gel = rep(gels, each = 2), target = targets[(rep(seq_along(gels), each = 2) * 2 + 0:1) %% 10 + 1])
pairs$scale <- rlnorm(nrow(pairs), 3, 0.3) * rlnorm(50, 0, 0.4)[match(pairs$gel, gels)]
ad <- pairs[rep(seq_len(nrow(pairs)), each = 16), ]
ad$time <- rep(rep(tps, each = 2), nrow(pairs))
ad$value <- truth[cbind(ad$target, ad$time)] * ad$scale * rlnorm(nrow(ad), 0, 0.1)
ad$experiment <- paste(ad$gel, ad$target, sep = "_")
am <- alignmentModel(ad, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ target * time,
                     scaling = sj ~ experiment + (1 | gel), error = sigmaR ~ 1,
                     trafo = c("yi <- exp(yi)", "sj <- exp(sj)", "sigmaR <- exp(sigmaR)"))
a0 <- align(am, fits = 1)
startA <- a0$fit$argument + 0.2
alignFits <- lapply(setNames(nm = modes), function(mode) {
  fit <- NULL
  secs <- min(vapply(seq_len(reps), function(i)
    system.time(fit <<- em(a0$objfn, a0$randomEffects, startA, control = list(derivMode = mode),
                           verbose = FALSE))[["elapsed"]], 0))
  list(fit = fit, seconds = secs)
})
data.frame(mode = modes, seconds = sapply(alignFits, `[[`, "seconds"),
           value = sapply(alignFits, function(r) r$fit$value), iterations = sapply(alignFits, function(r) r$fit$iterations))
